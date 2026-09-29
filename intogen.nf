
// Set here a list of files or directories to use. E.g. Channel.fromPath(["/path/*", "/path2/file"], type: 'any')
INPUT = Channel.fromPath(params.input.tokenize())
OUTPUT = file(params.output)
STEPS_FOLDER = file(params.stepsFolder)
ANNOTATIONS = Channel.value(params.annotations)
REGIONS = Channel.value("${params.datasets}/regions/cds.regions.gz")

// Optional omics data: DNA methylation (e.g. array beta values) and RNA-seq.
// Files are assigned to cohorts by their name: <COHORT>.<anything> (e.g. TCGA_WXS_BRCA.beta.tsv.gz)
METHYLATION_INPUT = params.methylation ?
	Channel.fromPath(params.methylation.tokenize()).map{ it -> [it.baseName.split('\\.')[0], it] } : Channel.empty()
EXPRESSION_INPUT = params.expression ?
	Channel.fromPath(params.expression.tokenize()).map{ it -> [it.baseName.split('\\.')[0], it] } : Channel.empty()

// Empty files standing for missing optional inputs (e.g. cohorts without omics data)
PLACEHOLDERS = "${baseDir}/assets/placeholders"
OMICS_SAMPLES = Channel.value(file(params.omics_samples ?: "${PLACEHOLDERS}/NO_OMICS_SAMPLES"))

def placeholder(value, name) {
	value != null ? value : file("${PLACEHOLDERS}/${name}")
}

// Command line option for an optional input file (empty if it is a placeholder)
def optArg(flag, path) {
	path.name.startsWith('NO_') ? '' : "${flag} ${path}"
}

// Command line option for a list of optional input files (quoted, separated by spaces)
def optFiles(flag, files) {
	def names = (files instanceof java.nio.file.Path ? [files] : files.collect()).findAll { !it.name.startsWith('NO_') }
	names ? "${flag} \"${names.join(' ')}\"" : ''
}



process ParseInput {
	tag "Parse input ${input}"
	label "core"
	publishDir "${STEPS_FOLDER}/inputs", mode: "copy"
	errorStrategy 'finish'

	input:
		path input from INPUT
		path annotations from ANNOTATIONS

	output:
		path("*.parsed.tsv.gz") into COHORTS

	script:
		if ( input.toRealPath().toFile().isDirectory() || input.endsWith(".bginfo" )) {
			"""
			openvar groupby ${input} --header -g DATASET -q -s 'gzip > \${GROUP_KEY}.parsed.tsv.gz'
			"""
		}
		else {
			cohort = input.baseName.split('\\.')[0]
			"""
			openvar cat ${input} --header | gzip > ${cohort}.parsed.tsv.gz
			"""
		}
}


COHORTS
	.flatten()
	.map{it -> [it.baseName.split('\\.')[0], it]}
	.into{ COHORTS1; COHORTS2; COHORTS3; COHORTS4; COHORTS5 }

process LoadCancer {
	tag "Load cancer type ${cohort}"
	label "core"

	input:
		tuple val(cohort), path(input) from COHORTS1

	output:
		tuple val(cohort), stdout into CANCERS

	script:
		"""
		get_field.sh ${input} CANCER
		"""
}

CANCERS.into { CANCERS1; CANCERS2; CANCERS3 }


process LoadPlatform {
	tag "Load sequencing platform ${cohort}"
	label "core"

	input:
		tuple val(cohort), path(input) from COHORTS2

	output:
		tuple val(cohort), stdout into PLATFORMS

	script:
		"""
		get_field.sh ${input} PLATFORM
		"""
}

PLATFORMS.into { PLATFORMS1; PLATFORMS2; PLATFORMS3 }

process LoadGenome {
	tag "Load reference genome ${cohort}"
	label "core"

	input:
		tuple val(cohort), path(input) from COHORTS3

	output:
		tuple val(cohort), stdout into GENOMES

	script:
		"""
		get_field.sh ${input} GENOMEREF
		"""
}

CUTOFFS = ['WXS': 1000, 'WGS': 10000]

process ProcessVariants {
	tag "Process variants ${cohort}"
	label "core"
	errorStrategy 'ignore'  // if a cohort does not pass the filters, do not proceed with it
	publishDir "${STEPS_FOLDER}/variants", mode: "copy"

	input:
		tuple val(cohort), path(input), val(platform), val(genome) from COHORTS4.join(PLATFORMS1).join(GENOMES)

	output:
		tuple val(cohort), path(output) into VARIANTS
		tuple val(cohort), path("${output}.stats.json") into STATS_VARIANTS

	script:
		cutoff = CUTOFFS[platform]
		output = "${cohort}.tsv.gz"
		if (cutoff)
			"""
			parse-variants --input ${input} --output ${output} \
				--genome ${genome.toLowerCase()} \
				--cutoff ${cutoff}
			"""
		else
			error "Invalid cutoff. Check platform: $platform"

}

VARIANTS.into { VARIANTS1; VARIANTS2; VARIANTS3; VARIANTS4; VARIANTS5 }


process FormatSignature {
	tag "Prepare for signatures ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/signature", mode: "copy"

	input:
		tuple val(cohort), path(input) from VARIANTS1

	output:
		tuple val(cohort), path(output) into VARIANTS_SIG

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format signature
		"""

}

REGIONS_PREFIX = ['WXS': 'cds', 'WGS': 'wg']

process ComputeProfile {
	tag "ComputeProfile ${cohort}"
	label "bgsignature"
	publishDir "${STEPS_FOLDER}/signature", mode: "copy"

	input:
		tuple val(cohort), path(input), val(platform) from VARIANTS_SIG.join(PLATFORMS2)

	output:
		tuple val(cohort), path(output) into SIGNATURES

	script:
		prefix = REGIONS_PREFIX[platform]
		output = "${cohort}.sig.json"
		if (prefix)
			"""
			bgsignature normalize -m ${input} \
				-r ${params.datasets}/regions/${prefix}.regions.gz \
				--normalize ${params.datasets}/signature/${prefix}.counts.gz \
				-s 3 -g hg38 --collapse \
				--cores ${task.cpus} \
				-o ${output}
			"""
		else
			error "Invalid prefix. Check platform: $platform"

}

SIGNATURES.into{ SIGNATURES1; SIGNATURES2; SIGNATURES3; SIGNATURES4; SIGNATURES5 }


process FormatFML {
	tag "Prepare for FML ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/oncodrivefml", mode: "copy"

	input:
		tuple val(cohort), path(input) from VARIANTS2

	output:
		tuple val(cohort), path(output) into VARIANTS_FML

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format fml
		"""

}

process OncodriveFML {
    tag "OncodriveFML ${cohort}"
    publishDir "${STEPS_FOLDER}/oncodrivefml", mode: "copy"

    input:
        tuple val(cohort), path(input), path(signature)  from VARIANTS_FML.join(SIGNATURES1)
        path regions from REGIONS

    output:
        tuple val(cohort), path("out/*.tsv.gz") into OUT_ONCODRIVEFML

	script:
		seedOpt = (params.seed == null)? '': "--seed ${params.seed}"
		debugOpt =  (params.debug)? '--debug': ''
		"""
		oncodrivefml -i ${input} -e ${regions} --signature ${signature} \
			-c /oncodrivefml/oncodrivefml_v2.conf  --cores ${task.cpus} \
			-o out ${seedOpt} ${debugOpt}
		"""
}


process FormatCLUSTL {
	tag "Prepare for CLUSTL ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/oncodriveclustl", mode: "copy"

	input:
		tuple val(cohort), path(input) from VARIANTS3

	output:
		tuple val(cohort), path(output) into VARIANTS_CLUSTL

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format clustl
		"""

}

process OncodriveCLUSTL {
    tag "OncodriveCLUSTL ${cohort}"
    publishDir "${STEPS_FOLDER}/oncodriveclustl", mode: "copy"

    input:
        tuple val(cohort), path(input), path(signature), val(cancer) from VARIANTS_CLUSTL.join(SIGNATURES2).join(CANCERS1)
        path regions from REGIONS

    output:
        tuple val(cohort), path("${cohort}.elements_results.txt") into OUT_ONCODRIVECLUSTL
        tuple val(cohort), path("${cohort}.clusters_results.tsv") into OUT_ONCODRIVECLUSTL_CLUSTERS

	script:
		seedOpt = (params.seed == null)? '': "--seed ${params.seed}"
		debugOpt =  (params.debug)? '--log-level debug': ''
		if (['CM', 'SBCC', 'SSCC'].contains(cancer))
			"""
			oncodriveclustl -i ${input} -r ${regions} \
				-g hg38 -sim region_restricted -n 1000 -kmer 5 \
				-sigcalc region_normalized \
				--concatenate \
				-c ${task.cpus} \
				-o ${cohort} ${seedOpt} ${debugOpt}
			
			mv ${cohort}/elements_results.txt ${cohort}.elements_results.txt
			mv ${cohort}/clusters_results.tsv ${cohort}.clusters_results.tsv
			"""
		else
			"""
			oncodriveclustl -i ${input} -r ${regions} \
				-g hg38 -sim region_restricted -n 1000 -kmer 3 \
				-sig ${signature} --concatenate \
				-c ${task.cpus} \
				-o ${cohort} ${seedOpt} ${debugOpt}

			mv ${cohort}/elements_results.txt ${cohort}.elements_results.txt
			mv ${cohort}/clusters_results.tsv ${cohort}.clusters_results.tsv	
			"""
}


process FormatDNDSCV {
	tag "Prepare for DNDSCV ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/dndscv", mode: "copy"

	input:
		tuple val(cohort), path(input) from VARIANTS4

	output:
		tuple val(cohort), path(output) into VARIANTS_DNDSCV

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format dndscv
		"""
}

process dNdScv {
    tag "dNdScv ${cohort}"
    publishDir "${STEPS_FOLDER}/dndscv", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_DNDSCV

    output:
        tuple val(cohort), path("${cohort}.dndscv.tsv.gz") into OUT_DNDSCV
        tuple val(cohort), path("${cohort}.dndscv_annotmuts.tsv.gz") into OUT_DNDSCV_ANNOTMUTS
        tuple val(cohort), path("${cohort}.dndscv_genemuts.tsv.gz") into OUT_DNDSCV_GENEMUTS

	script:
		"""
		Rscript /dndscv/dndscv.R \
			${input} ${cohort}.dndscv.tsv.gz \
			${cohort}.dndscv_annotmuts.tsv.gz \
			${cohort}.dndscv_genemuts.tsv.gz
		"""
}

OUT_DNDSCV.into{ OUT_DNDSCV1; OUT_DNDSCV2 }

process FormatVEP {
	tag "Prepare for VEP ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/vep", mode: "copy"

	input:
		tuple val(cohort), path(input) from VARIANTS5

	output:
		tuple val(cohort), path(output) into VARIANTS_VEP

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output}.tmp \
			--format vep
		sort -k1V -k2n -k3n ${output}.tmp |gzip > ${output}
		"""

}

process VEP {
	tag "VEP ${cohort}"
	publishDir "${STEPS_FOLDER}/vep", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_VEP

    output:
        tuple val(cohort), path(output) into OUT_VEP

	script:
		output = "${cohort}.vep.tsv.gz"
		"""
		vep -i ${input} -o STDOUT --assembly GRCh38 \
			--no_stats --cache --offline --symbol \
			--protein --tab --canonical --mane \
			--dir ${params.datasets}/vep \
			| grep -v ^## | gzip > ${output}
		"""
}


process ProcessVEPoutput {
	tag "Process vep output ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/vep", mode: "copy"

    input:
        tuple val(cohort), path(input) from OUT_VEP

    output:
        tuple val(cohort), path(output) into PARSED_VEP
        tuple val(cohort), path("${output}.stats.json") into STATS_VEP

	script:
		output = "${cohort}.tsv.gz"
		"""
		parse-vep --input ${input} --output ${output}
		"""
}


PARSED_VEP.into { PARSED_VEP1; PARSED_VEP2; PARSED_VEP3; PARSED_VEP4; PARSED_VEP5; PARSED_VEP6; PARSED_VEP7; PARSED_VEP8; PARSED_VEP9 }

process FilterNonSynonymous {
	tag "Filter non synonymus ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/nonsynonymous", mode: "copy"

    input:
        tuple val(cohort), path(input) from PARSED_VEP1

    output:
        tuple val(cohort), path(output) into PARSED_VEP_NONSYNONYMOUS

	script:
		output = "${cohort}.vep_nonsyn.tsv.gz"
		"""
		parse-nonsynonymous --input ${input} --output ${output}
		"""
}


process FormatSMRegions {
	tag "Prepare for SMRegions ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/smregions", mode: "copy"

	input:
		tuple val(cohort), path(input) from PARSED_VEP_NONSYNONYMOUS

	output:
		tuple val(cohort), path(output) into VARIANTS_SMREGIONS

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format smregions
		"""
}

process SMRegions {
	tag "SMRegions ${cohort}"
	publishDir "${STEPS_FOLDER}/smregions", mode: "copy"

    input:
        tuple val(cohort), path(input), path(signature)  from VARIANTS_SMREGIONS.join(SIGNATURES3)
        path regions from REGIONS

    output:
        tuple val(cohort), path(output) into OUT_SMREGIONS

	script:
		output = "${cohort}.smregions.tsv.gz"
		seedOpt = (params.seed == null)? '': "--seed ${params.seed}"
		debugOpt =  (params.debug)? '--debug': ''
		"""
		smregions -m ${input} -e ${regions} \
			-r ${params.datasets}/smregions/regions_pfam.tsv \
			-s ${signature} --cores ${task.cpus} \
			-c /smregions/smregions.conf \
			-o ${output} ${seedOpt} ${debugOpt}
		"""
}

OUT_SMREGIONS.into { OUT_SMREGIONS1; OUT_SMREGIONS2 }


process FormatCBaSE {
	tag "Prepare for CBaSE ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/cbase", mode: "copy"

	input:
		tuple val(cohort), path(input) from PARSED_VEP2

	output:
		tuple val(cohort), path(output) into VARIANTS_CBASE

	script:
		output = "${cohort}.in.tsv"
		"""
		format-variants --input ${input} --output ${output} \
			--format cbase
		"""
}

process CBaSE {
	tag "CBaSE ${cohort}"
	publishDir "${STEPS_FOLDER}/cbase", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_CBASE

    output:
        tuple val(cohort), path(output) into OUT_CBASE

	script:
		output = "${cohort}.cbase.tsv.gz"
		"""
		mkdir -p Output/

		python /cbase/cbase.py ${input} ${params.datasets}/cbase 0 output
		tail -n+2 Output/q_values_output.txt | gzip > ${output}
		"""
}

process FormatMutPanning {
	tag "Prepare for MutPanning ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/mutpanning", mode: "copy"

	input:
		tuple val(cohort), path(input) from PARSED_VEP3

	output:
		tuple val(cohort), path(muts), path(samples) into VARIANTS_MUTPANNING

	script:
		muts = "${cohort}.in_muts.tsv"
		samples = "${cohort}.in_samples.tsv"
		"""
		format-variants --input ${input} --output ${muts} \
			--format mutpanning-mutations
		format-variants --input ${input} --output ${samples} \
			--format mutpanning-samples
		"""
}

process MutPanning {
	tag "MutPanning ${cohort}"
	publishDir "${STEPS_FOLDER}/mutpanning", mode: "copy"

    input:
        tuple val(cohort), path(mutations), path(samples) from VARIANTS_MUTPANNING

    output:
        tuple val(cohort), path("out/SignificanceFiltered/Significance${cohort}.txt") into OUT_MUTPANNING

	script:
		// TODO remove the creation of the out file or move to the container
		"""
		mkdir -p out/SignificanceFiltered
		echo "Name\tTargetSize\tTargetSizeSyn\tCount\tCountSyn\tSignificance\tFDR\n" \
			> out/SignificanceFiltered/Significance${cohort}.txt
		java -cp /mutpanning/MutPanning.jar MutPanning \
			out ${mutations} ${samples} ${params.datasets}/mutpanning/Hg19/
		"""
}


process FormatHotMAPS {
	tag "Prepare for HotMAPS ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/hotmaps", mode: "copy"

	input:
		tuple val(cohort), path(input) from PARSED_VEP4

	output:
		tuple val(cohort), path(output) into VARIANTS_HOTMAPS

	script:
		output = "${cohort}.in.maf"
		"""
		format-variants --input ${input} --output ${output} \
			--format hotmaps
		"""
}

process HotMAPS {
	tag "HotMAPS ${cohort}"
	publishDir "${STEPS_FOLDER}/hotmaps", mode: "copy"

    input:
        tuple val(cohort), path(input), path(signatures) from VARIANTS_HOTMAPS.join(SIGNATURES4)

    output:
        tuple val(cohort), path("*.out.gz") into OUT_HOTMAPS
        tuple val(cohort), path("*.clusters.gz") into OUT_HOTMAPS_CLUSTERS

	script:
		"""
		/bin/sh /hotmaps/hotmaps.sh ${input} . ${signatures} \
			${params.datasets}/hotmaps ${task.cpus}
		"""
}


/*
 * Omics layers (optional): DNA methylation and RNA-seq
 */

process ParseMethylation {
	tag "Parse methylation ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/methylation", mode: "copy"

	input:
		tuple val(cohort), path(input) from METHYLATION_INPUT
		path samples from OMICS_SAMPLES

	output:
		tuple val(cohort), path(output), path(output_normal) into METHYLATION_MATRIX
		path("${output}.stats.json") into STATS_METHYLATION_MATRIX

	script:
		output = "${cohort}.promoter_methylation.tsv.gz"
		output_normal = "${cohort}.promoter_methylation.normal.tsv.gz"
		probes = params.methylation_probes ? file(params.methylation_probes) : "${params.datasets}/methylation/promoter_probes.tsv.gz"
		"""
		parse-methylation --input ${input} \
			--output ${output} --output-normal ${output_normal} \
			--probes ${probes} \
			--values ${params.methylation_values} \
			${optArg('--samples', samples)}
		"""
}

METHYLATION_MATRIX.into { METHYLATION_MATRIX1; METHYLATION_MATRIX2; METHYLATION_MATRIX3 }

process ParseExpression {
	tag "Parse expression ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/expression", mode: "copy"

	input:
		tuple val(cohort), path(input) from EXPRESSION_INPUT
		path samples from OMICS_SAMPLES

	output:
		tuple val(cohort), path(output) into EXPRESSION_MATRIX
		path("${output}.stats.json") into STATS_EXPRESSION_MATRIX

	script:
		output = "${cohort}.expression_matrix.tsv.gz"
		"""
		parse-expression --input ${input} --output ${output} \
			--units ${params.expression_units} \
			${optArg('--samples', samples)}
		"""
}

EXPRESSION_MATRIX.into { EXPRESSION_MATRIX1; EXPRESSION_MATRIX2; EXPRESSION_MATRIX3; EXPRESSION_MATRIX4 }

process ExpressionAnalysis {
	tag "Expression analysis ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/expression", mode: "copy"

	input:
		tuple val(cohort), path(input) from EXPRESSION_MATRIX1

	output:
		tuple val(cohort), path(output) into OUT_EXPRESSION
		tuple val(cohort), path(events) into OUT_EXPRESSION_EVENTS
		path("${output}.stats.json") into STATS_EXPRESSION

	script:
		output = "${cohort}.expression.tsv.gz"
		events = "${cohort}.expression_events.tsv.gz"
		"""
		expression-analysis --input ${input} --output ${output} --events ${events}
		"""
}

OUT_EXPRESSION.into { OUT_EXPRESSION1; OUT_EXPRESSION2; OUT_EXPRESSION3 }

// Methylation of every cohort, with its expression when available
METHYLATION_ANALYSIS_INPUTS = METHYLATION_MATRIX1
	.join(EXPRESSION_MATRIX2, remainder: true)
	.filter { it[1] != null }
	.map { it -> [it[0], it[1], it[2], placeholder(it[3], 'NO_EXPRESSION_MATRIX')] }

process MethylationAnalysis {
	tag "Methylation analysis ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/methylation", mode: "copy"

	input:
		tuple val(cohort), path(input), path(normal), path(expression) from METHYLATION_ANALYSIS_INPUTS

	output:
		tuple val(cohort), path(output) into OUT_METHYLATION
		tuple val(cohort), path(events) into OUT_METHYLATION_EVENTS
		path("${output}.stats.json") into STATS_METHYLATION

	script:
		output = "${cohort}.methylation.tsv.gz"
		events = "${cohort}.methylation_events.tsv.gz"
		"""
		methylation-analysis --input ${input} --normal ${normal} \
			${optArg('--expression', expression)} \
			--output ${output} --events ${events}
		"""
}

OUT_METHYLATION.into { OUT_METHYLATION1; OUT_METHYLATION2; OUT_METHYLATION3 }
OUT_METHYLATION_EVENTS.into { OUT_METHYLATION_EVENTS1; OUT_METHYLATION_EVENTS2 }

// Cohorts with mutations and at least one omics layer
OMICS_FEATURES_INPUTS = PARSED_VEP8
	.join(EXPRESSION_MATRIX3, remainder: true)
	.join(OUT_EXPRESSION1, remainder: true)
	.join(METHYLATION_MATRIX2.map { it -> [it[0], it[1]] }, remainder: true)
	.join(OUT_METHYLATION1, remainder: true)
	.join(OUT_METHYLATION_EVENTS1, remainder: true)
	.filter {
		if (it[1] == null) {
			log.warn "Omics data of cohort ${it[0]} are not integrated: there are no mutation results for it " +
				"(omics file names must start with the cohort ID followed by a dot)"
		}
		it[1] != null && it[2..-1].any { x -> x != null }
	}
	.map { it -> [it[0], it[1],
		placeholder(it[2], 'NO_EXPRESSION_MATRIX'), placeholder(it[3], 'NO_EXPRESSION'),
		placeholder(it[4], 'NO_METHYLATION_MATRIX'), placeholder(it[5], 'NO_METHYLATION'),
		placeholder(it[6], 'NO_METHYLATION_EVENTS')] }

process OmicsFeatures {
	tag "Omics features ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/omics", mode: "copy"

	input:
		tuple val(cohort), path(mutations), path(expression_matrix), path(expression), path(methylation_matrix), path(methylation), path(methylation_events) from OMICS_FEATURES_INPUTS

	output:
		tuple val(cohort), path(output) into OUT_OMICS
		path("${output}.stats.json") into STATS_OMICS

	script:
		output = "${cohort}.omics.tsv.gz"
		"""
		omics-features --mutations ${mutations} --output ${output} \
			${optArg('--expression-matrix', expression_matrix)} \
			${optArg('--expression', expression)} \
			${optArg('--methylation-matrix', methylation_matrix)} \
			${optArg('--methylation', methylation)} \
			${optArg('--methylation-events', methylation_events)}
		"""
}

OUT_OMICS.into { OUT_OMICS1; OUT_OMICS2 }
OMICS_LIST = OUT_OMICS2.map { it -> it[1] }

process OmicsSummary {
	tag "Omics summary"
	label "core"
	publishDir "${OUTPUT}", mode: "copy"

	input:
		path(input) from OMICS_LIST.collect()

	output:
		path("omics.tsv") into OMICS_SUMMARY

	script:
		"""
		omics-summary --output omics.tsv ${input}
		"""
}


COMBINATION_INPUTS = OUT_ONCODRIVEFML.join(OUT_ONCODRIVECLUSTL).join(OUT_DNDSCV1).join(OUT_SMREGIONS1).join(OUT_CBASE).join(OUT_MUTPANNING).join(OUT_HOTMAPS)
if (params.integrate_omics) {
	// omics-based evidence is added to the combination of the cohorts that have it
	COMBINATION_INPUTS = COMBINATION_INPUTS
		.join(OUT_METHYLATION2, remainder: true)
		.join(OUT_EXPRESSION2, remainder: true)
		.filter { it[1] != null }
		.map { it -> it[0..-3] + [placeholder(it[-2], 'NO_METHYLATION'), placeholder(it[-1], 'NO_EXPRESSION')] }
} else {
	COMBINATION_INPUTS = COMBINATION_INPUTS
		.map { it -> it + [placeholder(null, 'NO_METHYLATION'), placeholder(null, 'NO_EXPRESSION')] }
}

process Combination {
	tag "Combination ${cohort}"
	publishDir "${STEPS_FOLDER}/combination", mode: "copy"

    input:
        tuple val(cohort), path(fml), path(clustl), path(dndscv), path(smregions), path(cbase), path(mutpanning), path(hotmaps), path(methylation), path(expression) from COMBINATION_INPUTS

    output:
        tuple val(cohort), path("${cohort}.05.out.gz") into OUT_COMBINATION

	script:
		"""
		intogen-combine -o ${cohort} \
			--oncodrivefml ${fml} \
			--oncodriveclustl ${clustl} \
			--dndscv ${dndscv} \
			--smregions ${smregions} \
			--cbase ${cbase} \
			--mutpanning ${mutpanning} \
			--hotmaps ${hotmaps} \
			${optArg('--methylation', methylation)} \
			${optArg('--expression', expression)}
		"""

}


process FormatdeconstructSigs {
	tag "Prepare for deconstructSigs ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/deconstructSigs", mode: "copy"

	input:
		tuple val(cohort), path(input) from PARSED_VEP5

	output:
		tuple val(cohort), path(output) into VARIANTS_DECONSTRUCTSIGS

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format deconstructsigs
		"""
}

VARIANTS_DECONSTRUCTSIGS.into{ VARIANTS_DECONSTRUCTSIGS1; VARIANTS_DECONSTRUCTSIGS2 }

process deconstructSigs {
	tag "deconstructSigs ${cohort}"
	publishDir "${STEPS_FOLDER}/deconstructSigs", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_DECONSTRUCTSIGS1

    output:
        tuple val(cohort), path(output) into OUT_DECONSTRUCTSIGS
        tuple val(cohort), path("*.signature_likelihood") into OUT_DECONSTRUCTSIGS_SIGLIKELIHOOD

	script:
		output = "${cohort}.deconstructsigs.tsv.gz"
		likelihood = "${cohort}.signature_likelihood"
		"""
		python3 /deconstructsig/run_deconstruct.py \
			--input-file ${input} --weights ${output} \
			--build hg38
		python3 /deconstructsig/signature_assignment.py \
			--input-file ${output} \
			--output-file ${likelihood}
		"""
}

// TODO add a process to combine all stats

process CohortCounts {
	tag "Count variants ${cohort}"

    input:
        tuple val(cohort), path(input), val(cancer), val(platform) from COHORTS5.join(CANCERS2).join(PLATFORMS3)

    output:
		tuple val(cohort), path("*.counts") into COHORT_COUNTS

	script:
		"""
		variants=`zcat ${input} |tail -n+2 |wc -l`
		samples=`zcat ${input} |tail -n+2 |cut -f1 |uniq |sort -u| wc -l`
		echo "${cohort}\t${cancer}\t${platform}\t\$variants\t\$samples" > ${cohort}.counts
		"""
}

COHORT_COUNTS.into{ COHORT_COUNTS1; COHORT_COUNTS2 }
COHORT_COUNTS_LIST = COHORT_COUNTS1.map{ it -> it[1] }

process CohortSummary {
	tag "Count variants"
	publishDir "${OUTPUT}", mode: "copy"

    input:
        path(input) from COHORT_COUNTS_LIST.collect()

    output:
		path(output) into COHORT_SUMMARY

	script:
		output="cohorts.tsv"
		"""
		echo 'COHORT\tCANCER_TYPE\tPLATFORM\tMUTATIONS\tSAMPLES' > ${output}
		cat ${input} >> ${output}
		"""
}


MUTATIONS_INPUTS = PARSED_VEP6.map { it -> it[1] }

process MutationsSummary {
	tag "Mutations"
	publishDir "${OUTPUT}", mode: "copy"
	label "core"

    input:
        path(input) from MUTATIONS_INPUTS.collect()

    output:
		path(output) into MUTATIONS_SUMMARY

	script:
		output="mutations.tsv"
		"""
		mutations-summary --output ${output} \
			${input}
		"""
}


DRIVER_DISCOVERY_INPUTS = OUT_COMBINATION.join(VARIANTS_DECONSTRUCTSIGS2).join(OUT_DECONSTRUCTSIGS_SIGLIKELIHOOD).join(OUT_SMREGIONS2).join(OUT_ONCODRIVECLUSTL_CLUSTERS).join(OUT_HOTMAPS_CLUSTERS).join(OUT_DNDSCV2).join(CANCERS3)
if (params.methylation || params.expression) {
	// omics features are added to the cohorts that have them
	DRIVER_DISCOVERY_INPUTS = DRIVER_DISCOVERY_INPUTS
		.join(OUT_OMICS1, remainder: true)
		.filter { it[1] != null }
		.map { it -> it[0..-2] + [placeholder(it[-1], 'NO_OMICS')] }
} else {
	DRIVER_DISCOVERY_INPUTS = DRIVER_DISCOVERY_INPUTS.map { it -> it + [placeholder(null, 'NO_OMICS')] }
}

process DriverDiscovery {
	tag "Driver discovery ${cohort}"
	publishDir "${STEPS_FOLDER}/drivers", mode: "copy"
	label "core"

    input:
        tuple val(cohort), path(combination), path(deconstruct_in), path(sig_likelihood), path(smregions), path(clustl_clusters), path(hotmaps_clusters), path(dndscv), val(cancer), path(omics) from DRIVER_DISCOVERY_INPUTS

    output:
		path(output_drivers) into DRIVERS
		path(output_vet) into VET
		tuple val(cohort), path(output_vet) into VET_COHORT

	script:
		output_drivers = "${cohort}.drivers.tsv"
		output_vet = "${cohort}.vet.tsv"
		"""
		drivers-discovery --output_drivers ${output_drivers} \
			--output_vet ${output_vet} \
			--combination ${combination} \
			--mutations ${deconstruct_in} \
			--sig_likelihood ${sig_likelihood} \
			--smregions ${smregions} \
			--clustl_clusters ${clustl_clusters} \
			--hotmaps ${hotmaps_clusters} \
			--dndscv ${dndscv} \
			--ctype ${cancer} \
			--cohort ${cohort} \
			${optArg('--omics', omics)}
		"""
}


/*
 * Pathway analysis (optional): selection of gene sets through mutations and omics layers,
 * beyond their individually significant genes, and co-occurrence of dysregulation events
 */

if (params.pathways) {
	PATHWAY_INPUTS = OUT_DNDSCV_GENEMUTS.join(PARSED_VEP9).join(VET_COHORT)
		.join(OUT_METHYLATION3, remainder: true)
		.join(OUT_METHYLATION_EVENTS2, remainder: true)
		.join(METHYLATION_MATRIX3.map { it -> [it[0], it[1]] }, remainder: true)
		.join(OUT_EXPRESSION3, remainder: true)
		.join(OUT_EXPRESSION_EVENTS, remainder: true)
		.join(EXPRESSION_MATRIX4, remainder: true)
		.filter { it[1] != null }
		.map { it -> it[0..3] + [
			placeholder(it[4], 'NO_METHYLATION'), placeholder(it[5], 'NO_METHYLATION_EVENTS'),
			placeholder(it[6], 'NO_METHYLATION_MATRIX'), placeholder(it[7], 'NO_EXPRESSION'),
			placeholder(it[8], 'NO_EXPRESSION_EVENTS'), placeholder(it[9], 'NO_EXPRESSION_MATRIX')] }
} else {
	PATHWAY_INPUTS = Channel.empty()
}

process PathwayAnalysis {
	tag "Pathway analysis ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/pathways", mode: "copy"

	input:
		tuple val(cohort), path(genemuts), path(mutations), path(vet), path(methylation), path(methylation_events), path(methylation_matrix), path(expression), path(expression_events), path(expression_matrix) from PATHWAY_INPUTS

	output:
		tuple val(cohort), path(output) into OUT_PATHWAYS
		tuple val(cohort), path(cooccurrence) into OUT_PATHWAY_COOCCURRENCE
		tuple val(cohort), path(genes) into OUT_PATHWAY_GENES
		path("${output}.stats.json") into STATS_PATHWAYS

	script:
		output = "${cohort}.pathways.tsv.gz"
		cooccurrence = "${cohort}.pathway_cooccurrence.tsv.gz"
		genes = "${cohort}.pathway_genes.tsv.gz"
		geneSets = params.gene_sets ? file(params.gene_sets) : "${params.datasets}/pathways/gene_sets.tsv.gz"
		"""
		pathway-analysis --genemuts ${genemuts} --mutations ${mutations} --vet ${vet} \
			--gene-sets ${geneSets} \
			--output ${output} --cooccurrence ${cooccurrence} --genes ${genes} \
			${optArg('--methylation', methylation)} \
			${optArg('--methylation-events', methylation_events)} \
			${optArg('--methylation-matrix', methylation_matrix)} \
			${optArg('--expression', expression)} \
			${optArg('--expression-events', expression_events)} \
			${optArg('--expression-matrix', expression_matrix)}
		"""
}

PATHWAY_FILES = OUT_PATHWAYS.map { it -> it[1] }.mix(OUT_PATHWAY_COOCCURRENCE.map { it -> it[1] })

process PathwaySummary {
	tag "Pathway summary"
	label "core"
	publishDir "${OUTPUT}", mode: "copy"

	input:
		path(input) from PATHWAY_FILES.collect()

	output:
		path("pathways.tsv") into PATHWAYS_SUMMARY
		path("pathway_cooccurrence.tsv") into PATHWAY_COOCCURRENCE_SUMMARY
		path("pathway_modules.tsv") into PATHWAY_MODULES_SUMMARY

	script:
		"""
		pathway-summary --pathways pathways.tsv \
			--cooccurrence pathway_cooccurrence.tsv \
			--modules pathway_modules.tsv \
			${input}
		"""
}

// Significant gene sets of the genes of each cohort, added to the drivers
PATHWAY_GENES = params.pathways ?
	OUT_PATHWAY_GENES.map { it -> it[1] }.collect().ifEmpty(file("${PLACEHOLDERS}/NO_PATHWAYS")) :
	Channel.value(file("${PLACEHOLDERS}/NO_PATHWAYS"))

process DriverSummary {
	tag "Driver summary"
	publishDir "${OUTPUT}", mode: "copy"
	label "core"

    input:
        path (input) from DRIVERS.collect()
        path (input_vet) from VET.collect()
        path (mutations) from MUTATIONS_SUMMARY
        path (cohortsSummary) from COHORT_SUMMARY
        path (pathwayGenes) from PATHWAY_GENES

    output:
		path("drivers.tsv") into DRIVERS_SUMMARY
		path("unique_drivers.tsv") into UNIQUE_DRIVERS
		path("unfiltered_drivers.tsv") into UNFILTER_DRIVERS

	script:
		"""
		drivers-summary \
			--mutations ${mutations} \
			--cohorts ${cohortsSummary} \
			${optFiles('--pathways', pathwayGenes)} \
			${input} "${input_vet}"
		"""
}


process ParseProfile {
	tag "Parsing profile ${cohort}"
	publishDir "${STEPS_FOLDER}/boostDM/mutrate", mode: "copy"
	label "core"

    input:
        tuple val(cohort), path(signature) from SIGNATURES5

    output:
		tuple val(cohort), path("*.mutrate.json") into OUT_MUTRATE

	script:
		"""
		parse-profile -i ${signature} -o ${cohort}
		"""

}


process DriverSaturation {
	tag "Driver saturation"
	publishDir "${STEPS_FOLDER}/boostDM/saturation", mode: "copy"
	label "core"

    input:
        path (drivers) from DRIVERS_SUMMARY

    output:
		path("*.vep.gz") into DRIVERS_SATURATION

	script:
		"""

		drivers-saturation --drivers ${drivers}
		
		"""
}

FILT_MNVS_INPUTS = PARSED_VEP7.map { it -> it[1] }
process FilterMNVS {
	tag "MNVs filter"
	publishDir "${STEPS_FOLDER}/boostDM", mode: "copy"
	label "core"

	input:
		path(input) from FILT_MNVS_INPUTS.collect()

	output:
		path("mnvs.tsv.gz") into MNVS_FILTER
	
	script:
		"""
		
		parse-mnvs ${input}

		"""
}



