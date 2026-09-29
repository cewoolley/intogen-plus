
methylation_data_srcdir = ${src_datasets}/methylation

methylation_dir = $(INTOGEN_DATASETS)/methylation

$(methylation_dir): | $(INTOGEN_DATASETS)
	mkdir $@

# Illumina methylation arrays annotation (hg38) from https://zwdzwd.github.io/InfiniumAnnotation
methylation_annotation_url = https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno
methylation_arrays ?= HM450 EPIC EPICv2
methylation_gencode ?= v41
# promoter window around the TSS of the MANE transcript
methylation_promoter_upstream ?= 1500
methylation_promoter_downstream ?= 500

METHYLATION_MANIFESTS = $(foreach array,$(methylation_arrays),$(methylation_dir)/$(array).hg$(genome).manifest.gencode.$(methylation_gencode).tsv.gz)
METHYLATION_MASKS = $(foreach array,$(methylation_arrays),$(methylation_dir)/$(array).hg$(genome).mask.tsv.gz)

$(methylation_dir)/%.hg$(genome).manifest.gencode.$(methylation_gencode).tsv.gz: | $(methylation_dir)
	@echo Downloading $* annotation
	curl -L -s -f -o $@ "${methylation_annotation_url}/$*/$*.hg$(genome).manifest.gencode.$(methylation_gencode).tsv.gz"

$(methylation_dir)/%.hg$(genome).mask.tsv.gz: | $(methylation_dir)
	@echo Downloading $* masks
	curl -L -s -f -o $@ "${methylation_annotation_url}/$*/$*.hg$(genome).mask.tsv.gz"


METHYLATION_PROBES = $(methylation_dir)/promoter_probes.tsv.gz

$(METHYLATION_PROBES): ${methylation_data_srcdir}/promoter_probes.py $(METHYLATION_MANIFESTS) $(METHYLATION_MASKS) $$(BIOMART_CDS) | $(methylation_dir)
	@echo Building promoter probes of the methylation arrays
	python $< \
		--biomart $(BIOMART_CDS) \
		$(foreach f,$(METHYLATION_MANIFESTS),--manifest $(f)) \
		$(foreach f,$(METHYLATION_MASKS),--mask $(f)) \
		--upstream $(methylation_promoter_upstream) \
		--downstream $(methylation_promoter_downstream) \
		-o $@


DATASETS += $(METHYLATION_PROBES)
