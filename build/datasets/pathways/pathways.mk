
pathways_data_srcdir = ${src_datasets}/pathways

pathways_dir = $(INTOGEN_DATASETS)/pathways

$(pathways_dir): | $(INTOGEN_DATASETS)
	mkdir $@

# Reactome (CC0). The current release is downloaded; set reactome_url to pin one
reactome_url ?= https://reactome.org/download/current/ReactomePathways.gmt.zip
# MSigDB hallmark gene sets (CC BY 4.0)
msigdb_version ?= 2024.1.Hs
hallmarks_url ?= https://data.broadinstitute.org/gsea-msigdb/msigdb/release/$(msigdb_version)/h.all.v$(msigdb_version).symbols.gmt

REACTOME_GMT = $(pathways_dir)/ReactomePathways.gmt.zip
$(REACTOME_GMT): | $(pathways_dir)
	@echo Downloading Reactome pathways
	curl -L -s -f -o $@ "$(reactome_url)"

HALLMARKS_GMT = $(pathways_dir)/h.all.v$(msigdb_version).symbols.gmt
$(HALLMARKS_GMT): | $(pathways_dir)
	@echo Downloading MSigDB hallmarks
	curl -L -s -f -o $@ "$(hallmarks_url)"

GENE_SETS = $(pathways_dir)/gene_sets.tsv.gz
$(GENE_SETS): ${pathways_data_srcdir}/gene_sets.py $(REACTOME_GMT) $(HALLMARKS_GMT) $$(BIOMART_CDS) $$(SYMBOLS_MAP) | $(pathways_dir)
	@echo Building gene sets
	python $< \
		--reactome $(REACTOME_GMT) \
		--hallmarks $(HALLMARKS_GMT) \
		--biomart $(BIOMART_CDS) \
		--symbols $(SYMBOLS_MAP) \
		-o $@


DATASETS += $(GENE_SETS)
