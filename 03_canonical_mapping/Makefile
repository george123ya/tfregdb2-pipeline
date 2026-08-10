# tfregdb-blast pipeline — see README.md for the high-level flow.
# Each stage's output is idempotent: rerun anytime, only the steps whose
# inputs are newer than their outputs will actually do work.

PY := python3
INPUT  := data/input
OUTPUT := data/output

HYDRATED       := $(OUTPUT)/hydrated_HumanizedTable.xlsx
HYDRATED_ALL   := $(OUTPUT)/hydrated_HumanizedTable_with_animals.xlsx
ANIMAL_CACHE   := $(OUTPUT)/.animal_cache.json
QUERIES        := $(OUTPUT)/queries.fasta
DB             := $(OUTPUT)/custom_db/human_proteome_db.phr
BLAST_TSV      := $(OUTPUT)/blast_results.tsv
UNIFIED        := $(OUTPUT)/Unified_Valid_Domains.xlsx

.PHONY: all enrich fetch-animals queries db blast project clean help

help:
	@echo "Targets:"
	@echo "  make enrich        — hydrate HumanizedTable HUMAN rows via TF_completeIDs"
	@echo "  make fetch-animals — hydrate non-human rows via UniProt REST API (network)"
	@echo "  make db            — build the human-proteome BLAST database"
	@echo "  make queries       — emit queries.fasta from (fully-hydrated) tables"
	@echo "  make blast         — run blastp queries → results.tsv"
	@echo "  make project       — project hits onto canonical proteins → Unified xlsx"
	@echo "  make all           — full pipeline (enrich + fetch-animals + db + queries + blast + project)"
	@echo "  make clean         — remove output/ (keeps input/)"

all: $(UNIFIED)

enrich: $(HYDRATED)
$(HYDRATED): $(INPUT)/HumanizedTable.xlsx $(INPUT)/TF_completeIDs_with_ensembl_canonical.csv scripts/00_enrich_canonical.py
	$(PY) scripts/00_enrich_canonical.py \
	    --humanized $(INPUT)/HumanizedTable.xlsx \
	    --canonical $(INPUT)/TF_completeIDs_with_ensembl_canonical.csv \
	    --out       $(HYDRATED) \
	    --failures  $(OUTPUT)/enrichment_failures.tsv

fetch-animals: $(HYDRATED_ALL)
$(HYDRATED_ALL): $(HYDRATED) scripts/00b_fetch_animal_sequences.py
	$(PY) scripts/00b_fetch_animal_sequences.py \
	    --hydrated $(HYDRATED) \
	    --out      $(HYDRATED_ALL) \
	    --cache    $(ANIMAL_CACHE) \
	    --failures $(OUTPUT)/animal_fetch_failures.tsv

# Build the BLAST DB once. makeblastdb writes .phr/.pin/.psq alongside.
db: $(DB)
$(DB): $(INPUT)/human_proteome.fasta
	mkdir -p $(OUTPUT)/custom_db
	makeblastdb -in $< -dbtype prot -out $(OUTPUT)/custom_db/human_proteome_db

queries: $(QUERIES)
$(QUERIES): $(HYDRATED_ALL) $(INPUT)/TFRegDB1.xlsx scripts/01_prepare_queries.py
	$(PY) scripts/01_prepare_queries.py \
	    --humanized $(HYDRATED_ALL) \
	    --tfregdb1  $(INPUT)/TFRegDB1.xlsx \
	    --out       $(QUERIES)

blast: $(BLAST_TSV)
$(BLAST_TSV): $(QUERIES) $(DB) scripts/02_run_blast.py
	$(PY) scripts/02_run_blast.py \
	    --query $(QUERIES) \
	    --db    $(OUTPUT)/custom_db/human_proteome_db \
	    --out   $(BLAST_TSV)

project: $(UNIFIED)
$(UNIFIED): $(BLAST_TSV) $(INPUT)/TF_completeIDs_with_ensembl_canonical.csv scripts/03_project_canonical.py
	$(PY) scripts/03_project_canonical.py \
	    --blast      $(BLAST_TSV) \
	    --canonical  $(INPUT)/TF_completeIDs_with_ensembl_canonical.csv \
	    --humanized  $(HYDRATED_ALL) \
	    --tfregdb1   $(INPUT)/TFRegDB1.xlsx \
	    --out        $(UNIFIED) \
	    --failures   $(OUTPUT)/blast_failures.tsv

clean:
	rm -rf $(OUTPUT)
	@echo "Cleaned $(OUTPUT) (kept $(INPUT))"
