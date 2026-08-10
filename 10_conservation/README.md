# 10 - Conservation

Nucleotide phyloP over CDS (`extract_conway_phylop.py`) and over the full gene
locus including introns (`extract_phylop_genomic.py`), from the UCSC 241-way and
470-way tracks. Per-isoform ConSurf runs on the cluster (`consurf_array.qsub`);
inputs are prepared by `make_consurf_isoform_inputs.py` and parsed by
`parse_consurf_isoforms.py`.
