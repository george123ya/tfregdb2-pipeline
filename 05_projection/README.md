# 05 - Domain projection

Project each curated domain onto every isoform. Exact sequence match first, then a
BLOSUM62 pairwise alignment; a domain is kept only when it matches at 100% identity
over 100% of its length (`map_canonical_to_isoforms.py`). DNA-binding domains come
directly from CIS-BP per isoform.
