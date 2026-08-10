# 05 - Domain projection

Project each curated domain onto every isoform. Exact sequence match first, then a
BLOSUM62 pairwise alignment; a domain is kept when it matches at 100% identity over
at least 95% of its length (`map_canonical_to_isoforms.py`). DNA-binding domains
come directly from CIS-BP per isoform.
