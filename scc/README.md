# HPC (SGE) jobs

Three stages were run on an SGE cluster (Sun Grid Engine / `qsub`) rather than a
workstation, because they need many GPUs, a large sequence database, or a lot of
memory. The job scripts live with their stage, not here; this file just lists
them and the environment they expect. Each script is a template: set the
placeholders (`<PROJECT>`, `PIPELINE_DIR`, `CONSURF_BASE`) for your own cluster,
since paths, project codes, and GPU names are site-specific.

| Stage | Script | Job | Resources |
| --- | --- | --- | --- |
| 06 structure | `06_structure/af2_features.qsub` | AlphaFold MSA / features (CPU array) | 8 cores, 8 GB/core |
| 06 structure | `06_structure/af2_inference.qsub` | AlphaFold GPU inference | 1 GPU (L40S; A100-80G for >2700 aa) |
| 06 structure | `06_structure/extract_pae.qsub` | Extract PAE matrices (CPU array) | 4 cores |
| 10 conservation | `10_conservation/consurf_array.qsub` | Per-isoform ConSurf (4820-task array) | 4 cores |
| 11 interactions | `11_interactions/extract_intermaps.qsub` | Quantize finches epsilon matrices | 8 cores, high memory |

## Environment

- **AlphaFold**: modules `alphafold/2.3.1` + `cuda/11.8`; ParallelFold
  `run_alphafold.sh` wrapper; preset `monomer_ptm`, `full_dbs`, template cutoff
  `2022-04-14`. Inference on L40S (`sm_89`); proteins over 2700 aa on A100-80G
  (`sm_80`, which can also GPU-relax). Relax runs on CPU on the L40S (`-G`).
- **ConSurf**: standalone ConSurf in a conda env named `consurf`; HMMER search
  against UniRef90, up to 150 homologs (35-95% identity), one iteration.
- **finches intermaps**: any env with numpy; loads the ~28 GB interaction pickle
  in one process, so it needs a high-memory node.

Only 475 sequences were folded in house (15 canonical proteins missing from the
AlphaFold Database plus 460 alternative isoforms); everything else was downloaded
from the AlphaFold Database.
