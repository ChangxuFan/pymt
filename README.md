# The core code for Changxu Fan & Jun Wu, et al. 2024
## Dependencies:
The analysis workflow heavily depends on the following packages I wrote:
* [scFanc](https://github.com/ChangxuFan/scFanc) Release v1.0.1
* [liteRnaSeqFanc](https://github.com/ChangxuFan/liteRnaSeqFanc) Release v1.0.1

To begin analysis, you need to copy/softlink the .R files from these packages to the corresponding directories (R_scFanc, R_liteRnaSeqFanc).

Although not critical, some functions in this repo might call functions from other packages I wrote, such as `utilsFanc`, `abaFanc2`, etc. These can also be found in the corresponding github repos under the account `ChangxuFan`. 
In addition, some stand-alone scripts (`*.R`, `*.sh`) are used, and can be found in the `R_for_bash` or `scripts` repos under the account `ChangxuFan`.

## Sequence of analyses:
* Data processing starts at `fastcheck/step0.1_cellranger_count.sh`, where cellranger was used to align the reads.
* Scripts under `fast_check` were run first to perform clustering on a single-sample basis
* Then, scripts under `sync_all` were used to generate R objects containing all samples with batch correction
* Then, scripts under `de_noNZ` and `da_noNZ` were run to perform cluster-specific differential expression and differential accessibility analyses.
* `publication` contains scripts that generated the figures in the paper. 
* In general, scripts should be run sequentially, with step1.1 before step1.2. 

## Stochasticity of the pipeline:
* We noticed that the doublet detection process has some stochasticity: if you run it twice, the doublet scores will not be exactly the same, although they are quite similar. Therefore, after doublet filtering, you are likely to end up with a slightly different set of cells. Unfortunately, such slight differences might cause ArchR to assign different cluster names. We offer the integrated Seurat/ArchR objects (referred to as `soi` and `aoi` in the code, respectively) at the corresponding GEO DataSet entry of this paper, which contain the cells obtained after our run of doublet filtering and their cluster assignments.

## About dead symlinks:
* R_scFanc and R_liteRnaSeqFanc link to the /R directories under scFanc and liteRnaSeqFanc, respectively. You can download them from the corresponding repositories, as mentioned above.
* Other dead symlinks point to large files/directories that github cannot handle. But you can generate them by running the code.
