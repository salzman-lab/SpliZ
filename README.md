## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**salzmanlab/spliz** is a bioinformatics best-practise analysis pipeline for calculating the splicing z-score for single cell RNA-seq analysis.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`) and [`conda`](https://docs.conda.io/en/latest/).

2. Download environment file.
    ```bash
    wget https://raw.githubusercontent.com/salzmanlab/SpliZ/main/environment.yml
    ```

3. Create conda environment and activate.
    ```bash
    conda env create --name spliz_env --file=environment.yml
    conda activate spliz_env
    ```

4. Run the pipeline on the test data set. 
You may need to modify the [executor scope](https://www.nextflow.io/docs/latest/executor.html) in the config file, in accordance to your compute needs.
    ```bash
    nextflow run salzmanlab/spliz \
        -r main \
        -latest \
        -profile small_test_data
    ```
    [Sherlock](https://www.sherlock.stanford.edu/) users should use the `sherlock` profile:


        nextflow run salzmanlab/spliz \
            -r main \
            -latest \
            -profile small_test_data,sherlock
            
 5. Run the pipeline on your own dataset.
    1. Edit your config file with the parameters below. (You can use `/small_data/small.config` as a template, be sure to include any memory or time paramters.)
    2. Run with your config file:
    ```
    nextflow run salzmanlab/spliz \
        -r main \
        -latest \
        -c YOUR_CONFIG_HERE.conf
    ```


See [usage docs](https://nf-co.re/spliz/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:
* Calculate the SpliZ scores for:
    * Identifying variable splice sites
    * Identifying differential splicing between cell types.

## Input Parameters

| Argument              | Description       |Example Usage  |
| --------------------- | ---------------- |-----------|
| `dataname`            | Descriptive name for SpliZ run        | "Tumor_5" |
| `run_analysis`        | If the pipeline will perform splice site identifcation and differential splicing analysis | `true`, `false` |
| `input_file`          | File to be used as SpliZ input | *tumor_5_with_postprocessing.txt* |
| `SICILIAN`            | If `input_file` is output from [SICILIAN](https://github.com/salzmanlab/SICILIAN) | `true`, `false` |
| `pin_S`               | Bound splice site residuals at this quantile (e.g. values in the lower `pin_S` quantile and the upper 1 - `pin_S` quantile will be rounded to the quantile limits) | 0.1 |
| `pin_z`               | Bound SpliZ scores at this quantile (e.g. values in the lower `pin_z` quantile and the upper 1 - `pin_z` quantile will be rounded to the quantile limits) | 0 |  
| `bounds`              | Only include cell/gene pairs that have more than this many junctional reads for the gene | 5 |
| `light`               | Only output the minimum number of columns | `true`, `false` |
| `svd_type`            | Type of SVD calculation | `normdonor`, `normgene` |
| `n_perms`             | Number of permutations | 100 |
| `grouping_level_1`    | Metadata column by which the data is intially partitioned  | "tissue" |
| `grouping_level_2`    | Metadata column by which the partitioned data is grouped | "compartment" |
| `libraryType`         | Library prepration method of the input data | `10X`, `SS2` |

## Optional Parameters for non-SICILIAN Inputs (`SICILIAN` = `false`)
| Argument              | Description       |Example Usage  |
| --------------------- | ---------------- |-----------|
| `samplesheet`         | If input files are in BAM format, this file specifies the locations of the input bam files. Samplesheet formatting is specified below. | *Tumor_5_samplesheet.csv* |
| `annotator_pickle`    | [Genome-specific annotation file for gene names](https://github.com/salzmanlab/SICILIAN#annotator-and-index-files-needed-for-running-sicilian) | *hg38_refseq.pkl* |
| `exon_pickle`         | [Genome-specific annotation file for exon boundaries](https://github.com/salzmanlab/SICILIAN#annotator-and-index-files-needed-for-running-sicilian) | *hg38_refseq_exon_bounds.pkl* |
| `splice_pickle`       | [Genome-specific annotation file for splice sites](https://github.com/salzmanlab/SICILIAN#annotator-and-index-files-needed-for-running-sicilian) | *hg38_refseq_splices.pkl* |
| `gtf`                 | GTF file used as the reference annotation file for the genome assembly | *GRCh38_genomic.gtf* |
| `meta`                | If input files are in BAM format, this file contains per-cell annotations. This file must contain columns for `grouping_level_1` and `grouping_level_2`. | *metadata_tumor_5.tsv* |

### Samplesheets

The samplesheet must be in comma-separated value(CSV) format. The file must be without a header. The sampleID must be a unique identifier for each bam file entry.

For non-SICILIAN samples, samplesheets must have 2 columns: sampleID and path to the bam file.
```
Tumor_5_S1,tumor_5_S1_L001.bam
Tumor_5_S2,tumor_5_S2_L002.bam
Tumor_5_S3,tumor_5_S3_L003.bam
```

For SICILIAN SS2 samples, amplesheets must have 3 columns: sampleID, read 1 bam file, and read 2 bam file.
```
Tumor_5_S1,tumor_5_S1_L001_R1.bam,tumor_5_S1_L001_R2.bam
Tumor_5_S2,tumor_5_S2_L002_R1.bam,tumor_5_S2_L002_R2.bam
Tumor_5_S3,tumor_5_S3_L003_R1.bam,tumor_5_S3_L003_R2.bam
```

## Credits

salzmanlab/spliz was originally written by Salzman Lab.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  salzmanlab/spliz for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->
This repositiory contains code to perform the analyses in this paper:

> **The SpliZ generalizes “Percent Spliced In” to reveal regulated splicing at single-cell resolution**
>
> Julia Eve Olivieri*, Roozbeh Dehghannasiri*, Julia Salzman.
>
> _Nature Methods_ 2022 Mar 3. doi: [10.1038/s41592-022-01400-x](https://doi.org/10.1038/s41592-022-01400-x).

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).


