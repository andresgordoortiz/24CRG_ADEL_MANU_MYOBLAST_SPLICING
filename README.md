# 24CRG_ADEL_MANU_MYOBLAST_SPLICING
**Author**: Andrés Gordo Ortiz
**Institution**: Centre for Genomic Regulation (CRG)
**Supervisors**: Adel Al Jord
**Project Type**: Master Thesis

## Project Overview

This repository consists on the data exploration, processing and plotting for the Myoblast Splicing Analysis at different stages of development to myotubes  for the ***Al Jord Lab @CRG*** during my master thesis.

## Samples

| Name              | Publisher                             | ENA Link  | Year | Study Specs                                                                                                                                                                                           | Details                                                 |
|-------------------|---------------------------------------|-----------|------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| Kuo Zhang         | Nature Communications                 | [SRP160431](https://www.ebi.ac.uk/ena/browser/view/SRP160431) | 2018 | Reads from knock-down and shControl cells in myoblast and myotube stages. Controls include 3 paired-end biological replicates each for myoblast and myotubes, totaling 12 fastq files.                | ~30M per biosample                                      |
| Lingjian Tao      | Molecular Biotechnology               | [SRP496253](https://www.ebi.ac.uk/ena/browser/view/SRP496253) | 2024 | Myotube formation observed after 4 days, full maturation at 7 days. 3 paired-end biological replicates per condition across 4 timepoints (0, 2, 4, 7 days), totaling 24 fastq files.                  | ~44M per biosample                                      |
| Pengcheng Lyu     | MDPI Cells (data quality verified)    | [SRP335878](https://www.ebi.ac.uk/ena/browser/view/SRP335878) | 2022 | C2C12 cells differentiated in media, with autophagy inhibition using chloroquine (10 µM or 20 µM). 2 paired-end biological replicates per condition for myoblast and myotube stages, 8 fastq files.   | ~45M per biosample                                      |
| Dominic W. Kolonay | Frontiers in Cell and Developmental Biology | [SRP471123](https://www.ebi.ac.uk/ena/browser/view/SRP471123) | 2024 | Myoblasts differentiated to myotubes in DMEM with horse serum and antibiotics over 5 days. 3 paired-end biological replicates per condition (0 days myoblast, 5 days myotube), totaling 12 fastq files. | ~90M per biosample in myotubes, ~140M per biosample in myoblasts |
| Christopher Azar  | Physiological Reports                 | [SRP198848](https://www.ebi.ac.uk/ena/browser/view/SRP198848) | 2021 | Cells differentiated for 6 days with daily media changes. 3 single-end biological replicates per condition (0 days myoblast, 6 days myotube), totaling 6 fastq files.                                 | ~30M per biosample


## Repository Structure

- **config/**: Configuration files for Docker, and HPC cluster settings.
- **data/**: Raw, processed and metadata files. Raw and processed data are excluded from version control.
- **notebooks/**: Exploratory and analysis notebooks for VSCode and RStudio.
- **scripts/**: Custom scripts for data processing in Bash and downstream analysis in R.
- **workflows/**: Complete Pipelines to process the data from zero (no data in local) up to the finl Inclusion tables with all events. We recommend using *SLURM* for batch processing, otherwise the pipelines and scripts will have to be tweaked in order to work in your own implementation of work management.
- **docs/**: Documentation files, including installation instructions and workflow details.
- **.github/**: GitHub-specific workflows for continuous integration and issue templates.

## Processing Workflow

1. All samples *fastq.gz* files were quality checked using fastQC and MultiQC [[1]](#1)
2. Reads were then **mapped** to the latest [mm10 build](https://vastdb.crg.eu/libs/vastdb.mm2.23.06.20.tar.gz) of the mouse genome using the Bowtie2 implementation of **Vast Tools**, which also allows for curated identification and quantification of splicing events. See [Vast-tools](https://github.com/vastgroup/vast-tools) for indications on how to manually download the builds [[2]](#2).
3. Inclusion Tables from each study were then transferred to RStudio for downstream processing. Simulations through the **Beta distribution** and statistical significance (*p-value <= 0.05*) for each splicing event were calculated using the [betAS](https://github.com/DiseaseTranscriptomicsLab/betAS/) package [[3]](#3).

## Installation and Setup
In order to run this repo without forking and modifying it, you will need access to a HPC Cluster which uses *SLURM* as the work manager and any Linux-based distribution as OS, although *AlmaLinux* was used in our case.

1. Clone the repository:
   ```bash
   # Clone the repository to your HPC folder
   git clone https://github.com/andresgordoortiz/24CRG_ADEL_MANU_MYOBLAST_SPLICING.git
   cd 24CRG_ADEL_MANU_MYOBLAST_SPLICING
   ```

2. Run the pipelines using *SLURM*. All code runs under online *Docker* images translated into Singularity *SIF* files. Therefore, no specific software needs to be installed beforehand.
   ```bash
   # Important: you must pass a suitable VASTDB database as absolute path to run the pipelines
   sbatch workflows/full_processing_pipeline_christopher.sh /PATH_TO_VASTDB
   sbatch workflows/full_processing_pipeline_dominic.sh /PATH_TO_VASTDB
   sbatch workflows/full_processing_pipeline_lyu_cells.sh /PATH_TO_VASTDB
   sbatch workflows/full_processing_pipeline_tao_molbiotech.sh /PATH_TO_VASTDB
   sbatch workflows/full_processing_pipeline_zhangnature.sh /PATH_TO_VASTDB
   ```

3. Run the RMarkdown Notebook
   ```bash
   sbatch scripts/R/run_notebook.sh

## Docker Containers used
1. Rmarkdown Processing: [andresgordoortiz/splicing_analysis_r_crg:v1.2](https://hub.docker.com/layers/andresgordoortiz/splicing_analysis_r_crg/v1.2/images/sha256-67fd933eb88fbb9e3fe099c18eef9926bcaf916f92ff0f1fd5f9e876f78fd726?context=repo)
2. Vast-Tools Alternative Splicing Analysis: [andresgordoortiz/vast-tools:latest](https://hub.docker.com/layers/andresgordoortiz/vast-tools/latest/images/sha256-e760bb36d7383ad9d9447035d09d1b282f52e8d44acf6a14ffc23ffcc3d7d383?context=repo)

## References
<a id="1">[1]</a>
Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354

<a id="2">[2]</a>
Gohr, A., Mantica, F., Hermoso-Pulido, A., Tapial, J., Márquez, Y., Irimia, M. (2022). Computational Analysis of Alternative Splicing Using VAST-TOOLS and the VastDB Framework. In: Scheiffele, P., Mauger, O. (eds) Alternative Splicing. Methods in Molecular Biology, vol 2537. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-2521-7_7

<a id="3">[3]</a>
Ascensão-Ferreira, M., Martins-Silva, R., Saraiva-Agostinho, N. & Barbosa-Morais, N. L. betAS: intuitive analysis and visualization of differential alternative splicing using beta distributions. RNA 30, 337 (2024).


