# Systemic Activation Analysis Pipeline

This repository contains all scripts and resources used to generate the data and figures for the preprint:

**“Adrenergic signaling coordinates distant and local responses to amputation in axolotl”** *(preprint link TBD)*

---

## Quick Start Instructions

To generate the figures without reprocessing raw data:

1. Download the preprocessed Seurat object: [seurat_filtered_publication.rds](https://doi.org/10.7910/DVN/NW0ZN1)

   ** note : currently this doi link does not work, but can be accessed with limited permissions via this link (https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NW0ZN1)

2. Run the figure generation script:  `./scripts/4_producing_figures.Rmd`


---

## Full Reproduction Instructions

To reproduce the analysis end-to-end, run the scripts in `./scripts/` in the following order:

   ** note : prepare scanpy_env using `./bin/scanpy_env.yml` prior to running scripts 1 - 2.2

| Step | Script                                   | Purpose                                                                 |
|------|------------------------------------------|-------------------------------------------------------------------------|
| 0    | `0_dataverse_download.sh`                | Downloads raw FASTQ files (`./data/`), reference files (`./ref/`), and Seurat objects (`./output/seurat/`) from Dataverse |
| 1    | `1_indexing.sh`                          | Builds the transcriptome index using `kb ref`                          |
| 2.1  | `2.1_kb_count_2C.sh`                     | Aligns and quantifies the 2C sequencing data                           |
| 2.2  | `2.2_kb_count_4C.sh`                     | Aligns and quantifies the 4C sequencing data                           |
| 3    | `3_creating_and_filtering_seurat_object.Rmd` | Creates and filters a unified Seurat object for both 2C and 4C       |
| 4    | `4_producing_figures.Rmd`                | Generates all figures for the manuscript                               |

---

## Computational Environment

All analyses were conducted using the [Harvard FASRC Cluster](https://www.rc.fas.harvard.edu/) in 2025.

Software versions:

- Python v3.9.21 - see `./bin/scanpy_env.yml` for environment info

- - R v4.4.1 - see `./bin/session_info.txt` for package list

- kb-python v0.29.1

- kallisto v0.51.1 — used for reference index construction and pseudoalignment

- bustools v0.44.1 — used for quantification

---

## Data Availability

All raw and processed data — including FASTQ files, references, and Seurat objects — are available via: Harvard Dataverse https://doi.org/10.7910/DVN/NW0ZN1

---

## Contact & Authorship

This repository is maintained by members of the **Whited Lab** at the **Harvard University Department of Stem Cell and Regenerative Biology**.

**Creator**  
- **Name:** Tim Froitzheim  
- **Email:** [tfroitzheim@fas.harvard.edu](mailto:tfroitzheim@fas.harvard.edu)  
- **Role:** Lead Bioinformatician


**Author Contact**  
- **Name:** Dr. Duygu Payzin-Dogru
- **Email:** [duygupayzin@g.harvard.edu](mailto:duygupayzin@g.harvard.edu)  
- **Role:** Project Lead Scientist


**Principal Investigator**  
- **Name:** Dr. Jessica L. Whited  
- **Lab Website:** [www.whitedlab.com](http://www.whitedlab.com)  
- **Lab Email:** [whitedlab@gmail.com](mailto:whitedlab@gmail.com)


**Repository Contact**  
- **Name:** Hani Singer
- **Email:** [hani_singer@fas.harvard.edu](mailto:hani_singer@fas.harvard.edu)  
- **Role:** Research Laboratory Manager

---

## Issues & Contributions

For questions, bug reports, or to contribute:
- Open an Issue at [https://github.com/Whited-Lab/project-systemic-activation/issues](https://github.com/Whited-Lab/project-systemic-activation/issues)
- Or contact us directly by email
