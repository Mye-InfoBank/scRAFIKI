# SIMBA - **S**ingle-cell **I**ntegration **M**ethods pipeline for **B**ig **A**tlases

<img src="./www/SIMBA_sticker.png" width="250">

SIMBA isvgs a pipeline designed to integrate single-cell RNA-seq datasets to generate atlases. It incorporates several integration methods, cell filtering and automated annotation tools.

This pipeline is written in nextflow and based on the LUCA pipeline:
  
> Salcher, S., Sturm, G., Horvath, L., Untergasser, G., Kuempers, C., Fotakis, G., ... & Trajanoski, Z. (2022). High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer. Cancer Cell. [doi:10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008)

## Overview

![Metro map](./www/metro.png)
Created using [Tennessine](https://tennessine.co.uk/metro/)

## Usage

### Prerequisites

* [Nextflow](https://www.nextflow.io/index.html#GetStarted), version 21.10.6 or higher
* [Docker](https://docs.docker.com/get-docker/), [Apptainer](https://apptainer.org/docs/admin/main/installation.html) or one of the other [supported container engines](https://www.nextflow.io/docs/latest/container.html)

### 2. Prepare data

You have to create AnnData objects for each dataset, containing the UMI counts for each cell. The data need to be stored in the h5ad format.

The following metadata fields are required:
| Field | Description | Axis | Default |
| --- | --- | --- | --- |
| batch | Batch identifier, for integration | obs | *required* |
| celltype | Cell-type annotation | obs | unknown |
| condition | The condition of the tissue sample | obs | unknown |
| sex | The sex of the patient (`female` or `male`) | obs | unknown |
| patient | The patient identifier | obs | *required* |
| tissue | The tissue type | obs | *required* |

Patient and batch identifiers will be prepended with the dataset identifier to ensure uniqueness.
Only alphanumeric characters and underscores are allowed for all metadata fields.

### 3. Pipeline execution

```bash
nextflow run Mye-InfoBank/SIMBA -resume -profile <YOUR_PROFILE> --outdir "results" --samplesheet "samplesheet.csv"
```

#### Profiles
The following profiles are available:
- `standard`: Run the pipeline on your local machine without containerization, requires all dependencies to be installed (not recommended)
- `docker` (if running on an ARM-based machine, use `docker,arm`)
- `singularity`
- `podman`
- `shifter`
- `charliecloud`
- `apptainer`
- `gitpod`

#### Samplesheet
The samplesheet is a csv file with the following column names as header:
| Column | Description | Default |
| --- | --- | --- |
| id | The dataset identifier, needs to be unique across datasets | *required* |
| input_adata | Path to the input AnnData object | *required* |
| min_counts | Minimum number of counts per cell | 0 |
| max_counts | Maximum number of counts per cell | Infinity |
| min_genes | Minimum number of genes per cell | 0 |
| max_genes | Maximum number of genes per cell | Infinity |
| max_pct_mito | Maximum percentage of mitochondrial genes per cell | 100 |
| run_solo | Run SOLO for doublet detection | true |

Cell filtering will be handled by the pipeline if the respective columns are specified. If filtering has been done before, the columns can be omitted.

## Contact

If you encounter any problems with the pipeline or have a feature request, please use the [issue tracker](https://github.com/Mye-InfoBank/atlas-pipeline/issues).
