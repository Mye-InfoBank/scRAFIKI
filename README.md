# SIMBA - **S**ingle-cell **I**ntegration **M**ethods pipeline for **B**ig **A**tlases

<img src="./www/SIMBA_sticker.png" width="250">

SIMBA is a pipeline designed to integrate single-cell RNA-seq datasets to generate atlases. It incorporates several integration methods, cell filtering and automated annotation tools.

This pipeline is written in nextflow and based on the LUCA pipeline:
  
> Salcher, S., Sturm, G., Horvath, L., Untergasser, G., Kuempers, C., Fotakis, G., ... & Trajanoski, Z. (2022). High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer. Cancer Cell. [doi:10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008)

## Overview

The goal of the pipeline is to check the quality of the input data, integrate the datasets, and annotate the cell-types.

The following quality control steps are performed:
- Dataset-specific cutoff for number of detected genes and read counts as well as mitochondrial fraction
- Doublet detection with [SOLO](https://docs.scvi-tools.org/en/stable/api/reference/scvi.external.SOLO.html)
- Removal of ambient RNA

The following integration methods are available:
- scVI
- scANVI
- Harmony
- MNN

The resulting dataset can be annotated using cellypist. The result file is prepared for further investigation using cellxgene.

## Usage

### Prerequisites

* [Nextflow](https://www.nextflow.io/index.html#GetStarted), version 21.10.6 or higher
* [Docker](https://docs.docker.com/get-docker/), [Apptainer](https://apptainer.org/docs/admin/main/installation.html) or one of the other [supported container engines](https://www.nextflow.io/docs/latest/container.html)

### 2. Prepare data

You have to create AnnData objects for each dataset, containing the UMI counts for each cell. The data need to be stored in the h5ad format.

The following metadata fields are required:
| Field | Description | Axis | Default |
| --- | --- | --- | --- |
| batch | Batch identifier, for integration | obs | |
| cell_type | Cell-type annotation | obs | unknown |
| condition | The condition of the tissue sample | obs | unknown |
| sex | The sex of the patient (`female` or `male`) | obs | unknown |
| patient | The patient identifier, needs to be unique across datasets | obs |  |
| tissue | The tissue type | obs | |

### 3. Pipeline execution

```bash
nextflow run Mye-InfoBank/SIMBA -resume -profile <YOUR_PROFILE> --outdir "results" --samplesheet "samplesheet.csv"
```

## Contact

If you encounter any problems with the pipeline or have a feature request, please use the [issue tracker](https://github.com/Mye-InfoBank/atlas-pipeline/issues).
