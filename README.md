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

### 3. Configure nextflow

Depending on your HPC/cloud setup you will need to adjust the nextflow profile in `nextflow.config`, to tell
nextflow how to submit the jobs. Using a `withName:...` directive, special
resources may be assigned to GPU-jobs. You can get an idea by checking out the `icbi_lung` profile - which we used to run the
workflow on our on-premise cluster. Only the `build_atlas` workflow makes use of GPU processes.

### 4. Launch the workflows

```bash
# Run `build_atlas` workflow
nextflow run main.nf --workflow build_atlas -resume -profile <YOUR_PROFILE> \
    --outdir "./data/20_build_atlas"

# Run `downstream_analysis` workflow
nextflow run main.nf --workflow downstream_analyses -resume -profile <YOUR_PROFILE> \
    --build_atlas_dir "./data/20_build_atlas" \
    --outdir "./data/30_downstream_analyses"
```

As you can see, the `downstream_analysis` workflow requires the output of the `build_atlas` workflow as input.
The intermediate results from zenodo contain the output of the `build_atlas` workflow.

## Structure of this repository

* `analyses`: Place for e.g. jupyter/rmarkdown notebooks, grouped by their respective (sub-)workflows.
* `bin`: executable scripts called by the workflow
* `conf`: nextflow configuration files for all processes
* `lib`: custom libraries and helper functions
* `modules`: nextflow DSL2.0 modules
* `tables`: contains static content that should be under version control
* `workflows`: the main nextflow workflows

## Contact

If you encounter any problems with the pipeline or have a feature request, please use the [issue tracker](https://github.com/Mye-InfoBank/atlas-pipeline/issues).
