#!/bin/bash

nextflow run ../.. -resume -profile apptainer --samplesheet sampleSheetBRCA.csv --samplesheet2 sampleSheet2.csv --outdir out