#!/bin/bash

nextflow run ../.. -resume -profile apptainer --samplesheet sampleSheet.csv --samplesheet2 sampleSheet2.csv --outdir out