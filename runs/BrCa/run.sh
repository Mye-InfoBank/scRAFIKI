#!/bin/bash

nextflow run ../.. -resume -profile cluster --samplesheet sampleSheetBRCA.csv --samplesheet2 sampleSheet2.csv --outdir out