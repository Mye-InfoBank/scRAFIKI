#!/bin/bash

nextflow run ../.. -resume -profile standard --samplesheet sampleSheetBRCA.csv --samplesheet2 sampleSheet2.csv --outdir out