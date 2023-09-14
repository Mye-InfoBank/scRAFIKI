#!/bin/bash

nextflow run .. -resume -profile cluster --samplesheet sampleSheet.csv --samplesheet2 sampleSheet2.csv --outdir out