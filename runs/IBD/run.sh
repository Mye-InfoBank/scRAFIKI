#!/bin/bash

nextflow run ../.. -resume -profile standard --samplesheet sampleSheet.csv --samplesheet2 sampleSheet2.csv --outdir out