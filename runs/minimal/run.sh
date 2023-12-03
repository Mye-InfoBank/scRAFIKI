#!/bin/bash

nextflow run ../.. -resume -profile docker --samplesheet sampleSheet.csv --samplesheet2 sampleSheet2.csv --outdir out