#!/usr/bin/env bash


salmon quant -i transcripts_index -l <LIBTYPE> -r reads.fq --validateMappings -o transcripts_quant
