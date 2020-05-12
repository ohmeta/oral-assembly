#!/usr/bin/env bash

# metapi
# https://github.com/ohmeta/metapi

# Installation
conda install -c alienzj metapi

# Init project
metapi init -d ./ -s samples.tsv -b trimming

# see what metapi can do
metapi denovo_wf --list

# dry run
metapi denovo_wf classify_all --dry_run

# submmit jobs
metapi denovo_wf classify_all --qsub
