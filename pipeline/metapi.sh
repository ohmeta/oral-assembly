#!/usr/bin/env bash

conda install -c alienzj metapi
metapi init -d ./ -s samples.tsv -b raw -a metaspades
metapi run checkm_lineage_wf
