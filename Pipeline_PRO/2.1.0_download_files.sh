#!/bin/bash

cd /scratch/bhn9by/primary_transcript_annotation
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks_PRO/misc_scripts/pTA.functions.R
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks_PRO/misc_scripts/overlaps_remove.csv
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks_PRO/misc_scripts/overlaps_keepadjacent.csv
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks_PRO/misc_scripts/duplicate_remove.csv
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks_PRO/misc_scripts/duplicate_keepadjacent.csv

#remove DOS \r\n\ artifact from .csv (if applicable)
sed -i 's/\r$//' *csv
