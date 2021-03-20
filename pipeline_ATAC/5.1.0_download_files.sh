#!/bin/bash

cd /scratch/bhn9by/ATAC/
#for main figure
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks/misc_scripts/generate_composite_motif.SP.R
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks/misc_scripts/generate_composite_motif.KLF.R
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks/misc_scripts/prep.SP.KLF.fimo.R
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks/misc_scripts/SP_KLF_split.R
wget https://raw.githubusercontent.com/guertinlab/adipogenesis/master/Bao_code_chunks/misc_scripts/extract.motifs.from.combined.family.R

#for supplemental figure
wget https://github.com/guertinlab/adipogenesis/blob/master/Bao_code_chunks/misc_scripts/post.composite.fimo.supp.R
wget https://github.com/guertinlab/adipogenesis/blob/master/Bao_code_chunks/misc_scripts/plot.motif.enrichment.supp.R
