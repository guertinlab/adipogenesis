#!/bin/bash
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH -p standard
#SBATCH -A janeslab
#SBATCH -o ChIP.call.peaks.out

module purge
module load macs2/2.2.7.1

macs2 callpeak -t CEBPbeta.sam -c Input.sam -f SAM -n CEBPbeta_ChIP --outdir CEBPbeta_macs --keep-dup 5 -q 0.01 -m 10 50
macs2 callpeak -t GR.sam -c Input.sam -f SAM -n GR_ChIP --outdir GR_macs --keep-dup 5 -q 0.01 -m 10 50

macs2 callpeak -t cJun.sam -f SAM -n cJun_ChIP --outdir cJun_macs --keep-dup 5 -q 0.01 -m 10 50
macs2 callpeak -t Fosl2.sam -f SAM -n Fosl2_ChIP --outdir Fosl2_macs --keep-dup 5 -q 0.01 -m 10 50
macs2 callpeak -t JunB.sam -f SAM -n JunB_ChIP --outdir JunB_macs --keep-dup 5 -q 0.01 -m 10 50
macs2 callpeak -t KLF4.sam -f SAM -n KLF4_ChIP --outdir KLF4_macs --keep-dup 5 -q 0.01 -m 10 50
macs2 callpeak -t KLF5.sam -f SAM -n KLF5_ChIP --outdir KLF5_macs --keep-dup 5 -q 0.01 -m 10 50

macs2 callpeak -t Tw2_DM_rep1.sam Tw2_DM_rep2.sam -c Input_DM.sam -f SAM -n Tw2_DM_ChIP --outdir Tw2_DM_macs --keep-dup 5 -q 0.05 -m 10 50
macs2 callpeak -t Tw2_GM_rep1.sam Tw2_GM_rep2.sam -c Input_GM.sam -f SAM -n Tw2_GM_ChIP --outdir Tw2_GM_macs --keep-dup 5 -q 0.05 -m 10 50

ls -d *macs/ > peak_files_list.txt

echo DONE
