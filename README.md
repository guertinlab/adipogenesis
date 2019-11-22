# adipogenesis
LRT_clustering.R follows my post peaks calling workflow
3T3_atac_summit_200window.bed is the final peak calling file

Next I use LRT_rivanna_slurm.slurm to run LRT_largemem_clustering.R on rivanna

I load 191119_clusters.all.minc100.pval0.00000001.Rdata from rivanna and proceed with LRT_clustering.R line 108

at line 180 I ran all the initial clusters on rivanna using 191122_meme_all_cluster*.slurm and 191122_meme_cluster*.sh files


