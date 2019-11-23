# adipogenesis
LRT_clustering.R follows my post peaks calling workflow:

3T3_atac_summit_200window.bed is the final peak calling file
you can start with line 35 because I uploaded df.preadipo.Rdata 

At line 76 I use LRT_rivanna_slurm.slurm to run LRT_largemem_clustering.R on rivanna
the bed files for initial clustering are here: initial_clusters_bed.zip

I load 191119_clusters.all.minc100.pval0.00000001.Rdata from rivanna and proceed with LRT_clustering.R line 108
the bed files for the final six clusters are here: final_SIX_clusters_bed.zip

at line 180 I ran all the initial clusters on rivanna using 191122_meme_all_cluster*.slurm and 191122_meme_cluster*.sh files

at line 404 I ran all the final six clusters on rivanna using 191122_meme_all_cluster\*gradual\*slurm and  191122_meme_cluster\*gradual\*sh


