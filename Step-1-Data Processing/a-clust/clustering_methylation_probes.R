clustering_methylation_probes <-function(methylationSNMnorm_usr_prb_matrix, methylationCoord, corr_dist_thresh_name, corr_dist_thresh, bp_dist_thresh, chromosomes_array, chrs_start_ind, chrs_end_ind, run_name){
  
  library(readr)
  library(R.matlab)
  source('Dbp.merge.R')
  source('Acluster.R')
  source('calc.dist.d.neighbor.R')
  source('update.clust.indicator.R')
  source('calc.dist.clusters.R')
  source('last.R')
  source('first.R')
  source('calc.dist.clust.point.R')
  

  methylationSNMnorm_usr_prb_matrix_sampled = matrix(nrow = 0,ncol = dim(methylationSNMnorm_usr_prb_matrix)[2])
  

  methylationCoord_sampled = matrix(nrow = 0,ncol = dim(methylationCoord)[2])
  

  
 
  #prbs_chr_num = methylationCoord[, 2];
  #chrs_start_ind = match(chromosomes_array, prbs_chr_num)
  #chrs_end_ind = 
  which.clust.final = c();
  max_cluster_number = 0;
  for(chr_num in chromosomes_array){
    chr_start_ind = chrs_start_ind[chr_num];
    chr_end_ind = chrs_end_ind[chr_num];
    #print(paste("chromosome number:",as.character(chr_num),"/",as.character(length(chromosomes_array))))
    #chr_end_ind = dim(prbs_chr_num)[1] - match(chr_num, rev(t(prbs_chr_num))) +1;
    methylationCoord_for_chr = methylationCoord[chr_start_ind:chr_end_ind, ]
    methylationSNMnorm_usr_prb_matrix_for_chr = methylationSNMnorm_usr_prb_matrix[chr_start_ind:chr_end_ind, ];
  
    indices_in_order = order(as.numeric(unlist(methylationCoord_for_chr[,3])));
    locations_in_order = methylationCoord_for_chr[indices_in_order,3];
    locations_in_order = unlist(locations_in_order);
    locations_in_order = as.numeric(locations_in_order);
    data_in_order = t(methylationSNMnorm_usr_prb_matrix_for_chr[indices_in_order,]);
  
    which.clust.initial_for_chr = Dbp.merge(ordr.vec = data_in_order, thresh.dist = corr_dist_thresh, bp.thresh.dist = bp_dist_thresh, location.vec = locations_in_order,  dist.type = "spearman")
    which.clust.final_for_chr = Acluster(ordr.vec = data_in_order, thresh.dist = corr_dist_thresh, which.clust = which.clust.initial_for_chr, location.vec = locations_in_order, max.dist = Inf, type = "average", dist.type = "spearman")

    samples_indices_sorted = match(unique(which.clust.final_for_chr), which.clust.final_for_chr);
    samples_indices_in_oder = indices_in_order[samples_indices_sorted];
    methylationCoord_for_chr_sampled = methylationCoord_for_chr[samples_indices_in_oder,];
    methylationSNMnorm_usr_prb_matrix_for_chr_sampled = methylationSNMnorm_usr_prb_matrix_for_chr [samples_indices_in_oder,]
    
    methylationCoord_sampled = rbind(methylationCoord_sampled, methylationCoord_for_chr_sampled);
    methylationSNMnorm_usr_prb_matrix_sampled = rbind(methylationSNMnorm_usr_prb_matrix_sampled, methylationSNMnorm_usr_prb_matrix_for_chr_sampled)
    
    #chr_start_ind = chr_end_ind + 1;
    
    which.clust.final_for_chr = which.clust.final_for_chr + max_cluster_number;
    
    which.clust.final = c(which.clust.final, which.clust.final_for_chr);
    
    max_cluster_number = which.clust.final[length(which.clust.final)];
    
  }
  
  dir.create("../../rosmapAD/data/coordinates/aclust_sampled")
  dir.create("../../rosmapAD/data/aclust_sampled")
  dir.create("../../rosmapAD/data/aclust_assignments")
  write.table(methylationCoord_sampled, file = paste("../../rosmapAD/data/coordinates/aclust_sampled/methylationCoord_mutual_sampled_corr_dist_thresh_",corr_dist_thresh_name,"_bp_dist_thresh_",as.character(bp_dist_thresh),"_",run_name,".csv", sep=""),row.names = FALSE,col.names = FALSE,sep=",")
  write.table(methylationSNMnorm_usr_prb_matrix_sampled, file = paste("../../rosmapAD/data/aclust_sampled/methylationSNMnorm_usr_prb_matrix_sampled_corr_dist_thresh_",corr_dist_thresh_name,"_bp_dist_thresh_",as.character(bp_dist_thresh),"_",run_name,".csv", sep=""),row.names = FALSE,col.names = FALSE,sep=",")
  write.table(which.clust.final, file = paste("../../rosmapAD/data/aclust_assignments/which_clust_final_corr_dist_thresh_",corr_dist_thresh_name,"_bp_dist_thresh_",as.character(bp_dist_thresh),"_",run_name,".csv", sep=""),row.names = FALSE,col.names = FALSE,sep=",")
}