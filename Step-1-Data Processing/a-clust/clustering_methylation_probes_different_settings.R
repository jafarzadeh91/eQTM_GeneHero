library(readr)
library(plotly)
library(doParallel)
library(foreach)
library(bigmemory)
library('log4r')
debugSource("clustering_methylation_probes.R")

logger <- create.logger()
logfile(logger) <- 'base.log'
level(logger) <- 'INFO'

msg = paste('number of available cores:',as.character(detectCores()));
print(msg)
info(logger, msg)

cl <- makeCluster(detectCores())
registerDoParallel(cl)


msg = paste('number of assigned cores:',as.character(detectCores()));
print(msg);
info(logger, msg);

msg = 'Running has been started....'
print(msg)
info(logger, msg)

#corr_dist_thresh_array <- eval(parse(text=args[1]));
#corr_dist_thresh_name_array <- eval(parse(text=args[2]));
#bp_dist_thresh_array <- eval(parse(text=args[3]));

corr_dist_thresh_array = c(0.1, 0.2, 0.4, 0.8);
corr_dist_thresh_name_array = c('point1', 'point2', 'point4', 'point8');
bp_dist_thresh_array = c(500, 1000, 2000, 4000);

chromosome_array_whole = 1:22;

chromosome_array_train = 1:18;
chromosome_array_test = 19:22;

train_test_array = c("train", "test")

run = TRUE;
report = FALSE;
ploting = FALSE;

  
if(run == TRUE || report == TRUE){
  
  methylationCoord <- read_csv("../../rosmapAD/data/coordinates/methylationCoord_mutual.csv", col_names = FALSE)
  ## we need to sort matrices according to the probes positions first of all
  probes_indices_in_order=order( unlist(methylationCoord[,2]), unlist(methylationCoord[,3]) );
  methylationCoord = methylationCoord[probes_indices_in_order,];
  
  
  methylationSNMnorm_usr_prb_matrix <- read_csv("../../rosmapAD/data/methylationSNMnorm_usr_prb_matrix.csv", col_names = FALSE)
  methylationSNMnorm_usr_prb_matrix= methylationSNMnorm_usr_prb_matrix[probes_indices_in_order,];
  methylationSNMnorm_usr_prb_matrix=as.big.matrix(as.matrix(methylationSNMnorm_usr_prb_matrix,type = double));
  methylationSNMnorm_usr_prb_matrix_desciption <-describe(methylationSNMnorm_usr_prb_matrix);
  methylationSNMnorm_usr_prb_matrix_pointer = attach.big.matrix(methylationSNMnorm_usr_prb_matrix_desciption);


  probes_chromosome_number = as.integer(unlist(methylationCoord[, 2]));
  chrs_start_ind = match(chromosome_array_whole, t(methylationCoord[, 2]))
  chrs_end_ind = length(probes_chromosome_number)-match(chromosome_array_whole,rev(probes_chromosome_number))+1
}

if(run){
  foreach (corr_dist_thresh_ind = 1:length(corr_dist_thresh_array), .packages=c("foreach","bigmemory", "log4r"))%:%
    foreach(bp_dist_thresh_ind = 1:length(bp_dist_thresh_array),.packages=c("foreach","bigmemory", "log4r"))%dopar%{
      ## outer loop code
      corr_dist_thresh_name = corr_dist_thresh_name_array[corr_dist_thresh_ind];
      corr_dist_thresh = corr_dist_thresh_array[corr_dist_thresh_ind]
      bp_dist_thresh = bp_dist_thresh_array[bp_dist_thresh_ind];

      
      msg = paste('corr_tresh:',as.character(corr_dist_thresh), '&','bp_dist_thresh:',as.character(bp_dist_thresh));
      print(msg)
      info(logger, msg)
      methylationSNMnorm_usr_prb_matrix_pointer = attach.big.matrix(methylationSNMnorm_usr_prb_matrix_desciption)
      #print(paste('time before run:', as.character(Sys.time())));
      clustering_methylation_probes(methylationSNMnorm_usr_prb_matrix_pointer, methylationCoord, corr_dist_thresh_name, corr_dist_thresh, bp_dist_thresh, chromosome_array_train, chrs_start_ind, chrs_end_ind, train_test_array[1]); # train data
      clustering_methylation_probes(methylationSNMnorm_usr_prb_matrix_pointer, methylationCoord, corr_dist_thresh_name, corr_dist_thresh, bp_dist_thresh, chromosome_array_test, chrs_start_ind, chrs_end_ind, train_test_array[2]); # test data
      
      #print(paste('time after run:', as.character(Sys.time())));
    }
  
}

msg = "Running has been finished. Reporting is starting..."
print(msg)
info(logger, msg)

if(report){
  
  #print(paste('time before report:', as.character(Sys.time())))
for(type in train_test_array){
  msg = paste("reporting for:",type,sep = "")
  print(msg);
  info(logger, msg);
  
  report_matrix_corrs = matrix(0, nrow= length(corr_dist_thresh_array), ncol = length(bp_dist_thresh_array))
  report_matrix_clusts_number = matrix(0, nrow = length(corr_dist_thresh_array), ncol = length(bp_dist_thresh_array));


  foreach (corr_dist_thresh_ind = 1:length(corr_dist_thresh_array),.packages=c("foreach","log4r","readr","bigmemory"))%:%
    foreach (bp_dist_thresh_ind = 1:length(bp_dist_thresh_array), .packages=c("foreach","log4r","readr","bigmemory"))%do%{
      # outer loop code
      methylationSNMnorm_usr_prb_matrix_pointer = attach.big.matrix(methylationSNMnorm_usr_prb_matrix_desciption)
      corr_dist_thresh_name = corr_dist_thresh_name_array[corr_dist_thresh_ind];
      corr_dist_thresh = corr_dist_thresh_array[corr_dist_thresh_ind]
      #
      
      
      
      bp_dist_thresh = bp_dist_thresh_array[bp_dist_thresh_ind];
      
      msg = paste('corr_tresh:',as.character(corr_dist_thresh), '&','bp_dist_thresh:',as.character(bp_dist_thresh));
      print(msg)
      info(logger, msg);
      
      which_clust_final = read.csv(file = paste("../../rosmapAD/data/aclust_assignments/which_clust_final_corr_dist_thresh_",corr_dist_thresh_name,"_bp_dist_thresh_",as.character(bp_dist_thresh),"_",type,".csv", sep=""),header = FALSE,sep = ',')
      which_clust_final = unlist(which_clust_final)

      number_of_clusters = length(unique(which_clust_final));

      clusters_start_ind = match(unique(which_clust_final),which_clust_final);
      clusters_end_ind = c(clusters_start_ind[2:length(clusters_start_ind)]-1,length(which_clust_final));

      for (i in 1:number_of_clusters){

            this_clust_start_ind = clusters_start_ind[i];
            this_clust_end_ind = clusters_end_ind[i];
            if(this_clust_end_ind - this_clust_start_ind == 0){
              sub_matrix = as.matrix(methylationSNMnorm_usr_prb_matrix_pointer[this_clust_start_ind:this_clust_end_ind,]);
            } else{
              sub_matrix = methylationSNMnorm_usr_prb_matrix_pointer[this_clust_start_ind:this_clust_end_ind,];
              sub_matrix = t(sub_matrix);
            }
            methSNMnorm_usr_prb_mtrx_for_clust = sub_matrix
            corr_matrix_for_clust = cor(methSNMnorm_usr_prb_mtrx_for_clust);
            report_matrix_corrs[corr_dist_thresh_ind, bp_dist_thresh_ind] = report_matrix_corrs[corr_dist_thresh_ind, bp_dist_thresh_ind] + mean(corr_matrix_for_clust)

      }

      report_matrix_clusts_number[corr_dist_thresh_ind, bp_dist_thresh_ind] =  number_of_clusters;
      report_matrix_corrs[corr_dist_thresh_ind, bp_dist_thresh_ind] = report_matrix_corrs[corr_dist_thresh_ind, bp_dist_thresh_ind] / number_of_clusters;

    
  }

  #print(paste('time after report:', as.character(Sys.time())))
  write.table(file = paste('report_matrix_corrs_',type,'.csv',sep = ""), report_matrix_corrs,sep = ',',row.names = FALSE,col.names = FALSE);
  write.table(file = paste('report_matrix_clusts_number_',type,'.csv', sep = ""), report_matrix_clusts_number,sep = ',', row.names = FALSE, col.names = FALSE);
  
  }

}

msg = "Reporting has been finished. Ploting is starting..."
print(msg)
info(logger, msg)


if(ploting){
  
  
  f_plotly <- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )
  x_plotly <- list(
    title = "BP dist threshold",
    titlefont = f_plotly
  )
  y_plotly <- list(
    title = "Corr threshold",
    titlefont = f_plotly
  )
  
  for (type in train_test_array){
    report_matrix_corrs = read.csv(file = paste('report_matrix_corrs_',type,'.csv',sep = ""), header = FALSE,sep = ',')
    report_matrix_corrs = as.matrix(report_matrix_corrs);
  
    report_matrix_clusts_number = read.csv(file = paste('report_matrix_clusts_number_',type,'.csv',sep = ""), header = FALSE,sep = ',')
    report_matrix_clusts_number = as.matrix(report_matrix_clusts_number);
  
  
    report_matrix_corrs_plotly = plot_ly(x = as.character(bp_dist_thresh_array), y = 1-corr_dist_thresh_array ,z = report_matrix_corrs, type ="heatmap") %>% layout(xaxis = x_plotly, yaxis = y_plotly)
    
    report_matrix_clusts_number_plotly = plot_ly(x = as.character(bp_dist_thresh_array), y = 1-corr_dist_thresh_array,z = report_matrix_clusts_number, type ="heatmap") %>% layout(xaxis = x_plotly, yaxis = y_plotly)
    
    export(report_matrix_corrs_plotly, file = paste("report_matrix_corrs_plotly_",type,".png",sep = ""))
    export(report_matrix_clusts_number_plotly, file = paste("report_matrix_clusts_number_plotly_",type,".png",sep = ""))
  }
}

msg = "Done."
print(msg)
info(logger, msg)

stopCluster(cl)