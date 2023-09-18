# set up R
packages_needed <- c("dplyr", "tidyverse", "ggplot2", "readxl", "RColorBrewer", "ggpubr", "stats")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

###########################################################################################################################

# provide files, etc
#to run on the cluster
  #I have not tested this part, but I'm happy to help debug if you try and it fails.

#args = commandArgs(trailingOnly=TRUE)
#BASEDIR <- args[1]
#PREFIX <- args[2]
#METADATAFILE <- args[3]
#SHEET <- args[4]
#FEATURE_NAME <- args[5]

#to run locally
BASEDIR <- "C:/Users/laura.timm/Desktop/Sebastes_ruberrimus/" #the directory where the metadata file lives
PREFIX <- "SRUB-SMIN" #the prefix for the lcWGS run
METADATAFILE <- paste0(BASEDIR, "yelloweye_metadata.xlsx") #metadatafile name (excel file)
SHEET <- "sheet1" #sheet containing the pertinent metadata
FEATURE_NAME <- "loc" #column name (in the metadata file) that contains the categorical variable you want to compare depth across (batch, region, sex, etc)
DEPTH_LOWER_BOUND <- 1 #this is the default, but feel free to tune as desired

###########################################################################################################################

#define functions
  
#get_factors <- function(df, colnm) {
#  outlist <- unique(df[[colnm]])
#  return(outlist)
#}

depths_t_test <- function(df, colnm, f1, f2) {
  first_factor_depths <- df[df[[colnm]] == f1,
                            which( colnames(df) == "mean_depth" )]
  second_factor_depths <- df[df[[colnm]] == f2, 
                             which( colnames(df) == "mean_depth" )]
  if (length(first_factor_depths$mean_depth) < 5) {
    print(paste0("WARNING: ", f1, " is represented by fewer than five individuals. Consider removing these from the depths analysis and adding ", f1, " to the downsample list (factors_to_downsample) manually."))
  }
  if (length(second_factor_depths$mean_depth) < 5) {
    print(paste0("WARNING: ", f2, " is represented by fewer than five individuals. Consider removing these from the depths analysis and adding ", f2, " to the downsample list (factors_to_downsample) manually."))
  }
  t_out <- t.test(x = as.numeric(first_factor_depths$mean_depth), y = as.numeric(second_factor_depths$mean_depth), var.equal = FALSE)
  outlist <- list(pval = t_out$p.value, f1_depths = as.numeric(first_factor_depths$mean_depth), f2_depths = as.numeric(second_factor_depths$mean_depth))
  return(outlist)
}

###########################################################################################################################

# prepare data

full_metadata <- read_xlsx(METADATAFILE, sheet = SHEET)

features_df <- full_metadata[c("ABLG", FEATURE_NAME)]

just_depths <- read.csv(paste0(BASEDIR, PREFIX, "_depths.csv"), header = FALSE, sep = '\t')
colnames(just_depths) <- c("ABLG", "mean_depth")

depths_df <- left_join(features_df, just_depths, by = "ABLG")
colnames(depths_df) <- c("ABLG", "feature", "mean_depth")

###########################################################################################################################

# plot depth distribution, colored by feature of interest

depths_plot <- ggplot(data = depths_df, aes(x = reorder(ABLG, -mean_depth), y = mean_depth, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("individual") +
  ylab("mean depth") +
  geom_hline(yintercept = DEPTH_LOWER_BOUND) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

depths_plot

ggsave(paste0(BASEDIR, PREFIX, "-", FEATURE_NAME, "_mean_depths.jpg"), width = 8, height = 5, units = "in", dpi = 300)


###########################################################################################################################

# compare depth distributions (with a series of t-tests)

whitelist_df <- depths_df[depths_df$mean_depth >= DEPTH_LOWER_BOUND,]
#for shortraker
#whitelist_df <- whitelist_df %>%
#  filter(feature == "Chignik" | feature == "Chirikof")

factors_list <- unique(whitelist_df[["feature"]])

## get a list of factors (batch1, batch2, etc) that require downsampling
downsampled_factors <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(downsampled_factors) <- c("lo_factor", "hi_factor", "targeted_mean", "targeted_sd")

for (i in 1:length(factors_list)) {
  for (j in 2:length(factors_list)) {
    if (i < j) {
      first_factor <- as.character(factors_list[i])
      second_factor <- as.character(factors_list[j])
      results <- depths_t_test(whitelist_df, "feature", first_factor, second_factor)
      if (results$pval < 0.05) {
        if (mean(results$f1_depths) < mean(results$f2_depths)) {
          lo_mean <- mean(results$f1_depths)
          lo_sd <- sd(results$f1_depths)
          target_factor <- first_factor
          downsample_factor <- second_factor
        } else {
          lo_mean <- mean(results$f2_depths)
          lo_sd <- sd(results$f2_depths)
          target_factor <- second_factor
          downsample_factor <- first_factor
        }
        downsampled_factors <- rbind(downsampled_factors, 
                                     list("lo_factor" = as.character(target_factor),
                                          "hi_factor" = as.character(downsample_factor), 
                                          "targeted_mean" = as.double(lo_mean),
                                          "targeted_sd" = as.double(lo_sd)))
      }
    }
  }
}
#did it work?
downsampled_factors

#draw new depth values, if needed
if (nrow(downsampled_factors) > 0) {
  target_mean <- NULL
  target_sd <- NULL
  to_downsample <- c()
  
  ## determine the downsample targets (mean and sd)
  for (dsf in 1:nrow(downsampled_factors)) {
    if (downsampled_factors[[dsf, "targeted_mean"]] == min(downsampled_factors$targeted_mean)) {
      trendsetter <- downsampled_factors[[dsf, "lo_factor"]]
      target_mean <- downsampled_factors[[dsf, "targeted_mean"]]
      target_sd <- downsampled_factors[[dsf, "targeted_sd"]]
      needs_downsampling <- as.character(downsampled_factors[[dsf, "hi_factor"]])
      if (!(needs_downsampling %in% to_downsample)) {
        to_downsample <- c(to_downsample, needs_downsampling)
      }
      
    }
  }
  #did it work?
  trendsetter
  target_mean
  target_sd
  to_downsample
  
  downsample_df <- NULL
  
  for (high_depth_factor in to_downsample) {
    depths_pval <- 0 #to get the ball rolling
    downsampling_round <- 1
    
    while (depths_pval < 0.05) {
      print(paste0("downsampling round ", downsampling_round, " for ", high_depth_factor, " has resulted in..."))
      check_df <- whitelist_df[whitelist_df$feature == trendsetter,]
      check_df$target_prop <- 1
      high_depth_factor_df <- whitelist_df[whitelist_df$feature == high_depth_factor,]
      for (sample in 1:nrow(high_depth_factor_df)) {
        id <- high_depth_factor_df[[sample, "ABLG"]]
        current_depth <- high_depth_factor_df[[sample, "mean_depth"]]
        target_proportion <- 1
        target_depth <- 0.5
        while (target_depth <= DEPTH_LOWER_BOUND) {
          target_depth <- rnorm(1, mean = target_mean, sd = target_sd)
          }
        if (target_depth < current_depth) {
          target_proportion <- 42 + (target_depth / current_depth)
        } else {
            target_depth <- current_depth
        }
        check_df <- rbind(check_df,
                          as_tibble(list("ABLG" = as.character(id), 
                                         "mean_depth" = as.double(target_depth), 
                                         "feature" = as.character(high_depth_factor), 
                                         "target_prop" = as.double(target_proportion))))
      }
      
      depths_pval <- depths_t_test(check_df, "feature", trendsetter, high_depth_factor)$pval
      
      if (depths_pval < 0.05) {
        print("failure.")
        downsampling_round <- downsampling_round + 1
      }
    }
    print("success!")
    downsample_df <- rbind(downsample_df, check_df[check_df$feature == high_depth_factor,])
  }
}

#this step may have some unexpected behavior if some sample sizes are very small; check the output to make sure things look right before running the downsample step

###########################################################################################################################

#plot new depth values, if needed
if (nrow(downsampled_factors) > 0) {
  plotting_df <- downsample_df[,1:3]
  for (f in factors_list) {
    if (!(f %in% to_downsample)) {
      plotting_df <- rbind(plotting_df, whitelist_df[whitelist_df$feature == f,])
    }
  }
  
  ds_depths_plot <- ggplot(data = plotting_df, aes(x = reorder(ABLG, -mean_depth), y = mean_depth, fill = feature)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab("individual") +
    ylab("mean depth") +
    geom_hline(yintercept = 1) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(paste0(BASEDIR, PREFIX, "-", FEATURE_NAME, "_downsampled_mean_depths.jpg"), width = 8, height = 5, units = "in", dpi = 300)
  
  ggarrange(depths_plot, ds_depths_plot, ncol = 2, nrow = 1, common.legend = TRUE)
  
}


###########################################################################################################################

# write out subsetBAMS_array_input

if (exists("downsample_df")) {
  outfile <- file(paste0(BASEDIR, PREFIX, "-", FEATURE_NAME, "_downsample-bamsARRAY_inputTEST.txt"), "wb")
  iterator <- 1
  for (sample in 1:nrow(downsample_df)) {
    if (downsample_df[sample, "target_prop", 1] != 1) {
      outline <- paste0(as.character(iterator), ":", downsample_df[sample, "ABLG"], ":", as.character(round(as.numeric(downsample_df[sample, "target_prop"]), 3)))
      write(outline, file = outfile, append = TRUE)
      iterator <- iterator + 1
    }
  }
  
  close(outfile)
  
  #the file that was written is uploaded to sedna and used as the input array file for lcWGSpipeline_opt-downsample.py.
  #Be sure -a flag in opt-downsample matches the name of the file written here
  #unless you change it, the filename is "PREFIX-FEATURENAME_downsample-bamsARRAY_input.txt":
  print(paste0(PREFIX, "-", FEATURE_NAME, "_downsample-bamsARRAY_input.txt"))
}
