##########################################################################################################
## Generate the files to visualise the graphs/results tables with a different p-value and/or paire-data ##
##########################################################################################################


source("/lipidmaps/biopan/R/lib.r")
source("/lipidmaps/biopan/R/lib_pathway_analysis.r")


args <- commandArgs(TRUE)
# retrieve the values from the compute_pathway.php file
out_file_path <- args[1]
pw_level <- args[2]
p_value <- as.numeric(args[3])
is_paired <- as.logical(args[4])


# to avoid problems converting the p_value into exponential value
write_p_value <- as.character(format(p_value, scientific=F))

paired <- write_isPaired(is_paired)

# get the number of file names containing the combination [p_value + paired_data]
# if the [p_value + paired_data] is found in at least one file -> the pathways have already been computed and we don't need to compute it again
find_files1 <- length(grep(paste0(write_p_value, "_", paired),list.files(path = out_file_path)))

# if no file exists
if(find_files1 == 0){
  lp_data <- NULL
  fa_data <- NULL
  node_lp_df <- NULL
  node_species_df <- NULL
  fa_node_df <- NULL

  # retrieve saved data
  node_lp_df_path <- paste0(out_file_path,"node_lp_df.RData")
  fa_node_df_path <- paste0(out_file_path,"fa_node_df.RData")

  lp_data_path <- paste0(out_file_path,"lp_data.RData")
  fa_data_path <- paste0(out_file_path,"fa_data.RData")

  node_species_df_path <- paste0(out_file_path,"node_species_df.RData")

  if(file.exists(node_lp_df_path)){
    node_lp_df <- readRDS(node_lp_df_path)
  }
  if(file.exists(fa_node_df_path)){
    fa_node_df <- readRDS(fa_node_df_path)
  }

  if(file.exists(lp_data_path)){
    lp_data <- readRDS(lp_data_path)
  }
  if(file.exists(fa_data_path)){
    fa_data <- readRDS(fa_data_path)
  }

  if(file.exists(node_species_df_path)){
    node_species_df <- readRDS(node_species_df_path)
  }

  # read json msg file to get group data
  msg_file <- paste0(out_file_path,'msg2.json')
  msg_json <- fromJSON(msg_file, flatten=TRUE)
  groups <- msg_json$valid$groups


  # if the [paired_data] value is found in at least one file -> launch all the functions: process
  find_files2 <- length(grep(paste0("_", paired),list.files(path = out_file_path)))
  if(find_files2 == 0){
    react_keys_file <- paste0(out_file_path,"reaction/reaction.RData")
    if(file.exists(react_keys_file)){
      react_keys <- readRDS(react_keys_file)
    }else{
      react_keys <- NULL
    }
    info <- create_files("process", lp_data, fa_data, node_lp_df, node_species_df, fa_node_df, groups, is_paired, p_value, react_keys, out_file_path)

  }else{  # only the p_value files does not exist -> use existing files and do not recaculate scores 
    info <- create_files("compute", lp_data, fa_data, node_lp_df, node_species_df, fa_node_df, groups, is_paired, p_value, react_keys = NULL, out_file_path)
  }
}
