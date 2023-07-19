###################################################################################################################
## * Retrieve information on the conditions (e.g. conditions with less than 2 samples are not valid)             ##
## * Create and save list of reactions, lipid data, FA data..                                                    ##
## * Generate the files to visualise the graphs/results tables (calculate z-scores, find "the most" pathways...) ##
###################################################################################################################


source("/lipidmaps/biopan/R/lib.r")
source("/lipidmaps/biopan/R/lib_process_data.r")
source("/lipidmaps/biopan/R/lib_pathway_analysis.r")


args <- commandArgs(TRUE)
# retrieve the values from the processing.php file
in_file <- args[1]
out_file_path <- args[2]
summary_json_file <- args[3]
samples <- strsplit(args[4], split = ',')[[1]]
groups <- strsplit(args[5], split = ',')[[1]]
p_value <- as.numeric(args[6])
is_paired <- as.logical(args[7])


# data frame with the samples and the condition name associated
sample_group_df <- data.frame(sample = samples, group = groups)

# get the groups condition
groups <- unique(groups)

# read data from file
data <- csv_parser(in_file, sample_group_df)

# read json summary file to extract data
sum_json <- fromJSON(summary_json_file, flatten=TRUE)
processed_species <- sum_json$pathway$processed
processed_class <- names(get_species_class(processed_species))

filtered_gr <- sample_group_df %>% group_by(group) %>% summarise(n_calls = sum(!is.na(sample)))

# groups conditions which contain less than two samples
filtered_small_gr <- filtered_gr[which(filtered_gr$n_calls <= 1),]
# groups conditions which contain at least two samples
filtered_gr <- filtered_gr[which(filtered_gr$n_calls > 1),]

small_groups <- groups[which(groups %in% as.character(filtered_small_gr$group))]

# remove groups condition which contain less that two samples
data <- data[which(groups %in% as.character(filtered_gr$group))]
groups <- groups[which(groups %in% as.character(filtered_gr$group))]

small_sample_group_df <- sample_group_df[which(sample_group_df$group %in% small_groups),]
sample_group_df <- sample_group_df[which(sample_group_df$group %in% groups),]

group_freq <- group_by(sample_group_df, group)
group_freq <- summarise(group_freq, length(group))
colnames(group_freq) <- c("group", "freq")

# check if there is at least 2 conditions with the same number of samples
dup <- duplicated(group_freq$freq) | duplicated(group_freq$freq, fromLast = TRUE)
is_dup <- FALSE
if(TRUE %in% dup){
  is_dup <- TRUE
  is_dup_group <- which(dup %in% TRUE)
}

# extract lipid data
lp_data <- get_lipid_data(data, groups)

# extract fatty acid data
fa_data <- get_fa_data(data, groups)

fa_names <- c()
fa_processed_species <- c()
fa_names <- colnames(fa_data[[1]])
fa_node_df <- data.frame()
if(length(fa_names)> 0){
  fa_processed_species <- fa_names[which(fa_names %in% processed_species)]
  if(length(fa_processed_species) > 0){
    fa_node_df <- get_species_df(fa_processed_species) # id | lm_id | label | name
  }
}

fa_reactions <- fa_reaction$reaction
node_species_df <- get_species_df(processed_species)


#### Create list of reaction groups and save to files ####
node_lp_df <- NULL
react_keys_file <- paste0(out_file_path,"reaction/reaction.RData")

if(file.exists(react_keys_file)){
  react_keys <- readRDS(react_keys_file)

  reaction_gr_list <- readRDS(paste0(out_file_path,"reaction_gr_list.RData"))
  for(i in 1:length(react_keys)){
    dt <- reaction_gr_list[[react_keys[i]]]
    write.csv(dt, file=paste0(out_file_path,"reaction/",tolower(react_keys[i]),".csv"))
  }
  node_lp_df <- get_lipid_node_data_frame(processed_class)
}


#### Find reactions/pathways - calculation - create files ####
# by default, p-value = 0.05
# if at least 2 conditions with the same number of samples and at least one molecular species
if(length(groups) > 1 && nrow(node_species_df) > 0 && is_dup){
  info <- create_files("process", lp_data, fa_data, node_lp_df, node_species_df, fa_node_df, groups, is_paired, p_value, react_keys, out_file_path)

  # get json object with the classification of the lipids in reactions
  json_lp_class_react <- get_json_reactions_lp("reaction", "class", groups[is_dup_group[1]], groups[is_dup_group[2]])
  json_lp_species_react <- get_json_reactions_lp("reaction", "species", groups[is_dup_group[1]], groups[is_dup_group[2]])

  json_lp_class_pathway <- get_json_reactions_lp("pathway", "class", groups[is_dup_group[1]], groups[is_dup_group[2]])
  json_lp_species_pathway<- get_json_reactions_lp("pathway", "species", groups[is_dup_group[1]], groups[is_dup_group[2]])

  if(!is.null(json_lp_class_react)){
    out_file <- paste0(out_file_path,"lp_class_reaction.json")
    save_json_to_file(json_lp_class_react, out_file)
  }
  if(!is.null(json_lp_species_react)){
    out_file <- paste0(out_file_path,"lp_species_reaction.json")
    save_json_to_file(json_lp_species_react, out_file)
  }
  if(!is.null(json_lp_class_pathway)){
    out_file <- paste0(out_file_path,"lp_class_pathway.json")
    save_json_to_file(json_lp_class_pathway, out_file)
  }
  if(!is.null(json_lp_species_pathway)){
    out_file <- paste0(out_file_path,"lp_species_pathway.json")
    save_json_to_file(json_lp_species_pathway, out_file)
  }
}else{
  info <- NULL
}


#### Summary data ####
if(nrow(small_sample_group_df) > 0){
  small_group_freq <- group_by(small_sample_group_df, group)
  small_group_freq <- summarise(small_group_freq, length(group))
  colnames(small_group_freq) <- c("group", "freq")
  notvalid <- list(groups = small_group_freq$group, freqs = small_group_freq$freq)
}else{
  notvalid <- data.frame(groups = NULL, freqs = NULL)
}

if(length(unique(group_freq$freq)) > 1){ # if all groups don't have the same number of samples
  group_sample <- list(equal = FALSE)
}else{
  group_sample <- list(equal = TRUE)
}

valid <- list(groups = group_freq$group, freqs = group_freq$freq)
valid_reaction <- list(lp = info$valid_lp_reaction, fa = info$valid_fa_reaction)
valid_subset_reaction <- list(reaction = info$subset_reaction, pathway = info$subset_pathway)

msg <- list(valid = valid, notvalid = notvalid, reaction = valid_reaction, subset = valid_subset_reaction, group_sample = group_sample, is_dup = is_dup)
msg_json <- toJSON(msg)
out_file <- paste0(out_file_path,"msg2.json")
save_json_to_file(msg_json, out_file)


#### Save data frame for later use (calculate results for other p_values) ####
saveRDS(lp_data, paste0(out_file_path,"lp_data.RData"))

if(!is.null(node_lp_df)){
  saveRDS(node_lp_df, paste0(out_file_path,"node_lp_df.RData"))
}

if(length(fa_data) > 0){
  saveRDS(fa_data, paste0(out_file_path,"fa_data.RData"))
  if(nrow(fa_node_df) > 0){
    saveRDS(fa_node_df, paste0(out_file_path,"fa_node_df.RData"))
  }
}

if(nrow(node_species_df) > 0){
  saveRDS(node_species_df, paste0(out_file_path,"node_species_df.RData"))
}
