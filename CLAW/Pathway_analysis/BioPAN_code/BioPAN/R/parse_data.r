###################################################################################################################################
## * Reformat the loaded dataset (remove LipidLynxX unused columns, replace NAs values, identify duplicates and sum them,... )   ##
## * Get a list of reactions at the lipid molecular species level for each reaction at the subclass level                        ##
##    e.g. reactions: PC(30:0) -> PA(30:0), PC(32:0) -> PA(32:0) and PC(32:1) -> PA(32:1) for the reaction PC -> PA              ##
## * Categories each molecular species (unrecognised/processed/unprocessed)                                                      ##
###################################################################################################################################


# source("/lipidmaps/biopan/R/lib.r")
# source("/lipidmaps/biopan/R/lib_parse_data.r")
# 

source("R/lib.r")
source("R/lib_parse_data.r")

# 
# args <- commandArgs(TRUE)
# # retrieve the values from the parse_data.php file
# input_path <- args[1]
# input_path_clean <- args[2]
# output_path <- args[3]
# lipidlynxx <- args[4]
# output_lipidlynxx <- args[5]



# Define the file paths and values
input_path <- "input/input.csv"
input_path_clean <- "input/input_clean.csv" # Modify as needed
output_path <- "output/" # Modify as needed
lipidlynxx <- FALSE
output_lipidlynxx <- "output/output_lipidlynxx.csv" # Modify as needed

#### Reformat the file loaded by the user ####

# if the user import his own file LipidLynxX was launched
if(lipidlynxx == "yes"){
  # read data from file: output of LipidLynxX for the datasets processed by it
  if(file.exists(output_lipidlynxx)){
    data <- read.csv(file=output_lipidlynxx, check.names=FALSE, sep=",", na.strings=c("","NA"))

    
    # store lipid molecular species classified as UNPROCESSED and delete them from dataset file
    # these lipids will be categories later (unprocessed or unrecognised)
    to_categories <- data[data[,2] == "UNPROCESSED",3]
    if(length(to_categories) > 0){
      un_processed_lipidlynxx <- which(data[,2] == "UNPROCESSED")
      for(i in 1:length(un_processed_lipidlynxx)){
        data[un_processed_lipidlynxx[i],2] <- data[un_processed_lipidlynxx[i],3]
      }
    }
    
    colnames(data)[3] <- "Dataset"
    colnames(data)[2] <- "Sample"

    conversion <- data[c("Sample", "Dataset")]

    # remove unused columns to match with the BioPAN format file
    data <- data[,-c(1,3), drop=FALSE]

    raw <- data
    # remove empty rows
    data <- data[rowSums(is.na(data)) != ncol(data),]
    empty_rows <- FALSE
    if(!identical(raw,data)){
      empty_rows <- TRUE
    }

    raw <- data
    # replace NA values by 0
    data[is.na(data)] <- 0
    na_values <- FALSE
    if(!identical(raw,data)){
      na_values <- TRUE
    }

    # identify duplicated lipids
    nb_occur <- data.frame(table(data[,1]))
    nb_occur[nb_occur$Freq > 1,]
    duplicate_val <- data[data[,1] %in% nb_occur$Var1[nb_occur$Freq > 1],]
    dup_lp <- NULL
    if(nrow(duplicate_val) > 0){
      dup_lp <- unique(duplicate_val[,1])
      
      # sum the values of the lines with duplicate lipids (update the first line containing the lipid and delete the others)
      for(i in 1:length(dup_lp)){
        rows <- data[data[,1] == dup_lp[i], ]
        rowindex <- rownames(rows[2:ncol(rows)])
        for(j in 2:ncol(data[rowindex[1],])){
          data[rowindex[1],j] <- sum(rows[1:nrow(rows),j])
        }
        rows_delete <- paste0("c(", paste0(as.numeric(rowindex[2:length(rowindex)]), collapse = ","), ")")
        data <- data[!rownames(data) %in% rowindex[2:length(rowindex)],]
      }
    }

    msg <- list(na_values = na_values, empty_rows = empty_rows, dup_lp = dup_lp, error = FALSE)
  
  }else{ # issue with LipidLynxX
    msg <- list(na_values = FALSE, empty_rows = FALSE, dup_lp = NULL, error = TRUE)
  }



}else{ # if the user try BioPAN with a demo file
  data <- read.csv(file=input_path, check.names=FALSE, sep=",", na.strings=c("","NA"))
  to_categories <- c()
  msg <- list(na_values = FALSE, empty_rows = FALSE, dup_lp = NULL, error = FALSE)
}


# Save updates on the input file
msg_json <- toJSON(msg)
out_file <- paste0(output_path,"msg1.json")
save_json_to_file(msg_json, out_file)


# save the file without empty rows/values and duplicates lipids (replace the previous)
write.csv(data, file=input_path_clean, row.names=FALSE)



#### Extract data and create list of reaction groups (functions in lib_parse_data.r) ####

coln <- colnames(data)
samples <- as.character(coln[2:length(coln)])

# get the names of the samples
sample_groups <- get_sample_groups(samples)

# extract the lipid molecular species 
species <- as.vector(data[,1])
species_total = length(species)

# get mapper of all molecular species
#mapper <- get_species_standard_mapper(species)

# get standard species names
if(length(to_categories) > 0){
  species <- species[-which(species %in% to_categories)]
}

# extract fatty acyls [FA] lipid molecular species
fa_gr <- get_fas(species)

# list of lipid subclass containing the associated molecular species
species_groups <- get_species_class(species)


processed_species <- c()
processed_gr <- c()
fa_react_list <- c()
fa_processed_gr <- c()

if(length(fa_gr) > 0){
  fa_react_list <- get_reaction_list(fa_reaction)  

  valid_fa <- names(fa_react_list)
  # get all processed fatty acids
  fa_processed_gr <- get_processed_groups(fa_gr, valid_fa, fa_react_list)
  # add FaCoA species
  fa_processed_gr <- c(fa_processed_gr, fa_gr[which(grepl('facoa', tolower(fa_gr)))])

  # store all processed fatty acids molecular species
  processed_species <- c(processed_species, fa_processed_gr)
}

# extract lipid molecular species
lipids <- get_lipid_species(species)
lp_react_list <- c()
reaction_gr_list <- c()

# store temporaly the lipid processed molecular species (no FA)
processed_species_lp <- c()

# if at least one lipid group (no FA/FACoA)
if(length(lipids) > 0){
  lp_react_list <- get_reaction_list(lp_reaction)

  # get list of groups reactions
  reaction_gr_list <- get_group_reaction_list(lipids)
  if(length(reaction_gr_list) > 0){
    react_keys <- names(reaction_gr_list)
    saveRDS(react_keys, paste0(output_path,"reaction/reaction.RData"))

    valid_lp_class <- names(lp_react_list)

    # extract lipid subclasses from dataset
    gr <- names(species_groups)


    #### define the selected and non-selected molecular species for each reaction ####

    # get all processed lipid subclasses
    processed_gr <- get_processed_groups(gr, valid_lp_class, lp_react_list)
    
    # store all processed lipids (no fatty acids)
    for(i in 1:length(react_keys)){
      keys <- react_keys[i]
      tokens <- strsplit(keys, ",")[[1]]
      re <- tokens[1]
      pro <- tokens[2]
      selected_reaction <- reaction_gr_list[[keys]]
      
      selected_re <- gsub(" ", "",toString(selected_reaction$reactant))
      selected_pro <- gsub(" ", "",toString(selected_reaction$product))
      selected_re_arr <- strsplit(selected_re,",")[[1]]
      selected_pro_arr <- strsplit(selected_pro,",")[[1]]

      # remove duplicate lipids due to sphingoid base lipids (SPB/SPBP/dhSPB/dhSPBP)
      if(length(unique(selected_re_arr)) == 1){
        selected_re <- selected_re_arr[1]
      }
      if(length(unique(selected_pro_arr)) == 1){
        selected_pro <- selected_pro_arr[1]
      }

      processed_species_lp <- c(processed_species_lp, selected_re_arr, selected_pro_arr)

      non_selected_re <- species_groups[[re]][which(!species_groups[[re]] %in% selected_re_arr)]
      non_selected_pro <- species_groups[[pro]][which(!species_groups[[pro]] %in% selected_pro_arr)]

      non_selected_re_str <- c()
      non_selected_pro_str <- c()
      if(length(non_selected_re) > 0){
        non_selected_re_str <- paste0(non_selected_re, collapse=",")
      }
      if(length(non_selected_pro) > 0){
        non_selected_pro_str <- paste0(non_selected_pro, collapse=",")
      }

      data_list <- list(non_selected_re = non_selected_re_str, selected_re = selected_re, non_selected_pro = non_selected_pro_str, selected_pro = selected_pro)
      json_obj <- toJSON(data_list, pretty=T, auto_unbox=T)
      # save to file
      out_file <- paste0(output_path, tolower(keys),".json")
      save_json_to_file(json_obj, out_file)
    }
  }
}

# remove duplicate lipid molecular species (one species can be selected from more than one reaction)
if(length(processed_species_lp) > 1){
  processed_species <- c(processed_species, unique(processed_species_lp))
}

# get undefined molecular species (= subclass or molecular species (for FA) that are not recognised by BioPAN)
undef_species <- c()
undef_species <- get_undefined_species(species)

# species that are not processed - based on 2 others catergory already define
# correspond to lipid that are not processed and not undefined
if(!is.null(c(processed_species, undef_species))){
  un_processed_species <- species[-which(species %in% c(processed_species,undef_species))]
}else{
  un_processed_species <- species
}


# mapp the nomenclature use by BioPAN with the nomenclature from the input dataset
processed_species_dataset <- list()
undef_species_dataset <- list()
un_processed_species_dataset <- list()
if(lipidlynxx == "yes"){
  processed_species_dataset <- get_biopan_dataset_nom(processed_species, processed_species_dataset, conversion$Dataset, conversion$Sample)
  undef_species_dataset <- get_biopan_dataset_nom(undef_species, undef_species_dataset, conversion$Dataset, conversion$Sample)
  un_processed_species_dataset <- get_biopan_dataset_nom(un_processed_species, un_processed_species_dataset, conversion$Dataset, conversion$Sample)
}


# add LipidLynxX unprocessed molecular species to the correct catergory (unrecognised or unprocessed)
if(length(to_categories) > 0){
  to_categories_mapper <- get_species_standard_mapper(to_categories)
  to_categories_undef <- get_undefined_species(to_categories_mapper$new_name)
  to_categories_undef <- to_categories_mapper$old_name[which(to_categories_mapper$new_name %in% to_categories_undef)] #retrieve the nomenclature in the dataset because LipidLynxX didn't processed this lipid

  if(!is.null(to_categories_undef)){
    undef_species <- c(undef_species, to_categories_undef)
    to_categories_unprocessed <- setdiff(to_categories_mapper$old_name, to_categories_undef)
    un_processed_species <- c(un_processed_species, to_categories_unprocessed)

    # mapp the molecular species output UNPROCESSED by LipidLynxX
    undef_species_dataset <- get_biopan_dataset_nom_unprocessed(to_categories_undef, undef_species_dataset)
    un_processed_species_dataset <- get_biopan_dataset_nom_unprocessed(to_categories_unprocessed, un_processed_species_dataset)

  }else{
    un_processed_species <- c(un_processed_species, to_categories_mapper$old_name)    
    un_processed_species_dataset <- get_biopan_dataset_nom_unprocessed(to_categories_mapper$old_name, un_processed_species_dataset)
  }
}

summary <- list(total = species_total, groups = sample_groups, lipidlynxx = lipidlynxx, undef = undef_species, undef_dataset = undef_species_dataset, processed_dataset = processed_species_dataset, unprocessed_dataset = un_processed_species_dataset, pathway = list(name = 'Pathway anlysis', processed = processed_species, unprocessed = un_processed_species))
json_obj <- toJSON(summary, pretty=T, auto_unbox=T)
out_file <- paste0(output_path,"summary.json")
save_json_to_file(json_obj, out_file)


# save lipid reactions and pathways
if(length(reaction_gr_list) > 0){
  saveRDS(reaction_gr_list, paste0(output_path,"reaction_gr_list.RData"))
}

# save processed FA 
if(length(fa_gr) > 0){
  saveRDS(fa_processed_gr, paste0(output_path,"fa.RData"))
}
