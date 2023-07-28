################################################
## Functions used in the process_data.r file ##
################################################



# Normalise data for varibility across samples
# @matrix: matrix of the input file
# @return: normlised matrix
normalise_data <- function(matrix){
  for(i in 1:nrow(matrix)){
    matrix[i,] = matrix[i,]/sum(matrix[i,], na.rm=T)
  }
  return(matrix)
}


# Parser for csv file format given by defined groups of samples
# @file: path of the csv file to read
# @sampe_groups: groups for samples
# @return: list
csv_parser <- function(file, sample_groups){
  raw <- read.csv(file=file, check.names=FALSE, sep=",", na.strings=c("","NA"))
  coln <- colnames(raw)
  sample_names <- as.character(coln[2:length(coln)])

  # ignore the first column
  new_raw <- raw[,2:ncol(raw)]

  colnames(new_raw) <- sample_names
  # extract lipid species 
  lipids <- as.vector(raw[,1])
  data <- list()
  gr <- unique(as.character(sample_groups$group))
  
  for(i in 1:length(gr)){
    samples <- as.character(sample_groups[which(sample_groups$group == gr[i]),]$sample)
    ids <- which(sample_names %in% samples)
    rawMt <- as.matrix(sapply(subset(new_raw, select = ids), as.numeric))
    rawMt <- t(rawMt)
    colnames(rawMt) <- lipids
    rownames(rawMt) <- sample_names[ids]
    data[[gr[i]]] <- normalise_data(rawMt)
  }
  return(data)
}


# Return lipid data from input file
# @data: input file (dataframe)
# @groups: groups of samples
# @return: dataframe
get_lipid_data <- function(data, groups){
  names <- colnames(data[[1]])
  pattern <- "^(FA|fa|Fa)(\\w+)?\\((\\d+):(\\d+)\\)"
  # extract only lipid species
  names <- names[!str_detect(names, pattern)]
  for(g in 1:length(groups)){
    data[[groups[g]]] <- data[[groups[g]]][,names]
    if(is.vector(data[[groups[g]]])){
      names(data[[groups[g]]]) = names
    }else{
      colnames(data[[groups[g]]]) = names
    }
  }
  return(data)
}


# Return fatty acid data from input file
# @data: input file (dataframe)
# @groups: groups of samples
# @return: dataframe
get_fa_data <- function(data, groups){
  names <- colnames(data[[1]])
  pattern <- "^(FA|fa|Fa)?\\((\\d+):(\\d+)\\)"
  names <- str_extract(names,pattern)
  names <- names[!is.na(names)]
  for(g in 1:length(groups)){
    data[[groups[g]]] = data[[groups[g]]][,names] 
    if(is.vector(data[[groups[g]]])){
      names(data[[groups[g]]]) = names
    }else{
      colnames(data[[groups[g]]]) = names
    }
  }
  return(data)
}


# Get a json object with the classification of the lipids in reactions
# e.g.: for reactions at the species level: "DG(30:0)" is classifies into "DG" into "Glycerolipids and Glycerophospholipids"
# @subset: reaction or pathway
# @level: class or species
# @return: json object
get_json_reactions_lp <- function(subset, level, group1, group2){
  # only take reactions which alter acyl chains
  if(subset == "reaction"){
      react_df <- lp_reaction
  }else{
    react_df <- pathway
    colnames(react_df)[3] <- "reaction"
  }

  # get lipids of a file per level/subset
  # by default: calculate "notpaired", so if lipids the file exist
  file_name <- paste0(out_file_path, paste0(c("lp", level, subset, group1, group2, "active", "notpaired.json"), collapse="_"))
  
  if(file.exists(file_name)){
    df <- fromJSON(file_name)$nodes[[1]]
    colnames(df) <- c("id", "lm_id", "label", "name", "shape")
    selected_nodes <- df[,3] # get label of selected lipids

    class <- unique(as.character(react_df$class))
    root <- list()

    tolower_react_keys <- tolower(react_keys)

    root <- list()
    for(i in 1:length(class)){
      cl <- class[i]
      parent <- list()
      valid_lp <- c()
      for(j in 1:nrow(react_df)){
        if(react_df$class[j] == cl){
          lp <- unique(strsplit(toString(react_df$reaction[j]), ",|->")[[1]])
          for(k in 1:length(lp)){
            if(!lp[k] %in% valid_lp){
              if(level == "species" || (level == "class" && lp[k] %in% selected_nodes)){
                valid_lp <- c(valid_lp, lp[k])
                
                if(level == "species"){ # look at the molecular species in the dataset
                  children <- list()
                  for(l in 1:nrow(node_species_df)){
                    node <- toString(node_species_df$label[l])
                    if(node %in% selected_nodes){
                      if(strsplit(node,"\\(")[[1]][1] == lp[k]){
                        children <- lappend(children, list(text = node))
                      }
                    }
                  }
                  if(length(children) > 0){
                    p <- list(text = lp[k], children = children)
                    parent <- lappend(parent, p)
                  }
                }else{
                  parent <- lappend(parent, list(text = lp[k]))
                }
              }
            }
          }
        }
      }
      if(length(parent) > 0){
          r <- list(text = cl, children = parent)
          root <- lappend(root, r)
      }
    }
    json_obj <- toJSON(root, pretty=T, auto_unbox=T)
    
  }else{ #if no lipids
    json_obj <- NULL
  }

  return(json_obj)
}
