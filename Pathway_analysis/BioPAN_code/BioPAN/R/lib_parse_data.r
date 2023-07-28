################################################
## Functions used in the parse_data.r file    ##
################################################



# Get the longest common string starting from potition 0
# @word1: string
# @word2: string
# @return: string
get_longest_common_string <- function(word1,word2){
  n <- min(nchar(word1), nchar(word2))  # length of the shorter word
  # two vectors of characters of the same length n
  c1 <- strsplit(word1, "", fixed = TRUE)[[1]][1:n]
  c2 <- strsplit(word2, "", fixed = TRUE)[[1]][1:n]

  m <- as.logical(cumprod(c1 == c2))  # a vector that is TRUE as long as the characters match  
  return(paste(c1[m], collapse = ""))
}


# Return the longest string
# @str_v: string
get_longest_string <- function(str_v){
  l <- length(str_v)
  if(l == 1){
    return(str_v[1])
  }
  max <- 1
  pos <- 1
  for(i in 1:l){
    l1 <- nchar(str_v[i])
    if(l1 >= max){
      max <- l1
      pos <- i
    }
  }
  return(str_v[pos])
}


# Extract sample groups from the sample names
# @sample_names: sample names
# @return: list
get_sample_groups <- function(sample_names) {
  s1_v <- c()
  s2_v <- c()
  com_str_v <- c()
  l <- length(sample_names)
  for(i in 1:(l-1)){
    for(j in (i+1):l){
      com_str <- get_longest_common_string(sample_names[i], sample_names[j])
      s1_v <- c(s1_v, sample_names[i])
      s2_v <- c(s2_v, sample_names[j])
      com_str_v <- c(com_str_v, com_str)
    }
  }
    
  res_df <- data.frame(s1 = s1_v, s2 = s2_v, com_str = com_str_v)
  # remove empty common string
  res_df <- res_df[res_df$com_str != "",]
  #sort
  res_df <- dplyr::arrange(res_df, com_str)
  
  com_str <- unique(as.character(res_df$com_str))
  gr <- list()
  for(i in 1:length(com_str)){
    tmp <- c()
    df <- res_df[res_df$com_str == com_str[i],]
    for(j in 1: nrow(df)){
      dt <- df[j,]
      tmp <- c(tmp, union(as.character(dt$s1), as.character(dt$s2)))
    }
    gr[[com_str[i]]] <- unique(tmp)
  }

  gr_names <- names(gr)
  sample_gr <- list()
  for(i in 1:length(gr_names)){
    s <- gr[[gr_names[i]]]
    for(j in 1:length(s)){
      if(s[j] %in% names(sample_gr)){
        sample_gr[[s[j]]] <- c(sample_gr[[s[j]]], gr_names[i])
      }else{
        sample_gr[[s[j]]] <- gr_names[i]
      } 
    }
  }
  
  for(i in 1:length(sample_names)){
    g <- sample_gr[[sample_names[i]]]
    sample_gr[[as.character(sample_names[i])]] <- get_longest_string(g)
  }

  s1 <- as.data.frame(sample_gr, check.names=FALSE)
  s1 <- s1[order(colnames(s1))]
  s1 <- as.list(s1)
  return(s1)
}


# Get a list of reactions (lipid or FA) for each class of lipid
# @react_df: dataframe of reactions
# @return: list
get_reaction_list <- function(react_df){
  react_list = list()
  for(i in 1:nrow(react_df)){
    react <- strsplit(react_df[i,]$reaction, ",")[[1]]
    re <- react[1]
    pro <- react[2]
    gr_names <- names(react_list)
    if(is.null(gr_names)){
        react_list[[re]] <- pro
    }else if(re %in% gr_names){
        react_list[[re]] <- c(react_list[[re]], pro)
    }else{
        react_list[[re]] <- pro
    }
  }
  return(react_list)
}


# Return fatty acids species
# @species: all the species names in the input file
# @return: vector
get_fas <- function(species){
  pattern <- "^(FA|fa|Fa)(\\w+)?\\((\\d+):(\\d+)\\)"
  names <- str_extract(species, pattern)
  names <- names[!is.na(names)]
  return(names)
}


# Return lipids species (without FA)
# @species: all the species names in the input file
# @return: vector
get_lipid_species <- function(species){
  pattern <- "^(FA|fa|Fa)(\\w+)?\\((\\d+):(\\d+)\\)"
  names <- species[!str_detect(species, pattern)]
  return(names)
}




#### Get a list of reactions at the lipid molecular species level for each reaction at the subclass level ####

header <- c("reactant_group","product_group")

# Create a matrix to store the reactions
create_matrix <- function(){
  mat <- matrix(,,ncol=2)
  colnames(mat) <- header
  mat <- mat[-1,]
  return(mat)
}


# Get the acyls chains of FA and FACoA in the dataset
fa_arr <- c()
facoa <- c()
facoa_acyl <- c()
fa_acyl <- c()
fa <- c()
get_fa_coa_acyl_chain <- function(){
  fa_arr <- fa_processed_gr
  if(length(fa_arr) > 0){
    facoa <- fa_arr[grepl("facoa", tolower(fa_arr))]
    if(length(facoa) > 0){
      fa <- fa_arr[!fa_arr %in% facoa]
      facoa_acyl <- get_acyl_chain(facoa)
    }
    if(length(fa) > 0){
      fa_acyl <- get_acyl_chain(fa)
    }
  }
  res <- list(fa_acyl = fa_acyl, facoa_acyl = facoa_acyl)
  return(res)
}


##For reactions involving sphingoid bases lipids [SP01] ##
# Get a list of reactions at the lipid molecular species level
# @react: reactions
# @res: list of reactions at the lipid molecular species level to complete
# @key: reaction
# @re: reactant
# @pro: product
# @r_match: reactant molecular species
# @p_match: product molecular species
# @fa_acyl: FA acyl chains
# @facoa_acyl: FACoA acyl chains
get_reaction_sphingo <- function(react, res, key, re, pro, r_match, p_match, fa_acyl, facoa_acyl){
  sphingoid_bases <- c("SPB", "SPBP", "dhSPB", "dhSPBP")

  # if reactant and product are sphingoid bases
  if(re %in% sphingoid_bases && pro %in% sphingoid_bases){
    res[[key]] <- list(reactant = re, product = pro)
  
  }else{
    acyl_chain <- lp_reaction[lp_reaction$reaction == react, "compound_require"]
    if(acyl_chain == "fa"){
      # a reactant molecular species can be converted into a product (SPB/SPBP/dhSPB/dhSPBP) if there is a FA with the same acyl chain in the dataset
      re_acyl_chain <- data.frame("species"=r_match , "chain"=get_acyl_chain(r_match))
      r_match <- re_acyl_chain$species[which(re_acyl_chain$chain %in% fa_acyl)]
      p_match <- c(rep(pro, length(r_match)))
    }else{
      pro_acyl_chain <- data.frame("species"=p_match , "chain"=get_acyl_chain(p_match))
      p_match <- pro_acyl_chain$species[which(pro_acyl_chain$chain %in% facoa_acyl)]
      r_match <- c(rep(re, length(p_match)))
    }

    if(!identical(r_match, character(0))){ #if there is at least a reaction between 2 molecular species
      mat <- create_matrix()
      for(i in 1:length(r_match)){
        mat <- rbind(mat, c(r_match[i], p_match[i]))
      }
      res[[key]] <- as.data.frame(mat)
    }
  }
  return(res)
}


## For reactions that do not involve FA or FACOA ##
# Get a list of reactions at the lipid molecular species level for the reactions that does not involved FA/FACoA
# @reactions: list of reactions
# @res: list of reactions at the lipid molecular species level to complete
# @keys: reactions already analysed
# @key: reaction
# @re: reactant
# @pro: product
# @r_match: reactant molecular species
# @p_match: product molecular species
get_reaction_same_struct <- function(reactions, res, keys, key, re, pro, r_match, p_match){
  r_fas <- get_acyl_chain(r_match)
  p_fas <- get_acyl_chain(p_match)

  shared_fa <- which(r_fas %in% p_fas)

  if(length(shared_fa) > 0){
    mat <- create_matrix()

    for(i in 1:length(shared_fa)){	
      fa <- r_fas[shared_fa[i]]
      rev <- paste0(re, "(", fa, ")")
      prov <- paste0(pro, "(", fa, ")")
      mat <- rbind(mat, c(rev, prov))
    }

    if(nrow(mat) > 0){
      if(is.null(keys) | !key %in% keys){
        res[[key]] <- as.data.frame(mat)
      }
    
      # if the reaction is reversable then add the reverted one as well
      rev_key <-  paste0(pro, ",", re)
      rev_react <- lp_reaction[which(lp_reaction$reaction == rev_key), ]
      if(nrow(rev_react) > 0){
        if(is.null(keys) | !rev_key %in% keys){
          mat <- mat[ , c("product_group","reactant_group")]
          if(!is.matrix(mat)){ # case: single row so reconstitute the matrix
            mat<-matrix(mat, nrow=1)
          }
          colnames(mat) <- header
          res[[rev_key]] <- as.data.frame(mat)
        }
        remove_item <- which(reactions == rev_key)
        if(length(remove_item)>0){
          reactions <- reactions[-remove_item]
        }
      }
    }
  }

  out <- list(res = res, reactions = reactions)
  return(out)
}


## For reactions involving FA or FACOA ##
# Check if a reaction at the molecular species level is valid
# @re: reactant (molecular species)
# @pro: product (molecular species)
# @fa_acyl: FA acyl chains
# @facoa_acyl: FACoA acyl chains
is_valid_reaction <- function(re, pro, fa_acyl, facoa_acyl){
  re_tokens <- strsplit(get_acyl_chain(re), split = ":")[[1]]
  pro_tokens <- strsplit(get_acyl_chain(pro), split = ":")[[1]]
  re_carbon <- as.numeric(re_tokens[1])
  re_bond <- as.numeric(re_tokens[2])
  pro_carbon <- as.numeric(pro_tokens[1])
  pro_bond <- as.numeric(pro_tokens[2])

  key <- paste0(get_species_class(re), ",", get_species_class(pro))
  item <- lp_reaction[which(lp_reaction$reaction==key), ]
  if(length(item) == 0){
    return(FALSE)
  }
  
  acyl_add <- as.logical(item$acyl_add)
  compound_require <- item$compound_require
  
  # if FACoA required in the reaction
  if(acyl_add){
    carbon_diff <- pro_carbon - re_carbon #calculate the number of carbons required by the FACoA for this reaction
    bond_diff <- pro_bond - re_bond
    acyl <- paste0(carbon_diff, ":", bond_diff)
    if(acyl %in% facoa_acyl){ #check if the FACoA exist in the dataset
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    # release FA
    if(!is.na(compound_require) & compound_require =='fa'){
      carbon_diff <- re_carbon - pro_carbon
      bond_diff <- re_bond - pro_bond
      acyl <- paste0(carbon_diff, ":",bond_diff)
      if(acyl %in% fa_acyl){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
  }
}


# Return groups of reactant and product that are transformable
# @r_match: reactant molecular species
# @p_match: product molecular species
# @fa_acyl: FA acyl chains
# @facoa_acyl: FACoA acyl chains
get_valid_reaction_groups <- function(r_match, p_match, fa_acyl, facoa_acyl){
  mat <- matrix(,,ncol=2)
  colnames(mat) <- header
  mat <- mat[-1,]

  while(length(r_match) > 0 & length(p_match) > 0){
    r <- r_match[1] # for each lipid molecular species in the reactant list
    p_v <- c()
    for(i in 1:length(p_match)){ # for each l. m. in the product list
      p <- p_match[i]
      if(is_valid_reaction(r, p, fa_acyl, facoa_acyl)){ #check if the reaction between these 2 lipids is possible
        p_v <- c(p_v, p)
      }
    }

    # if at least a reaction exist: add to the matrix and remove the associated product(s) species from the list of products
    if(length(p_v) > 0){
      mat <- rbind(mat,c(r,paste0(p_v,collapse = ",")))
      p_match <- p_match[-which(p_match %in% p_v)]
    }
    r_match <- r_match[-which(r_match == r)]
  }

  if(nrow(mat) == 0){
    return(NULL)
  }

  mat <- as.data.frame(mat)
  return(mat)
}


# Get the final list of reactions at the lipid molecular species, grouped at the subclass level
# @lipids: all lipid molecular species
get_group_reaction_list <- function(lipids){
  res <- list()
  reactions <- lp_reaction$reaction
  lp_group <- names(get_species_class(lipids))
  
  fa_facoa <- get_fa_coa_acyl_chain()
  fa_acyl <- fa_facoa$fa_acyl
  facoa_acyl <- fa_facoa$facoa_acyl
 
  # look at each reactions stored by BioPAN
  while(length(reactions) > 0){
    keys <- names(res)
    react <- reactions[1]
    dt <- lp_reaction[which(lp_reaction$reaction == react),]
    reactions <- reactions[-1]
    tokens <- strsplit(react, ",")[[1]]
    re <- tokens[1]
    pro <- tokens[2]
    
    # if the reactant and product subclasses of the reaction exist in the dataset
    if(re %in% lp_group & pro %in% lp_group){
      key <- react
      r_match <- as.vector(lipids[str_detect(lipids, regex(paste0("^", re,"\\((\\w+)|:\\)"), ignore_case = T))])
      p_match <- as.vector(lipids[str_detect(lipids, regex(paste0("^", pro,"\\((\\w+)|:\\)"), ignore_case = T))])

      # CL (cardiolipins) involved: BioPAN define as selected all reactant and product lipid molecular species
      if(grepl("CL", react)){
        res[[key]] <- list(reactant = paste0(r_match, collapse = ","), product = paste0(p_match, collapse = ","))
      
      # compute SPB/SPBP/dhSPB/dhSPBP differently because they do not have a structure
      }else if(grepl("SPB|SPBP|dhSPB|dhSPBP", react)){
        result <- get_reaction_sphingo(react, res, key, re, pro, r_match, p_match, fa_acyl, facoa_acyl)
        res <- result

      }else{
        # if the reaction does not include a FA/FACoA: search for reactions of the same structure
        if(is.na(dt$compound_require)){
          result <- get_reaction_same_struct(reactions, res, keys, key, re, pro, r_match, p_match)
          reactions <- result$reactions
          res <- result$res

        # if the reaction required a compound (FA/FACoA)
        }else{
          reaction_gr <- get_valid_reaction_groups(r_match, p_match, fa_acyl, facoa_acyl) 
          if(!is.null(reaction_gr)){
            if(is.null(keys) | !key %in% keys){
              res[[key]] <- reaction_gr
            }   
          }
        }
      }
    }
  }
  return(res)
}




#### Categories each molecular species (unrecognised/processed/unprocessed) ####

# Get molecular species which will be processed
# @gr: moelcular species
# @valid_lp_class: valid lipid classes or species for the FA
# @lp_react_list: list of reactions stored by BioPAN
# @return: vector
get_processed_groups <- function(gr, valid_lp_class, lp_react_list){
  processed_gr <- c()
  if(length(gr) == 0){
    return(processed_gr)
  }
  
  for(i in 1:length(gr)){
    if(gr[i] %in% valid_lp_class){
      # get all products
      prods <- lp_react_list[[gr[i]]]
      # if at least one product is found then the reactant will be processed
      inters <- intersect(prods, gr)
      if(length(inters) > 0){
        processed_gr <- c(processed_gr, gr[i])
        processed_gr <- union(processed_gr, inters)
      }
    }
  }
  return(processed_gr)
}


# Get a dataframe that maps molecular species names to standard format (old:FA 14:0 / new:FA(14:0))
# @species_names: species names
# @return: dataframe
get_species_standard_mapper <- function(species_names){
  new_col_name <- c()
  for(i in 1:length(species_names)){
    if(grepl(" ", species_names[i])){ # if there is a space between the lipid subclass and the structure
      splitter <- " "
    }else if(grepl("\\(", species_names[i])){ # there is brackets around the structure
      splitter <- "\\(|\\)"
    }else{
      new_col_name <- c(new_col_name, species_names[i])
      next
    }
    token_list <- strsplit(species_names[i], splitter)

    lp_name <- token_list[[1]][1]
    struct <- token_list[[1]][2]
    if(!grepl("SPB|SPBP|dhSPB|dhSPBP", lp_name)){ #lipid subclasses that don't have a structure
      new_name <- paste0(lp_name,"(", struct, ")")
      new_col_name <- c(new_col_name, new_name)
    }else{
      new_name <- paste0(lp_name)
      new_col_name <- c(new_col_name, new_name)
    }
    
  }
  mapper <- data.frame(old_name=species_names, new_name=new_col_name, check.names=FALSE)
  return(mapper)
}


# Get a list of species that are unrecognised
# @mapper: data frame
get_undefined_species <- function(species){ 
  undef <- c()
  database_cl <- database_classes
  database_cl <- sapply(database_cl, tolower)
  for(i in species){
    lipid <- strsplit(i,"\\(")
    specie <- tolower(lipid[[1]][1])
    if(!(specie %in% database_cl)){
      undef <- c(undef,i)
    }
  }
  return(undef)
}


# Get a mapping between BioPAN lipid molecular species nomenclature and the nomenclature of the dataset
# for molecular species output UNPROCESSED by LipidLynxX
# @species: unprocessed molecular species
# @species_dataset: molecular species correspondance with dataset nomenclature - to complete
get_biopan_dataset_nom_unprocessed <- function(species, species_dataset){
  if(length(species) > 0){
    for(i in 1:length(species)){
      species_dataset[[species[i]]] <- paste(species_dataset[[species[i]]], species[i], collapse = ",")
    }
  }
  return(species_dataset)
}

# for the other molecular species
# @species: unprocessed molecular species
# @species_dataset: molecular species correspondance with dataset nomenclature - to complete
# @dataset: molecular species in the dataset
# @biopan: molecular species converted into biopan nomenclature (or the same if unprocessed by LipidLynxX)
get_biopan_dataset_nom <- function(species, species_dataset, dataset, biopan){
  if(length(species) > 0){
    for(i in 1:length(species)){
      species_dataset[[species[i]]] <- paste(species_dataset[[species[i]]], dataset[which(biopan == species[i])], collapse = ",")
    }
  }
  return(species_dataset)
}
