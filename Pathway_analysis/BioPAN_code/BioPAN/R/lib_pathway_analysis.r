####################################################################
## Functions used in the processe_data and compute_pathways files ##
####################################################################



# Get lipid pathway reactant and product species
# @reactant: reactant (lipid class)
# @product: product (lipid class)
get_lipid_pathway_re_pro <- function(reactant, product){
  # get reaction file
  file_name <- tolower(paste0(reactant,",",product,".csv"))
  
  dt <- read.csv(paste0(out_file_path,"reaction/",file_name))

  if(!grepl("cl",file_name)){
    re <- c()
    pro <- c()
    for(i in 1:nrow(dt)){
      if(as.character(dt$reactant[i]) != "" | as.character(dt$product[i]) != ""){
        re <- c(re, as.character(dt$reactant[i]))
        if(grepl(",",dt$product[i])){
          pro <- c(pro, strsplit(as.character(dt$product[i]),",")[[1]])
        }else{
          pro <- c(pro,as.character(dt$product[i]))
        }
      }
    }
  }else{
    re <- str_split(as.character(dt$reactant), ",")[[1]]
    pro <- str_split(as.character(dt$product), ",")[[1]]
  }

  if(length(re[!is.na(re)]) == 0 | length(pro[!is.na(pro)]) == 0){
    return(list())
  }else{
    re_pro <- list(re=re, pro=pro)
    return(re_pro)
  }
}



# Filter lipid species by lipid class
# @lp_data: data
# @species: lipid species to be filtered
# @return: matrix
filter_lipid_class <- function(lp_data, species){
  if(is.vector(species)){
    if(is.vector(lp_data)){
      new_data <- lp_data[which(tolower(colnames(lp_data)) %in% tolower(species))]
    }else{
      new_data <- lp_data[,which(tolower(colnames(lp_data)) %in% tolower(species))]
    }
    
  }else{
    if(is.vector(lp_data)){
      new_data <- lp_data[which(tolower(colnames(lp_data)) == tolower(species))]
    }else{
      new_data <- lp_data[,which(tolower(colnames(lp_data)) == tolower(species))]
    }
  }
  
  return(new_data)
}



# Calculate t.test and catch errors
# @d_data_test: vector of calculated disease values (pro_disease/re_disease)
# @c_data_test: vector of calculated control values (pro_control/re_control)
# @alt: alternative hypothesis
# @is_paired: boolean, type of t-test
launch_t_test <- function(d_data_test, c_data_test, alt, is_paired){
  tryCatch({
      t.test(d_data_test, c_data_test, alternative=alt, paired=is_paired)
    },
    error = function(cond) cond
  )
}


# Extract z-score calculated in launch_t_test function
# @lp: count of 0 values in reactant
# @disease_sum: disease data (reactant and product)
# @control_sum: control data (reactant and product)
# @alt: alternative hypothesis
# @is_paired: boolean, type of t-test
extract_zscore <- function(l, ind, disease_sum, control_sum, alt, is_paired){
  z_score <- 0

  if(l > 0){
    d_data_test <- disease_sum[[1]][-ind]/disease_sum[[2]][-ind]  
    c_data_test <- control_sum[[1]][-ind]/control_sum[[2]][-ind]  
  }else{
    d_data_test <- disease_sum[[1]]/disease_sum[[2]]  
    c_data_test <- control_sum[[1]]/control_sum[[2]]
  }

  if(length(d_data_test) > 1 & length(c_data_test) > 1){
    t <- launch_t_test(d_data_test, c_data_test, alt, is_paired)

    if(grepl("Error", t)){
      if(!grepl("data are essentially constant", t)){
        return(NULL)
      } # else: values between reac an prod are the same. So no difference: assign 0 to z-score (line 90)
    }else{
      if(!is.nan(t$p.value) & !is.infinite(t$p.value)){
        z_score <- qnorm(1 - t$p.value)
        
      }
    }
  }

  return(z_score)
}



# Get lipid pathway zscore
# @re_pro: list of reactant and product species (lipid class)
# @control_dt: control samples data
# @disease_dt: disease samples data
# @alt: alternative hypothesis
# @is_paired: boolean, type of t-test
get_lipid_pathway_zscore <- function(re_pro, disease_dt, control_dt, alt = 'greater', is_paired = FALSE){
  re <- re_pro$re
  pro <- re_pro$pro

  r_disease_dt <- filter_lipid_class(disease_dt,re)
  p_disease_dt <- filter_lipid_class(disease_dt,pro)
  r_control_dt <- filter_lipid_class(control_dt,re)
  p_control_dt <- filter_lipid_class(control_dt,pro)

  disease_sum <- list()
  control_sum <- list()
 
  if(is.vector(p_disease_dt)){
    disease_sum[[1]] <- p_disease_dt
  }else{
    disease_sum[[1]] <- rowSums(p_disease_dt)  
  }
  if(is.vector(p_control_dt)){
    control_sum[[1]] <- p_control_dt
  }else{
    control_sum[[1]] <- rowSums(p_control_dt)
  }
  
  if(is.vector(r_disease_dt)){
    disease_sum[[2]] <- r_disease_dt
  }else{
    disease_sum[[2]] <- rowSums(r_disease_dt)  
  }
  if(is.vector(r_control_dt)){
    control_sum[[2]] <- r_control_dt
  }else{
    control_sum[[2]] <- rowSums(r_control_dt)
  }
  
  z_score <- 0
  ind1 <- which(disease_sum[[2]] == 0)
  ind2 <- which(control_sum[[2]] == 0)
  ind <- union(ind1,ind2)
  l <- length(ind)
  na_values <- sum(is.na(c(disease_sum[[1]], disease_sum[[2]], control_sum[[1]], control_sum[[2]])))

  if(l < length(disease_sum[[2]])-1 && na_values == 0){
    z_score <- extract_zscore(l, ind, disease_sum, control_sum, alt, is_paired)
    z_score <- format(round(z_score, 3), nsmall=3)
    if(is.null(z_score)){
      return(NULL)
    }
  }
  
  return(z_score)
}



# Get zscore from data to get "species pathways zscore" and "FA pathway zscore"
# @p_disease_dt: product disease samples data
# @p_control_dt: product control samples data
# @r_disease_dt: reactant disease samples data
# @r_control_dt: reactant control samples data
# @alt: alternative hypothesis
# @is_paired: boolean, type of t-test
get_zscore <- function(p_disease_dt, p_control_dt, r_disease_dt, r_control_dt, alt, is_paired){
  disease_sum <- list()
  control_sum <- list()
  disease_sum[[1]] <- p_disease_dt
  control_sum[[1]] <- p_control_dt
  disease_sum[[2]] <- r_disease_dt
  control_sum[[2]] <- r_control_dt
  z_score <- 0
  ind1 <- which(disease_sum[[2]] == 0)
  ind2 <- which(control_sum[[2]] == 0)
  ind <- union(ind1,ind2)
  l <- length(ind)

  if(l < length(disease_sum[[2]])-1){
    z_score <- extract_zscore(l, ind, disease_sum, control_sum, alt, is_paired)
    z_score <- format(round(z_score, 3), nsmall=3) # keep 3 values after decimal point
    if(is.null(z_score)){
      return(NULL)
    }
  }

  return(z_score)
}



# Get lipid species pathway zscore
# @re: reactant (lipid specie)
# @pro: product (lipid specie)
# @control_dt: control samples data
# @disease_dt: disease samples data
# @alt: alternative hypothesis
# @is_paired: boolean, type of t-test
get_lipid_species_pathway_zscore <- function(re, pro, disease_dt, control_dt, alt = 'greater', is_paired = FALSE){
  r_disease_dt <- filter_lipid_class(disease_dt,re)
  p_disease_dt <- filter_lipid_class(disease_dt,pro)
  r_control_dt <- filter_lipid_class(control_dt,re)
  p_control_dt <- filter_lipid_class(control_dt,pro)

  z_score <- get_zscore(p_disease_dt, p_control_dt, r_disease_dt, r_control_dt, alt, is_paired)
  z_score
}



# Get FA pathway zscore
# @reactant: reactant (FA specie)
# @product: product (FA specie)
# @control_dt: control samples data
# @disease_dt: disease samples data
# @is_paired: boolean, type of t-test
get_fa_pathway_zscore <- function(reactant, product, disease_dt, control_dt, alt, is_paired = TRUE){
  if(is.vector(disease_dt)){
    r_disease_dt <- disease_dt[which(tolower(colnames(disease_dt)) == tolower(reactant))]
    p_disease_dt <- disease_dt[which(tolower(colnames(disease_dt)) == tolower(product))] 
  }else{
     r_disease_dt <- disease_dt[,which(tolower(colnames(disease_dt)) == tolower(reactant))]
     p_disease_dt <- disease_dt[,which(tolower(colnames(disease_dt)) == tolower(product))] 
  }
  
  if(is.vector(control_dt)){
    r_control_dt <- control_dt[which(tolower(colnames(control_dt)) == tolower(reactant))]
    p_control_dt <- control_dt[which(tolower(colnames(control_dt)) == tolower(product))] 
  }else{
     r_control_dt <- control_dt[,which(tolower(colnames(control_dt)) == tolower(reactant))]
     p_control_dt <- control_dt[,which(tolower(colnames(control_dt)) == tolower(product))] 
  }
  
  z_score <- get_zscore(p_disease_dt, p_control_dt, r_disease_dt, r_control_dt, alt, is_paired)
  z_score
}



# Get cytoscape formatted ids
# @id: string
get_cyto_format_id <- function(id){
  res <- gsub(x = id, pattern = paste(c("\\(","\\)",":"), collapse = "|"), replacement = "")
  res
}


# Get color of the reaction depending on the score (positive or negative)
get_color <- function(color, score){
  if(score > 0){
    color <- c(color, '#24a19c')
  }else{
    color <- c(color, '#A634C7')
  }
  color
}


# Get reactions from pathways (= known in the litterature)
# @react_keys: lipid rections stored in the BioPAN reaction
get_reaction_from_pathway <- function(react_keys){
  pw_arr <- unique(as.character(pathway$pathway))
  res <- c()
  for(i in 1:length(pw_arr)){
    react <- strsplit(pw_arr[i], "->")[[1]]
    res <- c(res, react)
  }
  res <- unique(res)
  res[which(tolower(res) %in% tolower(react_keys))]
}



# Compute scored lipid class pathways
# @lp_data: lipid data
# @disease: disease samples name
# @control: control samples name
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @react_keys: lipid rections stored in the BioPAN reaction
# @return: data frame
get_scored_lipid_pathways <- function(lp_data, disease , control, alt = 'greater', is_paired = FALSE, type = 'reaction', react_keys){
  # get set of reactants and products
  reactant <- vector()
  product <- vector()
  weight <- vector()
  color <- vector()
  
  if(type == 'pathway'){
  react <- get_reaction_from_pathway(react_keys)
  }else{
    react <- react_keys
  }

  if(length(react) > 0){
    for(j in 1:length(react)){
      token <- strsplit(react[j], ",")[[1]]
      re <- token[1]
      pro <- token[2]
      
      re_pro <- get_lipid_pathway_re_pro(re, pro)

      if(length(re_pro) == 2){
        score <- get_lipid_pathway_zscore(re_pro, lp_data[[disease]], lp_data[[control]], alt, is_paired)
        if(!is.null(score)){
          reactant = c(reactant, re)
          product = c(product, pro)
          weight <- c(weight, score)
          
          color <- get_color(color, score)
          class <- lp_reaction[lp_reaction$reaction == react[j],]$class
        }
      }
    }

    if(length(reactant) > 0){
      id <- paste0(get_cyto_format_id(reactant), get_cyto_format_id(product))
      dt <- data.frame(id = id, reactant = reactant, product = product, score = weight, color = color)
    }else{
      dt <- data.frame()
    }
    
  }else{
    dt <- data.frame()
  }

  return(dt)
}



# Compute scored lipid species pathways
# @lp_data: lipid data
# @disease: disease samples name
# @control: control samples name
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @return: data frame
convert_scored_lipid_species_pathway <- function(react, disease, control, alt = 'greater', is_paired = FALSE){
  # get set of reactants and products
  reactant <- vector()
  product <- vector()
  weight <- vector()
  color <- vector()
  
  for(j in 1:length(react)){ 
    file_name <- paste0(out_file_path,"reaction/",tolower(react[j]),".csv")
    
    if(file.exists(file_name)){
      r <- read.csv(file_name)
      re_arr <- as.character(r$reactant)

      for(k in 1:length(re_arr)){
        pro_arr <- strsplit(as.character(r$product[k]),",")[[1]]
        for(i in 1:length(pro_arr)){
          score <- get_lipid_species_pathway_zscore(re_arr[k], pro_arr[i], lp_data[[disease]], lp_data[[control]], alt, is_paired)
          if(!is.null(score)){
            reactant = c(reactant, re_arr[k])
            product = c(product, pro_arr[i])
            weight <- c(weight, score)
            color <- get_color(color, score)
          }
        }
      }
    }
  }

  dt <- NULL
  if(length(weight) > 0){   
    id <- paste0(get_cyto_format_id(reactant), get_cyto_format_id(product))
    dt <- data.frame(id = id, reactant = reactant, product = product, score = weight, color = color)
  }

  return(dt)
}



# Build the pathway and write it in another way to match with the database. Remove the structure for species
# @pathway: string pathway
# @node: string, lipid to be added to the pathway
get_pw_struct <- function(pathway, node){
  path <- c()
  if(length(pathway) > 1){
    for(j in 1:(length(pathway)-1)){
      if(grepl(":",pathway[j])){
        pathway[j] <- strsplit(pathway[j],"\\(")[[1]]
      }
      if(grepl(":",pathway[j+1])){
        pathway[j+1] <- strsplit(pathway[j+1],"\\(")[[1]]
      }
      path <- paste0(path, pathway[j],",",pathway[j+1],"->")
    }
  }
  if(grepl(":",node)){
    node <- strsplit(node,"\\(")[[1]][1]
  }
  if(grepl(":",pathway[length(pathway)])){
    pathway[length(pathway)] <- strsplit(pathway[length(pathway)],"\\(")[[1]]
  }
  
  path <- paste0(path,pathway[length(pathway)],",",node)
  return(path)
}


# Compute significant lipid pathway
# @dataL: list to return
# @pathway: pathway to compute
# @edgeW: edges
# @root: root node
# @disease_dt: control data
# @control_dt: disease data
# @p_value: p_value
# @alt: alternative hypothesis (greater/less)
# @is_paired: type of t-test
# @typeL: lipid/species/fa
# @subset: reaction or pathway (for lipid/species "typeL") or null (for fa)
# @json_obj: json containing edges/nodes info
# @return: list
get_sig_pathways <- function(dataL, pathwayC, edgeW, root, disease_dt, control_dt, p_value = 0.05, alt = 'greater', is_paired = TRUE, typeL = 'lipid', subset, json_obj){
  pw <- as.vector(t(pathway$pathway))
  path <- c()
  
  visitted <- dataL[[2]]
  visitted[root] <- TRUE
  dataL[[2]] <- visitted
  edgeWL <- edgeW[[root]]
  # get list of its neighbors with decreasing order of ZScores
  neighborL <- sort(edgeWL, decreasing = TRUE)
  # get the next node
  nodeL <- names(neighborL)

  if(length(nodeL)>0){
    for(i in 1:length(nodeL)){
      if(!visitted[nodeL[i]]){
        if(subset == "pathway"){
          path <- get_pw_struct(pathwayC, nodeL[i])
        } 

        if((subset != "pathway") || (subset == "pathway" && length(grep(path,pw)) > 0)){
          zscore <- get_sub_pathway_zscore(c(pathwayC,nodeL[i]), disease_dt, control_dt, alt, is_paired, typeL, json_obj)

          # pathway is considered as signifiant if zscore is greater than a critical value
          if(as.numeric(zscore) > qnorm(1-p_value)){
            pathVec <- dataL[[1]]
            path <- paste(pathwayC, collapse = right_arrow)
            apathway <- c(pathwayC,nodeL[i])

            newPath <- paste(apathway, collapse = right_arrow)
            if(length(pathVec)>0){ 
              pathN <- names(pathVec)
              # if pathway has been computed then ignore it
              if(path %in% pathN){
                pathVec <- pathVec[-which(pathN==path)]
              }
            }
            pathVec[newPath] <- zscore
            dataL[[1]] <- pathVec

            dataL <- get_sig_pathways(dataL, apathway, edgeW, nodeL[i], disease_dt, control_dt,p_value, alt, is_paired, typeL, subset, json_obj)
          }
        }
      }
    }
  }
  return(dataL)
}



# Calculate Z-score of subpathway
# @pathway: string
# @disease_dt: control data
# @control_dt: disease data
# @p_value: p_value
# @alt: alternative hypothesis (greater/less)
# @is_paired: type of t-test
# @type: lipid/species/fa
# @json_obj: json containing edges/nodes info
get_sub_pathway_zscore <- function(pathway, disease_dt, control_dt, alt = 'greater', is_paired = FALSE, type = 'lipid', json_obj){
  size <- length(pathway)-1
  z <- vector()
  z_score <- 0
  if(type == 'species'){
    react <- json_obj$react_df
  }else{
    if(class(json_obj) == 'json'){
      react <- fromJSON(json_obj)$edges[[1]]
    }else{
      react <- json_obj$edges
      names(react) <- c("id","source","target","weight","color")
    } 
  }

  for(i in 1:size){
    # get z-score for each edge
    if(type == 'lipid'){
      z[i] <- react$weight[which(react$id == paste0(pathway[i],pathway[i+1]))]
    }else if(type == 'species'){
      z[i] <- react$score[which(react$id == gsub("[()]|:", "", paste0(pathway[i],pathway[i+1])))]
    }else{
      z[i] <- react$weight[which(react$id == gsub("[()]|:", "", paste0(pathway[i],pathway[i+1])))]
    }
    z_score <- sum(as.numeric(z))/sqrt(size)
    z_score <- format(round(z_score, 3), nsmall=3)
  }
  
  return(z_score)
}



# Compute all significant fatty acid pathway
# @disease: disease samples data
# @control: control samples data
# @p_value: p-value
# @alt: alternative hypothesis (greater/less)
# @is_paired: type of t-test
# @json_obj: json containing edges/nodes info
# @return: data frame
find_all_fatty_acid_sig_pathways <- function(disease_dt, control_dt , p_value = 0.05, alt = 'greater', is_paired = FALSE, json_obj){
  reactant <- vector()
  product <- vector()
  weight <- vector()
  color <- vector()

  fa_species <- tolower((colnames(disease_dt)))
  fa_reactions <- fa_reaction$reaction

  for(j in 1:length(fa_reactions)){
    token <- strsplit(fa_reactions[j],",")[[1]]
    re <- token[1]
    pro <- token[2]
    if(tolower(re) %in% fa_species & tolower(pro) %in% fa_species){
      score <- get_fa_pathway_zscore(re, pro, disease_dt, control_dt, alt, is_paired)
      if(!is.null(score)){
        reactant = c(reactant, re)
        product = c(product, pro)
        weight <- c(weight, as.numeric(score))
      }
      
    }else{
      next
    }  
  }
  
  # generate a list of found subpathways with its z_score
  vPathway <- c()
  vZscore <- c()
  # build a graph
  df <- cbind(reactant, product)
  gr <- ftM2graphNEL(df, W=weight, V=NULL, edgemode="directed")
  node <- nodes(gr)
  edgeW <- edgeWeights(gr)

  # if only a row in df: a dot is added in the reactant name in edgeW
  if(nrow(df) == 1){
    names(edgeW[[df[[1,1]]]]) <- df[[1,2]]
  }

  visitted <- vector()
  
  for(i in 1:length(node)){
    visitted[node[i]] <- FALSE  
  }
  dataL <- list()
  dataL[[1]] <- vector()
  for(i in 1:length(node)){
    visitted[which(visitted==TRUE)] <- FALSE
    dataL[[2]] <- visitted
    if(!visitted[node[i]]){
      pathway <- c(node[i])
      dataL <- get_sig_pathways(dataL,pathway, edgeW, node[i], disease_dt, control_dt, p_value, alt,is_paired, type = 'fa', "", json_obj)
    }
  }
  vZscore <- c(vZscore,dataL[[1]])
  vPathway <- c(vPathway,names(dataL[[1]]))
  
  dt <- data.frame(pathway = factor(vPathway), score = vZscore)
  dt <- dt[order(dt$score, decreasing = FALSE),]
  return(dt)
}



# Compute all significant lipid pathway either active or suppress
# @disease_dt: disease samples data
# @control_dt: control samples data
# @p_value: p-value
# @alt: alternative hypothesis (greater/less)
# @is_paired: type of t-test
# @type: subset of lipid data (reaction/pathway)
# @json_obj: json containing edges/nodes info
# @return: data frame
find_all_lipid_sig_pathways <- function(disease_dt, control_dt , p_value = 0.05, alt = 'greater', is_paired = FALSE, type = 'reaction', json_obj = NULL){
  reactant <- vector()
  product <- vector()
  weight <- vector()
  color <- vector()
  
  lp_names <- database_classes_upper
  obj <- fromJSON(json_obj, flatten=TRUE)
  edges <- obj$edges

  for(j in 1:length(edges$data.source)){
    reactant = c(reactant, as.character(lp_names$classes[which(edges$data.source[j] == lp_names$lower)]))
    product = c(product, as.character(lp_names$classes[which(edges$data.target[j] == lp_names$lower)]))
    weight <- c(weight, as.numeric(edges$data.weight[j]))
  }
  
  # generate a list of found subpathways with its z_score
  vPathway <- c()
  vZscore <- c()
  # build a graph
  df <- cbind(reactant,product)

  classes <- ftM2graphNEL(df, W=weight, V=NULL, edgemode="directed")
  node <- nodes(classes)
  edgeW <- edgeWeights(classes)

  # if only a row in df: a dot is added in the reactant name in edgeW
  if(nrow(df) == 1){
    names(edgeW[[df[[1,1]]]]) <- df[[1,2]]
  }

  visitted <- vector()
  
  for(i in 1:length(node)){
    visitted[node[i]] <- FALSE  
  }

  dataL <- list()
  dataL[[1]] <- vector()
  for(i in 1:length(node)){
    visitted[which(visitted==TRUE)] <- FALSE
    dataL[[2]] <- visitted
    if(!visitted[node[i]]){
      pathway <- c(node[i])
      dataL <- get_sig_pathways(dataL,pathway, edgeW, node[i], disease_dt, control_dt, p_value, alt,is_paired, type = 'lipid', type, obj)
    }
  }

  vZscore <- c(vZscore,dataL[[1]])
  vPathway <- c(vPathway,names(dataL[[1]]))
   
  dt <- data.frame(pathway = factor(vPathway), score = vZscore)
  dt <- dt[order(dt$score, decreasing = TRUE),]
  return(dt)
}



# Compute all significant lipid species pathway either active or suppress
# @disease_dt: disease samples data
# @control_dt: control samples data
# @p_value: p-value
# @alt: alternative hypothesis (greater/less)
# @is_paired: type of t-test
# @type: subset of lipid data (reaction/pathway)
# @list: list containing edges/nodes info
# @return: data frame
find_all_lipid_species_sig_pathways <- function(disease_dt, control_dt , p_value = 0.05, alt = 'greater', is_paired = FALSE, type = 'reaction', list){
  reactant <- vector()
  product <- vector()
  weight <- vector()
  color <- vector()
  id <- vector()
  label <- vector()
  
  react <- list$react_df
  for(i in 1:nrow(react)){
    reactant <- c(reactant, as.character(react$reactant[i]))
    product <- c(product, as.character(react$product[i]))
    weight <- c(weight, as.numeric(react$score[i]))
  }

  if(length(reactant) > 0 & length(product) > 0){
    # generate a list of found subpathways with its z_score
    vPathway <- c()
    vZscore <- c()
    # build a graph
    df <- cbind(reactant,product)
    gr <- ftM2graphNEL(df, W=weight, V=NULL, edgemode="directed")
    node <- nodes(gr)
    edgeW <- edgeWeights(gr)

    # if only a row in df: a dot is added in the reactant name in edgeW
    if(nrow(df) == 1){
      names(edgeW[[df[[1,1]]]]) <- df[[1,2]]
    }

    visitted <- vector()

    for(i in 1:length(node)){
      visitted[node[i]] <- FALSE  
    }

    dataL <- list()
    dataL[[1]] <- vector()
    for(i in 1:length(node)){
      visitted[which(visitted==TRUE)] <- FALSE
      dataL[[2]] <- visitted
      if(!visitted[node[i]]){
        pathway <- c(node[i])
        dataL <- get_sig_pathways(dataL,pathway, edgeW, node[i], disease_dt, control_dt, p_value, alt,is_paired, type = 'species', type, list)
      }
    }
    vZscore <- c(vZscore,dataL[[1]])
    vPathway <- c(vPathway,names(dataL[[1]]))

    dt <- data.frame(pathway = factor(vPathway), score = vZscore)
    dt <- dt[order(dt$score, decreasing = TRUE),]
    return(dt)
  }
}



# Compute the most significant pathways
# @disease: disease samples
# @control: control samples
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @type: lipid/species/fa
# @return: data frame
find_most_sig_pathways <- function(result, disease_dt, control_dt, alt, is_paired = FALSE, type)
{
  paths <- strsplit(as.character(result$pathway),right_arrow)

  selected <- c()
  vPathway <- c()
  vZscore <- c()
  for(i in 1:length(paths)){
    if(!paths[[i]][1] %in% selected){
      selected <- c(selected, paths[[i]][1])
      vPathway <- c(vPathway, as.character(result$pathway[i]))
      vZscore <- c(vZscore, result$score[i])
    }else{
      pos <- match(paths[[i]][1],selected)
      compare <- strsplit(vPathway[pos],right_arrow)[[1]]
      samePos <- NULL

      for(j in 1:length(paths[[i]])){
        if (paths[[i]][j] != compare[j]){
          samePos <- j - 1
          break
        }
      }

      if(samePos){
        if(type == 'lipid'){
          re_pro_Path <- get_lipid_pathway_re_pro(compare[samePos],compare[samePos+1])
          re_pro_New <- get_lipid_pathway_re_pro(paths[[i]][samePos],paths[[i]][samePos+1])
          if(length(re_pro_Path) == 2 & length(re_pro_New) == 2){
            z_Path <- get_lipid_pathway_zscore(re_pro_Path, disease_dt, control_dt, alt, is_paired)
            z_New <- get_lipid_pathway_zscore(re_pro_New,disease_dt, control_dt, alt, is_paired)
          }
        }else if(type == 'species'){
          z_Path <- get_lipid_species_pathway_zscore(compare[samePos],compare[samePos+1],disease_dt, control_dt, alt, is_paired)
          z_New <- get_lipid_species_pathway_zscore(paths[[i]][samePos],paths[[i]][samePos+1],disease_dt, control_dt, alt, is_paired)
        }else{ #FA
          z_Path <- get_fa_pathway_zscore(compare[samePos],compare[samePos+1],disease_dt, control_dt, alt, is_paired)
          z_New <- get_fa_pathway_zscore(paths[[i]][samePos],paths[[i]][samePos+1],disease_dt, control_dt, alt, is_paired)
        }
        
        if(!is.null(z_Path) && !is.null(z_New)){
          if(z_New > z_Path){
            vPathway[pos] <- as.character(result$pathway[i])
            vZscore[pos] <- result$score[i]
          }
        }
      }
    }
  }

  dt <- data.frame(pathway = factor(vPathway), score = vZscore)
  dt <- dt[order(dt$score, decreasing = TRUE),]
  return(dt)
}



# Compute scored fatty acids pathways
# @fa_data: fatty acid data
# @disease: disease samples
# @control: control samples
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @return: data frame
get_scored_fa_pathways <- function(fa_data, disease, control , alt, is_paired = FALSE){
  # get set of reactants and products
  reactant <- vector()
  product <- vector()
  weight <- vector()
  color <- vector()
  if(is.vector(fa_data[[1]])){
    fa <- tolower(names(fa_data[[1]]))
  }else{
    fa <- tolower(colnames(fa_data[[1]]))  
  }

  fa_reactions <- fa_reaction$reaction
  for(j in 1:length(fa_reactions)){
    token <- strsplit(fa_reactions[j],",")[[1]]
    re <- token[1]
    pro <- token[2]

    if(tolower(re) %in% fa & tolower(pro) %in% fa){
      score <- get_fa_pathway_zscore(re, pro, fa_data[[disease]], fa_data[[control]],alt, is_paired)
      if(!is.null(score)){
        reactant = c(reactant, re)
        product = c(product, pro)
        weight <- c(weight, score)
        color <- get_color(color, score)
      }
    }else{
      next
    }
  }

  id <- paste0(get_cyto_format_id(reactant), get_cyto_format_id(product))
  dt <- data.frame(id = id, reactant = reactant, product = product, score = weight, color = color)
  return(dt)
}



# Return shape of nodes in function of the lipid classification
# @class: classification of the lipid
get_lp_shape <- function(class){
  if(class == "Glycerolipids and Glycerophospholipids"){
    shape <- "ellipse"
  }else if(class == "Sphingolipids"){
    shape <- "rectangle"
  }else if(class == "FA"){
    shape <- "triangle"
  }else{
    shape <- "octagon"
  }
  return(shape)
}


# Get a node data frame from lipid subclass data
# @lp_names: lipid subclass names
# @return: data frame
get_lipid_node_data_frame <-function(processed_class){
  nodes <- as.data.frame(read.csv("resources/lipid_nodes.csv", header = T))
  nodes <- nodes[which(tolower(nodes$node) %in% tolower(processed_class)), ]

  # get systematic names if there are any
  names <- as.vector(nodes$systematic_name)
  lm_id <- as.vector(nodes$LM_id)
  lmid_vec <- c()
  label_vec <- c()
  name_vec <- c()
  shape <- c()
  lp_class <- lp_classification
  
  for(i in 1: nrow(nodes)){
    dt <- nodes[i,]
    node_name <- as.character(dt$node)
    label_vec <- c(label_vec, node_name)
    lmid_vec <- c(lmid_vec, lm_id[i])
    name_vec <- c(name_vec, names[i])
    class <- lp_class[lp_class$lp == node_name, "class"]
    shape <- c(shape, get_lp_shape(class))
    
  }
  dt <- data.frame(id = tolower(label_vec), lm_id = lm_id, label = label_vec, name = name_vec, shape = shape)
  return(dt)
}



# Convert scored lipid pathways from lipid data to json
# @lp_data: lipid data
# @node_df: nodes dataframe
# @control: control samples
# @disease: disease samples
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @type: subset of lipid data (reaction/pathway)
# @react_keys: lipid rections stored in the BioPAN reaction
# @return: json object
convert_scored_lipid_pathway_to_json <- function(lp_data, node_df, disease, control, alt='greater', is_paired = FALSE, type = 'reaction', react_keys){
  react_df <- get_scored_lipid_pathways(lp_data, disease, control, alt, is_paired, type, react_keys)
  
  if(length(react_df) > 0){
    json_obj <- convert_pathway_to_json(react_df, node_df)
    json_obj
  }else{
    NULL
  }
}


# Get reactions from pathways (= known in the litterature). Filtered to remove reactions including CL lipids because we don't know which specie gives which species
# @react_keys: lipid rections stored in the BioPAN reaction
get_reaction_from_pathway_filter <- function(react_keys){
  pw_arr <- unique(as.character(pathway$pathway))
  res <- c()
  for(i in 1:length(pw_arr)){
    react <- strsplit(pw_arr[i], "->")[[1]]
    if(!grepl("CL",pw_arr[i])){
      res <- c(res, react)
    }
  }
  res <- unique(res)
  res <- res[which(tolower(res) %in% tolower(react_keys))]
  return(res)
}


# Convert scored lipid species pathways from lipid data to json
# @lp_data: lipid data
# @control: control samples
# @disease: disease samples
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @type: subset of lipid data (reaction/pathway)
# @react_keys: lipid rections stored in the BioPAN reaction
# @return: list
convert_scored_lipid_species_pathway_to_json <- function(node_df, disease, control, alt = 'greater', is_paired = FALSE, type = 'reaction', react_keys){
  if(type == 'reaction'){
    react_df <- lp_reaction   
    react <- react_df$reaction[-which(grepl("CL",react_df$reaction))] #CL species don't appear on the lipid species graph as we don't know which species produces which species
  } else{
    react <- get_reaction_from_pathway_filter(react_keys)
  }
  react_df <- convert_scored_lipid_species_pathway(react, disease, control, alt, is_paired)

  list <- list("react_df" = react_df, "node_df" = node_df)
  return(list)
}


# Convert scored fatty acid pathways from fa data to json
# @fa_data: fatty acid data
# @node_df: node data frame
# @disease: disease samples
# @control: control samples
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @return: json object
convert_scored_fa_pathway_to_json <- function(fa_data, node_df, disease, control, alt = 'greater', is_paired = FALSE){
  react_df <- get_scored_fa_pathways(fa_data, disease, control, alt, is_paired)
  node_df$shape = rep("triangle",nrow(node_df))
  json_obj <- convert_pathway_to_json(react_df, node_df)
  return(json_obj)
}


# Convert scored fatty acid pathways from fa data to json
# @fa_data: fatty acid data
# @node_df: node data frame
# @disease: disease samples
# @control: control samples
# @alt: alternative hypothesis
# @is_paired: type of t-test
# @return: json object
convert_scored_lipid_pathway_to_json <- function(lp_data, node_df, disease, control, alt='greater', is_paired = FALSE, type = 'reaction', react_keys){
  react_df <- get_scored_lipid_pathways(lp_data, disease, control, alt, is_paired, type, react_keys)
  
  if(length(react_df) > 0){
    json_obj <- convert_pathway_to_json(react_df, node_df)
    return(json_obj)
  }else{
    return(NULL)
  }
}


# Add back bone to lipid species structures (Cer/dhCer/SM/dhSM)
# @s: lipid species
# @splitter: splitter character
# @back_bone: string to add
add_d_to_struct <- function(s, splitter = ":", back_bone = "d18:0"){
  acyls <- get_acyl_chain(s)
  gr <- get_species_class(s[1])
  res <- c()
  for(i in 1:length(acyls)){
    if(back_bone == "d18:0"){
      gr <- gsub("dh|Dh|DH|dH","",gr)
    }
    res <- c(res, paste0(gr,"(",back_bone, "/",acyls[i],")"))
  }
  return(res)
}


# Get infos on id, label and shape on each molecular species
# @species_name: species names
# @return: data frame (id | label | shape)
get_species_df <- function(species_names){
  id <- c()
  shape <- c()
  for(i in 1:length(species_names)){
    id <- c(id, get_cyto_format_id(tolower(species_names[i])))
    class <- strsplit(species_names[i], " |\\(")[[1]][1]
    classification <- lp_classification$class[toupper(lp_classification$lp) == toupper(class)]
    shape <- c(shape, get_lp_shape(classification))
  }
  new_df <- data.frame(id = id, label = species_names, shape = shape)
  return(new_df)
}



# Convert a data frame to json object
# @data_fr: id | reactant | product | score | color
# @node_fr: id | label | shape (lm_id for the subclass data frames)
# @return: json object
convert_pathway_to_json <- function(data_fr, node_fr){
  nodes <- list()
  edges <- list()
  data <- list()

  j = 1 #nodes counter
  k = 1 #edges counter

  valid_node <- c()

  data_fr_product <- as.character(data_fr$product)
  data_fr_reactant <- as.character(data_fr$reactant)
  node_fr_label <- as.character(node_fr$label)

  for(i in 1:length(data_fr_reactant)){
    reac_prod = c(data_fr_reactant[i], data_fr_product[i])
    
    for(value in reac_prod){
      if(!value %in% valid_node){ #if the reactant/product is not in the list of nodes
        if(tolower(value) %in% tolower(node_fr_label)){
          dt <- node_fr[which(tolower(node_fr_label) == tolower(value)),]
          node <- list(data = list(id = dt$id, lm_id = dt$lm_id, label = dt$label, name = dt$name, shape = dt$shape))
        }else{
          id <- gsub("[()]|:","",value)
          node <- list(data = list(id = tolower(id), lm_id = "", label = value, name = "", shape = dt$shape))
        }
        # add node to list
        nodes[[j]] <- node
        j <- j + 1
        valid_node <- c(valid_node, value)
      }
    }

    #Add edges
    edge <- list(data = list(id = data_fr$id[i], source = gsub("[()]|:","",tolower(data_fr_reactant[i])), target = gsub("[()]|:","",tolower(data_fr_product[i])), weight = data_fr$score[i], color = data_fr$color[i]))
    edges[[k]] <- edge
    k <- k + 1
  }

  gr <- list(nodes = nodes, edges = edges)
  json_obj <- toJSON(gr, pretty=T, auto_unbox=T)

  return(json_obj)
}


# Get highlighted significant pathways
# @dt: significant pathway data frame
# @return: list
get_highlighted_sig_pathways <- function(dt){
  nodeL <- c()
  edgeL <- c()
  for(i in 1: nrow(dt)){
    d <- dt[i,]
    pathway <- as.character(d$pathway)
    score <- d$score
    nodes <- strsplit(pathway, right_arrow)[[1]]
    nodeL <- union(nodeL, nodes)
    l <- length(nodes)
    for(j in 1:(l-1)){
      re <- nodes[j]
      pro <- nodes[j+1]
      edgeL = union(edgeL, paste0(get_cyto_format_id(re),get_cyto_format_id(pro)))
    }
  }
  format_ids <- tolower(nodeL)
  format_ids <- get_cyto_format_id(format_ids)
  cy_nodes <- paste(paste0("#", format_ids), collapse=",")
  cy_edges <- paste(paste0("#", edgeL), collapse=",")

  gr = list(nodes = cy_nodes, edges = cy_edges)
  return(gr)
}


# Return lipid subclass of a molecular species
# @species: string containg the lipid molecular species
get_class_species <- function(species){
  class <- strsplit(species, "\\(")[[1]][1]
  return(class)
}


# Get gene information
# @pathway_df: data frame
# @type: lipid or fatty acid
# @subset: subset of lipid data (reaction, pathway or fa)
# @res: list
get_genes_info <- function(pathway_df, type = 'lipid', subset){
  gene_col <- c()
  if(type == 'lipid'){
      react <- lp_reaction
  }else{
      react <- fa_reaction
  }
  
  for(i in 1: nrow(pathway_df)){
    dt <- pathway_df[i,]
    pw <- as.character(dt$pathway)
    pw_tk <- strsplit(pw, right_arrow)[[1]]
    len <- length(pw_tk)
    geneL <- c()
    for(j in 1:(len-1)){
      re <- pw_tk[j]
      pro <- pw_tk[j+1]

      if(type == 'lipid'){
        if(grepl("\\(|\\)",re)){
          re <- get_class_species(re)
        }
        if(grepl("\\(|\\)",pro)){
          pro <- get_class_species(pro)
        } 
      }
      key <- paste0(re, ",", pro)
      res <- react[which(tolower(react$reaction) == tolower(key)),]
  
      if(nrow(res)>0){
        gene <- res$genes
        geneL <- c(geneL, gene)
      }  
    }
    gene_col <- c(gene_col, paste(geneL, collapse=","))
  }
  
  pathway_df$gene = gene_col
  pathway_df = pathway_df[order(pathway_df$score, decreasing = TRUE),]
  pathways = list()
  for(i in 1:nrow(pathway_df)){
    dt <- pathway_df[i,]
    if(subset == "pathway"){
      pw <- list(data = list(pathway = dt$pathway, score = dt$score, gene = dt$gene, class = dt$class))
    }else{
      pw <- list(data = list(pathway = dt$pathway, score = dt$score, gene = dt$gene))
    }
    
    # add edge to list
    pathways[[i]] <- pw
  }
  res <- list(pathways = pathways)
  return(res)
}


# Get classification information (for subset of lipid data equal to pathway)
# @pathway_df: data frame 
# @res: data frame
get_class_info <- function(pathway_df){
  path <- as.character(pathway$pathway)
  cl <- as.character(pathway$class)
  sub_cl <- as.character(pathway$sub_class)
    
  pathway_df$pathway <- as.character(pathway_df$pathway)
  cl_subcl <- c()
  
  for(i in 1:length(pathway_df$pathway)){
    out <- c()
    info <- c()
    nb <- 0
    bool <- FALSE
    reactions <- strsplit(pathway_df$pathway[i], right_arrow)
    
    for(j in 1:(length(reactions[[1]])-1)){
      if(grepl(":",reactions[[1]][j])){
        reactions[[1]][j] <- strsplit(reactions[[1]][j], "\\(")[[1]][1]
      }
      if(grepl(":",reactions[[1]][j+1])){
        reactions[[1]][j+1] <- strsplit(reactions[[1]][j+1], "\\(")[[1]][1]
      }
      out <- paste0(out,reactions[[1]][j],",",reactions[[1]][j+1],"->")
    }
    out <- gsub(" ","",substr(out, 1,nchar(out)-2))    
    
    for(l in 1:length(path)){
      if(grepl(out,path[l])){
        if(nb > 0){
          if(grepl(sub_cl[l],info)){
            bool <- TRUE
          }else{
            bool <- FALSE
          }
        } 
        if(!bool){
          info <- c(paste0(info, ", ", sub_cl[l], " (", cl[l], ")"))
          nb <- nb + 1
        }      
      }
    }
    info <- substr(info, 3,nchar(info))
    cl_subcl <- c(cl_subcl, info)
  }

  pathway_df$class <- cl_subcl 
  return(pathway_df)
}



write_isPaired <- function(is_paired){
  if(!is_paired){
    paired <- "notpaired"
  }else{
    paired <- "paired"
  }
  return(paired)
}


# Create all the lipid files names
lp_files_names <- function(out_file_path, subset, disease, control, status, p_value, paired){
  type <- c("lp_class", "lp_species")
  status_detail <- c("", "_most")
  extension <- c("", "_tbl")
  files <- c()

  for(i in 1:length(type)){
    for(j in 1:length(status_detail)){
      for(k in 1:length(extension)){
        file <- paste0(out_file_path, paste0(c(type[i], subset, disease, control), collapse="_"), paste0(c(status_detail[j], status, p_value, paired), collapse="_"), extension[k], ".json")
        files <- c(files, file)
      }
    }
  }
  return(files)
}


# Create all the FA files names
fa_files_names <- function(out_file_path, disease, control, status, p_value, paired){
  status_detail <- c("", "_most")
  extension <- c("", "_tbl")
  files <- c()

  for(i in 1:length(status_detail)){
    for(j in 1:length(extension)){
      file <- paste0(out_file_path, paste0(c("fa", disease, control), collapse="_"), paste0(c(status_detail[i], status, p_value, paired), collapse="_"), extension[j], ".json")
      files <- c(files, file)
    }
  }
  return(files)
}



# Run function to get highlighted pathways and get gene data
run_highlighted_genes <- function(dt, sub_set, type, file_to_save1, file_to_save2){
  gr <- get_highlighted_sig_pathways(dt)
  json_obj <- toJSON(gr, pretty=T, auto_unbox=T) 
  save_json_to_file(json_obj, file_to_save1)
  
  if(sub_set == "pathway"){
    dt <- get_class_info(dt) # get classification of the pathway(s)
  }
  
  gene_info <- get_genes_info(dt, type, sub_set)
  json_obj <- toJSON(gene_info, pretty=T, auto_unbox=T)
  save_json_to_file(json_obj, file_to_save2)
  
  return(dt)
}



# Find reactions/pathways - calculation - create files
# @mode: process (1st calculation or change paired-data) or compute (change p-value)
create_files <- function(mode, lp_data, fa_data, node_lp_df, node_species_df, fa_node_df, groups, is_paired, p_value, react_keys, out_file_path){
 paired <- write_isPaired(is_paired)
 write_p_value <- as.character(format(p_value, scientific=F))

  alts <- c("greater", "less")
  status <- c("active", "suppressed")
  subset <- c("reaction","pathway")
  l <- length(groups)

  valid_lp_reaction <- FALSE
  valid_fa_reaction <- FALSE
  subset_reaction <- data.frame(reaction=FALSE, pathway=FALSE)

  for(k in 1:length(alts)){
    alt <- alts[k]
    # to compare each group with the others (i/j)
    for(i in seq(1,l)){
      for(j in seq(1,l)[-i]){
        control <- groups[i]
        disease <- groups[j]
        disease_dt <- lp_data[[disease]]
        control_dt <- lp_data[[control]]
        fa_disease_dt <- fa_data[[disease]]
        fa_control_dt <- fa_data[[control]]

        if(nrow(disease_dt) == nrow(control_dt)){ # if same number of samples in both conditions
          ############### For reactions and pathways subset of lipid data ###############
          for(sub_set in subset){
            
            if(!is.null(lp_data) && !is.null(node_lp_df)){
              valid_lp_reaction <- TRUE
              subset_reaction[subset] <- TRUE
              lp_files_to_save <- lp_files_names(out_file_path, sub_set, disease, control, status[k], write_p_value, paired)

              ############### All lipid subclass pathways ###############
              valid1 <- FALSE
              
              # Calculate scores or open the file if already exist
              if(mode == "process"){
                json_obj <- convert_scored_lipid_pathway_to_json(lp_data, node_lp_df, disease, control, alt = alt, is_paired, type = sub_set, react = react_keys)
                if(!is.null(json_obj)){
                  save_json_to_file(json_obj, gsub(paste0("_", p_value), "", lp_files_to_save[1]))
                  valid1 <- TRUE
                }
              }else{
                file <- gsub(paste0("_", p_value), "", lp_files_to_save[1])
                if(file.exists(file)){ 
                  json_obj <- gsub(paste0("_", p_value), "", lp_files_to_save[1])
                  valid1 <- TRUE
                }
              }
              
              if(valid1){ 
                dt <- find_all_lipid_sig_pathways(disease_dt, control_dt, p_value , alt = alt, is_paired, type = sub_set, json_obj = json_obj)              
                if(nrow(dt) > 0){
                  dt <- run_highlighted_genes(dt, sub_set, "lipid", lp_files_to_save[1], lp_files_to_save[2])

                  ############### Most significant lipid subclass pathways ###############
                  dt2 <- find_most_sig_pathways(dt, disease_dt, control_dt, alt, is_paired, "lipid")                
                  dt2 <- run_highlighted_genes(dt2, sub_set, "lipid", lp_files_to_save[3], lp_files_to_save[4])
                }
              
            
                ############### All lipid molecular species pathways ###############
                valid2 <- FALSE
                if(mode == "process"){
                  list <- convert_scored_lipid_species_pathway_to_json(node_species_df, disease, control, alt = alt, is_paired, type = sub_set, react = react_keys)
                  if(!is.null(list)){
                    file_to_save <- gsub(paste0("_", p_value), "", lp_files_to_save[5])
                    save(list, file = gsub("json", "RData", file_to_save))
                    json_obj <- convert_pathway_to_json(list$react_df, list$node_df)
                    save_json_to_file(json_obj, file_to_save)
                    valid2<- TRUE
                  }
                }else{
                  file <- paste0(out_file_path, paste0(c("lp_species",sub_set, disease, control, status[k], paired), collapse="_"), ".RData")
                  load(file)
                  if(!is.null(list)){
                    valid2 <- TRUE
                  }
                }
                
                if(valid2){
                  dt <- find_all_lipid_species_sig_pathways(disease_dt, control_dt, p_value , alt = alt, is_paired, type = sub_set, list = list)
                  if(nrow(dt) > 0){
                    dt <- run_highlighted_genes(dt, sub_set, "lipid", lp_files_to_save[5], lp_files_to_save[6])

                    ############### Most significant lipid molecular species pathways ###############
                    dt2 <- find_most_sig_pathways(dt, disease_dt, control_dt, alt, is_paired, 'species')
                    dt2 <- run_highlighted_genes(dt2, sub_set, "lipid", lp_files_to_save[7], lp_files_to_save[8])
                  }
                }
              }
            }
          }

  
          ############### Fatty acid pathways ###############
          if(!is.null(fa_node_df)){
            if(nrow(fa_node_df) > 0){
              valid_fa_reaction <- TRUE
              fa_files_to_save <- fa_files_names(out_file_path, disease, control, status[k], p_value, paired)

              valid3 <- FALSE
              file_to_save <- gsub(paste0("_", p_value), "", fa_files_to_save[1])
              if(mode == "process"){
                json_obj <- convert_scored_fa_pathway_to_json(fa_data, fa_node_df, disease, control, alt, is_paired)
                if(!is.null(json_obj)){
                  save_json_to_file(json_obj, file_to_save)
                  valid3 <- TRUE
                }
              }else{
                if(file.exists(file_to_save)){
                  json_obj <- fromJSON(file_to_save, flatten=TRUE)    
                  valid3 <- TRUE
                }
              }

              if(valid3){
                dt <- find_all_fatty_acid_sig_pathways(fa_disease_dt, fa_control_dt, p_value , alt = alt, is_paired, json_obj)
                if(nrow(dt) > 0){
                  dt <- run_highlighted_genes(dt, "fa", "fa", fa_files_to_save[1], fa_files_to_save[2])

                  ############### Most significant fatty acid pathways ###############
                  dt2 <- find_most_sig_pathways(dt, disease_dt, control_dt, alt, is_paired, 'fa')
                  dt2 <- run_highlighted_genes(dt2, "fa", "fa", fa_files_to_save[3], fa_files_to_save[4])
                }
              }
            }
          }
        }
      }
    }
  }
  
  info <- data.frame(valid_lp_reaction = valid_lp_reaction, valid_fa_reaction = valid_fa_reaction, subset_reaction = subset_reaction$reaction, subset_pathway = subset_reaction$pathway)
  return(info)
}
