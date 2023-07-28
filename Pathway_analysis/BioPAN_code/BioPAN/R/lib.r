########################################
## R file containing shared functions ##
########################################


#### Load libraries #### 
library(jsonlite)
library(stringr)
library(stringdist)
library(dplyr)
library(NMF)
library(graph)
library(RPostgreSQL)


#### Postgres connection ####
con <- dbConnect(PostgreSQL(), user= "biopan", password="access_biopan", dbname="biopan")
right_arrow = "&#8594;"


#### Get the reactions/pathways stored in the BioPAN database ####

# Get the reactions stored in BioPAN with associated genes
get_reaction_gene_df <- function(){
  query <- "SELECT t3.reaction, t3.class, t3.sub_class, t3.type, t3.acyl_add, t3.compound_require, t3.genes 
            FROM (( SELECT t.reaction_id, array_to_string(array_agg(t.gene_symbol), ',') genes 
                    FROM (select r.reaction_id, g.gene_symbol from biopan_reaction r, biopan_gene g, biopan_reaction_gene rg 
                    WHERE r.reaction_id = rg.reaction_id and g.gene_id = rg.gene_id) t 
                    GROUP BY t.reaction_id) t2 
                    LEFT JOIN (select reaction_id, reaction, class, sub_class, type, acyl_add, compound_require from biopan_reaction) t1 ON t2.reaction_id = t1.reaction_id) t3 "
  rs <- dbSendQuery(con,query)
  
  # fetch records from the resultSet into a data.frame
  df <- fetch(rs, n = -1)   # extract all rows
  return(df)
}

# Get only the reactions
get_reaction_df <- function(){
  query <- "SELECT reaction, class, sub_class, type, acyl_add, compound_require
            FROM biopan_reaction"
  rs <- dbSendQuery(con,query)
  
  # fetch records from the resultSet into a data.frame
  df <- fetch(rs, n = -1)   # extract all rows
  return(df)
}


# Get the pathways stored in BioPAN
get_pathway_df <- function(){
  query <- "SELECT DISTINCT tt1.pathway, tt2.name, tt2.descr
            FROM (SELECT t.pathway_batch_id, t.pathway_id, array_to_string(array_agg(t.reaction),'->') pathway 
                  FROM (SELECT t2.pathway_batch_id, t2.pathway_id, t1.reaction
                        FROM biopan_reaction t1, biopan_pathway t2 
                        WHERE t1.reaction_id = t2.reaction_id order by sort_id) t 
                  GROUP BY t.pathway_batch_id, pathway_id) tt1, biopan_pathway tt2
            WHERE tt1.pathway_batch_id = tt2.pathway_batch_id and tt1.pathway_id = tt2.pathway_id "
  rs <- dbSendQuery(con,query)
  
  df <- fetch(rs, n = -1)
  return(df)
}


# Reorder BioPAN pathways
# @pathway_df: dataframe containing pathways stored in the BioPAN database
get_pathway_frame <- function(pathway_df){
  res <- data.frame(class = pathway_df$name, sub_class = pathway_df$descr, pathway = pathway_df$pathway)
  res <- res[order(res$class),]
  return(res)
}


# Return lipid subclasses stored in the BioPAN database
get_database_classes <- function(){
  classes <- c("FaCoA")
  reaction <- get_reaction_gene_df()
  for(i in reaction$reaction){
    is_class <- str_detect(i,"\\(")
    class <- strsplit(i, ",")
    if (!is_class){
      classes <- c(classes,class[[1]][1],class[[1]][2])
    }
    else{
      class1 <- strsplit(class[[1]][1], "\\(")
      class2 <- strsplit(class[[1]][2], "\\(")
      classes <- c(classes,class1[[1]][1],class2[[1]][1])
    }
  }
  return(unique(classes))
}


# Return a dataframe with lipid classes stored in the BioPAN database in upper an lower case
get_database_classes_upper <- function(){
  classes <- database_classes
  lower <- c()
  for(i in 1:length(classes)){
    lower <- c(lower, tolower(classes[i]))
  }
  dt <- data.frame(classes = classes, lower = lower)
  return(dt)
}


# Return a dataframe with the classification of each lipid of the BioPAN database
get_lp_classification <- function(){
  lp <- c()
  class <- c()

  for(i in 1:nrow(lp_reaction)){
    reaction <- lp_reaction[i,"reaction"]
    lipids <- strsplit(reaction,",")
    for(j in 1:length(lipids[[1]])){
      if(!lipids[[1]][j] %in% lp){
        lp <- c(lp, lipids[[1]][j])
        class <- c(class, lp_reaction[i,"class"])
      }
    }
  }

  # add FA
  lp <- c(lp, "FA", "FACOA")
  class <- c(class, "Fatty acid", "Fatty acyl")

  dt <- data.frame(lp = lp, class = class)
  return(dt)
}



#### Store reactions from the BioPAN database ####
reaction <- get_reaction_gene_df()

# Add to the list of reactions the one whithout genes associated
reaction_no_gene <- get_reaction_df()
reaction_no_gene <- reaction_no_gene[which(!reaction_no_gene$reaction %in% reaction$reaction), ]
reaction_no_gene['genes'] <- c(rep(NA, nrow(reaction_no_gene)))
reaction <- rbind(reaction, reaction_no_gene)

pathway <- get_pathway_df()
pathway <- get_pathway_frame(pathway)

lp_reaction <- reaction[which(reaction$type=="lipid"),]
fa_reaction <- reaction[which(reaction$type=="fa"),]

database_classes <- get_database_classes()
database_classes_upper <- get_database_classes_upper()
lp_classification <- get_lp_classification()



# Save json object to file
# @json_object: json object
# @write: file name
save_json_to_file <- function(json_obj, where){
  write(json_obj, where)
}


# Append array to list
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}


# Get a list of lipid subclass containing the associated molecular species
# @species_names: list of species
# @return: list
get_species_class <- function(species_names){
  pattern <- "^(\\w+)\\([^()]+\\)" 
  if(length(species_names) == 1){
    class <- strsplit(as.character(species_names),"\\(")[[1]][1]
    return(class)
  }
  s_class <- list()
  for(i in 1:length(species_names)){
    s <- as.character(species_names[i])
    class <- strsplit(s,"\\(")[[1]][1]
    # if the list is empty then add element
    class_names <- names(s_class)
    if(is.null(class_names)){
      s_class[[class]] <- s
    }else if(class %in% class_names){
      # if the class has some elements then just add the new one
      s_class[[class]] <- c(s_class[[class]] , s)
    }else{
      # if the group has not been created then create it and add new element
      s_class[[class]] <- s
    }
  }
  return(s_class)
}


# Extract acyl chain structure
# @s: lipid species
get_acyl_chain <- function(s){
  if(is.vector(s)){
    fas <- c()
    for(i in 1:length(s)){
      if(str_count(s[i],"\\)") == 1){
        acyl <- str_extract_all(s[i], "\\([^()]+\\)")[[1]]
        acyl <- substring(acyl, 2, nchar(acyl)-1)
      }else{
        acyl <- strsplit(s[i],"\\(")[[1]][2]
      }
      fas <- c(fas, acyl) 
    }
    return(fas)
  }else{
    if(str_count(s[i] == 1)){
      acyl <- str_extract_all(s, "\\([^()]+\\)")[[1]]  
      acyl <- substring(acyl, 2, nchar(acyl)-1)
    }else{
      acyl <- strsplit(s[i],"\\(")[[1]][2]
    }
   return(acyl)
  }
}
