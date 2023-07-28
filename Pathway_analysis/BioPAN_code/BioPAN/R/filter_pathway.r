###################################################################################################
## Filter pathways according to the species selected by the user to only visualise these species ##
## Store the results in a "filter" file                                                          ##
###################################################################################################


source("/lipidmaps/biopan/R/lib.r")


args <- commandArgs(TRUE)
# retrieve the values from the filter_pathway.php file
out_file_path <- args[1]
path <- args[2]


# file containing the selected molecular species
selected_file <- paste0(out_file_path, "selected_items.csv")

# to check if the file is empty
filecontent <- lapply(selected_file, function(x) {
    tryCatch(read.csv(x, header=F, sep = ','), error=function(e) NULL)
})

if(!is.null(filecontent[[1]])){
    selected_items <- read.csv(selected_file, header=F, sep = ',')
    selected_items <- gsub("\\(|\\)|:", "", tolower(as.character(t(selected_items)[,1])))
    json_obj <- fromJSON(path, flatten=TRUE)

    nodes <- list()
    edges <- list()
    j = 1 #nodes counter
    k = 1 #edges counter

    for(i in 1:nrow(json_obj$nodes)){
        if(json_obj$nodes[i,"data.id"] %in% selected_items){
            node <- list(data = list(id = json_obj$nodes[i,"data.id"], lm_id = json_obj$nodes[i,"data.lm_id"], label = json_obj$nodes[i,"data.label"], name = json_obj$nodes[i,"data.name"], shape = json_obj$nodes[i,"data.shape"]))
            nodes[[j]] <- node
            j <- j + 1
        }
    }

    for(i in 1:nrow(json_obj$edges)){
        if(json_obj$edges[i,"data.source"] %in% selected_items && json_obj$edges[i,"data.target"] %in% selected_items){
            edge <- list(data = list(id = json_obj$edges[i,"data.id"], source = json_obj$edges[i,"data.source"], target = json_obj$edges[i,"data.target"], weight = json_obj$edges[i,"data.weight"], color = json_obj$edges[i,"data.color"]))
            edges[[k]] <- edge
            k <- k + 1
        }
    }

    gr <- list(nodes = nodes, edges = edges)
    json_obj_update <- toJSON(gr, pretty=T, auto_unbox=T)

    file_name <- paste0(gsub(".json", "", path), "_filter.json")
    save_json_to_file(json_obj_update, file_name)
}
