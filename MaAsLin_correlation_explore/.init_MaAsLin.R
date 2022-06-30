init_env <- function(){
  library(tidyverse)
  library(dplyr)
  library(readxl)
  library(readr)
  library(stringr)
  # library(pheatmap) # for drawing heat maps
  library(RColorBrewer)
  library(viridis)
  library(Maaslin2)
  # library(vegan) # for Principle Coordinate Analysis
  library(ggplot2)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
  getwd()
}

gene_abundance_file = "../data/Ustick2021-some-nutrient-genes-coverage-per-sample.csv"
metadata_file = "../organize_metadata/Bio-GO-SHIP_sample_metadata.csv"

fudge_df = function(df) {
  fudge_factor <- min(df[df>0])/2
  new_df = df + fudge_factor
  new_df = new_df[apply(new_df, 1, var) != 0, ]
  return(list(new_df, fudge_factor))
}

# Retrieval and integration of metadata
init_BGS_metadata = function(only_phosphate = FALSE){
  init_env()
  md = read_csv(metadata_file)
  geneabun = read_csv(gene_abundance_file)
  geneabun_md = geneabun[c("SRA_Accession", "Cruise", "Depth_m", "Month", "Local_Hour")]
  md = merge(md, geneabun_md, by="SRA_Accession")

  md = md %>% remove_rownames %>% column_to_rownames(var="SRA_Accession")
  
  # round up local hours to hourly, then cut a day into 6 periods
  md$time_in_day = round(as.numeric(md$Local_Hour))
  md$time_in_day = cut(md$time_in_day, breaks = c(0, 4, 8, 12, 16, 20, 24),
                       labels = c("0to4", "4to8", "8to12", "12to16","16to20","20to24"))
  
  if (only_phosphate) {
    md = filter(md, !is.na(phosphate))
    p_fudge_factor = min(md$phosphate[md$phosphate > 0])/2
    # n_fudge_factor = min(md$nitro[md$nitro > 0])/2 # too much NAs
    md$log_phosphate = log10(md$phosphate + p_fudge_factor)
    # md$log_nitro = log10(md$nitro + n_fudge_factor)
  }
  return(md)
}

# metadata = init_BGS_metadata()

get_geneabundance = function(metadata, to_fudge = TRUE) {
  geneabun = read_csv(gene_abundance_file)
  geneabun = geneabun %>% filter(SRA_Accession %in% rownames(metadata))
  geneabun = geneabun[-c(2:13,98:106)] # get rid of useless data
  geneabun = geneabun %>% remove_rownames %>% column_to_rownames(var="SRA_Accession")
  
  if (to_fudge) {
    return(fudge_df(geneabun)) # df + fudge_factor
  } else {
    return(geneabun)
  }
}
get_significant_res <- function(res_path, sep_by, strict = 0.25) {
  feature_sigs <- read_delim(file = paste(res_path,"significant_results.tsv",sep="/"), delim = "\t")
  
  if (sep_by[1] != "all") { feature_sigs = filter(feature_sigs, metadata %in% c(sep_by)) }
  
  feature_sigs_list = feature_sigs %>% filter(qval < strict) %>% arrange(coef, qval) %>% pull(feature)
  
  feature_sigs_list = gsub("X1CMET2","1CMET2", fixed=TRUE, feature_sigs_list)
  return(feature_sigs_list)
}

simplify_abundance_df <- function(abundance_df){
  simple_rn=gsub(pattern = "(", replacement = ".", fixed = TRUE, rownames(abundance_df))
  
  for (k in c(")","|",":"," ","-",",","&",";","'","+")) {
    simple_rn=gsub(k,".",fixed=TRUE,simple_rn)
  }
  
  abundance_df$new_rowname = simple_rn
  s_abundance_df = abundance_df %>% remove_rownames %>% column_to_rownames(var="new_rowname")
  return(s_abundance_df)
}

gen_df_for_graphing <- function(abundance_df, res_path, sep_by, 
                                strict = 0.25, return_sorted = FALSE) {
  #import sig features
  path_sigs_Feats = get_significant_res(res_path, sep_by, strict)
  abundance_df = simplify_abundance_df(abundance_df)
  
  msg = "Err (see init_share.R), not all significant features are plotted. Setdiff should return nothing."
  if (length(setdiff(path_sigs_Feats,rownames(abundance_df)) > 0)) print(msg)
  
  #Select features most significantly different between PBH and everyone else
  graph_df = abundance_df %>% filter(rownames(abundance_df) %in% path_sigs_Feats)
  
  if (return_sorted) {
    return(list(path_sigs_Feats, graph_df))
  } else {
    return(graph_df) 
  }
}

gen_heatmap = function(graph_df, metadata, outpath, g_width = 5, exclude_features = c(),
                       g_height=10, cluster_row = FALSE, annotation_colors = NA) {
  
  graph_df = graph_df[!row.names(graph_df) %in% exclude_features, ]
  
  #Plot pheatmap
  pheatmap(mat = graph_df[rownames(metadata)],
           annotation_col = metadata,
           cluster_cols = F,
           cluster_rows = cluster_row,
           show_rownames = T,
           show_colnames = F,
           annotation_colors = annotation_colors,
           cellwidth = g_width, cellheight = g_height, fontsize = 8, 
           filename = outpath, #"./analysis_graph/MaAsLin_sig_pheatmap_pathway_univar.png"
           scale = "none"
  )
  # dev.off()
}

# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# code taken from https://stackoverflow.com/questions/7519790/assign-multiple-new-variables-on-lhs-in-a-single-line
# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}


