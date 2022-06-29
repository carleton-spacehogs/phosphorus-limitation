init_env <- function(){
  library(tidyverse)
  library(dplyr)
  library(readxl)
  library(readr)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(Maaslin2)
  library(vegan)
  library(ggplot2)
  library("textclean")
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
  getwd()
}

# Retrieval and integration of metadata
init_BGS_metadata <- function(){
  init_env()
  md = read_csv("../organize_metadata/Bio-GO-SHIP_sample_metadata.csv")
  geneabun = read_csv("../data/Ustick2021-some-nutrient-genes-coverage-per-sample.csv")
  geneabun_md = geneabun[c("SRA_Accession", "Cruise", "Depth_m_")]
  
  md = merge(md, geneabun_md, by="SRA_Accession")
  md$log_phosphate = log(md$phosphate)
  
  md = md %>% remove_rownames %>% column_to_rownames(var="SRA_Accession")
  return(md)
}

get_geneabundance = function() {
  md = init_BGS_metadata()
  geneabun = read_csv("../data/Bio-Go-nutrient-genes-coverage-per-sample.csv")
  geneabun2 = geneabun %>% filter(SRA_Accession %in% rownames(md))
  geneabun2 = geneabun2 %>% remove_rownames %>% column_to_rownames(var="SRA_Accession")
  
  geneabun2 = geneabun2[-c(1:12,97:105)]
  
  nudge_factor <- min(geneabun2[geneabun2>0])/2
  geneabunN = geneabun2 + nudge_factor
  geneabunN = geneabunN[apply(geneabunN, 1, var) != 0, ]
  return(list(geneabunN, nudge_factor))
}

get_significant_res <- function(res_path, sep_by) {
  feature_sigs <- read_delim(file = paste(res_path,"significant_results.tsv",sep="/"), delim = "\t")
  #[1] "X1CMET2.PWY..folate.transformations.III..E..coli..g__Butyrivibrio.s__Butyrivibrio_crossotus"
  # manually remove the "X" before "X1C..." didn't know where the X came from.
  
  if (sep_by != "all") {
    feature_sigs_list <- feature_sigs %>%
      filter(metadata == sep_by) %>%
      arrange(coef, qval) %>%
      pull(feature)
  } else {
    feature_sigs_list <- feature_sigs %>%
      arrange(coef, qval) %>%
      pull(feature)
  }
  
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

gen_df_for_graphing <- function(abundance_df, res_path, sep_by) {
  #import sig features
  path_sigs_Feats = get_significant_res(res_path, sep_by)
  abundance_df = simplify_abundance_df(abundance_df)
  
  msg = "Err (see init_share.R), not all significant features are plotted. Setdiff should return nothing."
  if (length(setdiff(path_sigs_Feats,rownames(abundance_df)) > 0)) print(msg)
  
  #Select features most significantly different between PBH and everyone else
  graph_df = abundance_df %>% filter(rownames(abundance_df) %in% path_sigs_Feats)
  
  return(graph_df)
}

gen_heatmap = function(graph_df, metadata, outpath, g_width = 5, g_height=10) {
  #Plot pheatmap
  pheatmap(mat = graph_df[rownames(metadata)],
           annotation_col = metadata,
           cluster_cols = F,
           cluster_rows = F,
           show_rownames = T,
           show_colnames = F,
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


