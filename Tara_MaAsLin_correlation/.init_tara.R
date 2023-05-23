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
  # install.packages("ggVennDiagram")
  library("ggVennDiagram")
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
  getwd()
}

#from: Community-Level Responses to Iron Availability in Open Ocean Plankton Ecosystems
iron_f="../metadata/Tara-iron_gbc20832-sup-0009-table_s01_downloaded-April-30-2023.xlsx"

gene_abundance_file = "../data/Tara_139_COGonly.tsv" # download April 12th 2023
metadata_file = "../metadata/OM.CompanionTables-downloaded-April11-2022.xlsx"

fudge_df = function(df) {
  fudge_factor = min(df[df>0])/2
  new_df = df + fudge_factor
  # new_df = new_df[apply(new_df, 1, var) != 0, ]
  return(list(new_df, fudge_factor))
}

fudge_vec = function(vec) {
  tmp = vec[!is.na(vec)]
  tmp = tmp[tmp > 0]
  return(vec + min(tmp)/2)
}

remove_special_chars <- function(x) {
  x = gsub("\\[.*?\\]", "", x)
  x = gsub("\\*", "", x)
  x = trimws(x, "right")
  return(x)
}

# Retrieval and integration of metadata
init_tara_metadata = function(){
  init_env()
  md = read_excel(metadata_file, sheet = 9, col_names = TRUE, col_types = NULL, na = "NA", skip = 0)
  linker = read_excel(metadata_file, sheet = 2, col_names = TRUE, skip = 0)[c(1,5)]
  colnames(linker) = c("sample", "PANGAEA Sample ID")
  
  linker %>% filter(! str_detect(sample, "<"))
  
  colnames(md) = sapply(colnames(md), remove_special_chars)
  
  md2 = merge(linker, md, by = "PANGAEA Sample ID")
  
  md2$tmp = gsub("TARA_", "", md2$sample)
  md2$tmp = gsub("^0+", "", md2$tmp)

  md2 = md2 %>%
    separate(tmp, into = c("station", "layer", "size_fraction"), sep = "_") %>%
    mutate(size_fraction = str_replace(size_fraction, "<", "0.0"),
           log_PO4 = log10(fudge_vec(PO4)),
           log_NO2 = log10(fudge_vec(NO2)),
           log_NO2NO3 = log10(fudge_vec(NO2NO3)),
           log_depth = log10(Mean_Depth),
           log_FC_heterotrophs = log10(`FC - heterotrophs`),
           log_FC_bacteria = log10(`FC - bacteria`),
           Mean_Nitrates = ifelse(Mean_Nitrates > 0, Mean_Nitrates, 0),
           log_Nitrates = log10(fudge_vec(Mean_Nitrates))) 
  
  iron = get_iron_concentration()
  md3 = merge(iron, md2, by = c("station", "layer"), all.y = TRUE)  
  
  md3 = md3 %>%
    mutate(sqrt_iron = sqrt(PISCES2),
           size_fraction = as.factor(size_fraction),
           layer = factor(layer, levels=c('MIX', 'MES', 'DCM', 'SRF'))) %>%
    remove_rownames %>% 
    column_to_rownames(var="sample")
  
  return(md3)
}

get_iron_concentration = function() {
  Tara_iron=read_excel(iron_f, skip = 4, na = "NA",
                       col_names = c("station", "layer", "PISCES2", "ECCO2-DARWIN",
                      "Mean_Lat", "Mean_Long", "observed_iron", "Resolution"))
  Tara_iron$layer = gsub("SUR", "SRF", Tara_iron$layer)
  # most of the ECCO2-DARWIN were NA for the Tara samples we have
  
  return(Tara_iron[c("station", "layer", "PISCES2", "ECCO2-DARWIN", "observed_iron", "Resolution")])
}


get_COGabundance = function(to_fudge = TRUE) {
  COGabun = read.csv(gene_abundance_file, sep = "\t")
  colnames(COGabun) = str_replace(str_replace(colnames(COGabun), "\\.3", "-3"), "\\.1.6", "-1.6")
  return(COGabun)
}

get_significant_res = function(res_path, sep_by, strict = 0.05) {
  res = paste(res_path,"significant_results.tsv",sep="/")
  feature_sigs = read_delim(file=res, delim = "\t", col_types = cols())
  
  if (sep_by[1] != "all") { feature_sigs = filter(feature_sigs, metadata %in% c(sep_by)) }
  
  feature_sigs_list = feature_sigs %>% filter(qval < strict) %>% arrange(coef, qval) %>% pull(feature)
  return(feature_sigs_list)
}

resolve_COG_function = function(vec) {
  COG_dict = read_excel(metadata_file, sheet = 8, col_names = TRUE, col_types = NULL, na = "NA", skip = 0)
  COG_dict = COG_dict[c('COG', 'Name')] %>% 
    remove_rownames %>% 
    column_to_rownames(var='COG')
  return(COG_dict[vec, ])
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

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
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


