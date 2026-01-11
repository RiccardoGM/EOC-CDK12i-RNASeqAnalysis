######################################
#            Libraries               #
######################################

message("Loading required libraries...") 

{
  library(stringr)
  library(EnsDb.Hsapiens.v79)
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(scales)
  library(jsonlite)
  library(ashr)
  library(rtracklayer)
  library(patchwork)
  source("utils.R")
}


######################################
#            Set Seed                #
######################################

seed_number <- 1234
set.seed(seed_number)


######################################
#            Load Data               #
######################################

message("Loading data and metadata...")

{
  #--- Read parameters ---#
  params        <- fromJSON("./parameters.json")
  line_name     <- params$line_name
  path_data     <- params$counts_dir
  path_metadata <- params$metadata_dir
  path_results  <- paste(params$results_dir, "/", line_name, "/", sep="")
  path_figures  <- paste(params$figures_dir, "/", line_name, "/", sep="")
  gene_type     <- params$gene_type
  path_gtf      <- params$gtf_dir
  n_replicates  <- params$n_replicates
  dge_params    <- params$DGE
  gsea_params   <- params$gsea
  
  print(sprintf("Line: %s | Gene type: %s | Replicates: %d", line_name, gene_type, n_replicates))
  
  #--- Read contrasts ---#
  contrasts_json <- fromJSON(params$contrasts_file, simplifyVector=FALSE)
  contrasts <- lapply(contrasts_json$contrasts, unlist)

  #--- Create directories ---#
  if (!dir.exists(path_results)) {dir.create(path_results, recursive=TRUE, showWarnings=FALSE)}
  if (!dir.exists(path_figures)) {dir.create(path_figures, recursive=TRUE, showWarnings=FALSE)}
  
  #--- Load metadata ---#
  file_name   <- sprintf("%s_metadata.csv", line_name)
  path        <- paste(path_metadata, file_name, sep="/")
  metadata    <- read.csv(path, header=T, stringsAsFactors=F, row.names=1, sep=";") 
  sample_name <- metadata$name
  
  #--- Load counts matrix ---#
  file_name      <- sprintf("%s_counts_matrix.txt", line_name)
  path           <- paste(path_data, file_name, sep="/")
  data           <- read.table(path, header=F, row.names=1)
  colnames(data) <- sample_name
  
  #--- Load reference GTF ---#
  gtf_data <- import(path_gtf)
  gtf_df   <- as.data.frame(gtf_data)
}


######################################
#          Prepare Data              #
######################################

message("Preparing data: filtering and normalizing genes...")

{
  #--- Remove version from Ensembl IDs ---#
  genes_ID <- rownames(data)
  genes_ID_noversion <- stringr::str_replace(genes_ID, pattern=".[0-9]+$", replacement="")
  rownames(data) <- genes_ID_noversion
  
  #--- Set factors ---#
  metadata$treatment <- factor(metadata$treatment)
  metadata$baseline  <- factor(metadata$baseline)
  
  #--- DESeqDataSet ---#
  colData <- metadata
  rownames(colData) <- colData$name
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=data[,1:(ncol(data))], 
                                        colData=colData, design=~1)
  
  #--- Filter genes by type ---#
  # Protein coding genes
  protein_coding_genes_ensemblid <- gtf_df[gtf_df$gene_type == "protein_coding", "gene_id"]
  protein_coding_genes_ensemblid <- stringr::str_replace(protein_coding_genes_ensemblid, 
                                                         pattern=".[0-9]+$", replacement="")
  protein_coding_genes_ensemblid <- unique(protein_coding_genes_ensemblid)
  
  # lncRNA genes
  lncRNA_genes_ensemblid <- gtf_df[gtf_df$gene_type == "lncRNA", "gene_id"]
  lncRNA_genes_ensemblid <- stringr::str_replace(lncRNA_genes_ensemblid, pattern=".[0-9]+$", replacement="")
  lncRNA_genes_ensemblid <- unique(lncRNA_genes_ensemblid)

  # Filter based on gene_type parameter
  dds_genes_ensemblid <- rownames(dds)
  if (gene_type == "lncRNA_ProteinCoding") {
    selected_gene_indices <- match(c(protein_coding_genes_ensemblid, lncRNA_genes_ensemblid), dds_genes_ensemblid)
  } else if (gene_type == "ProteinCoding") {
    selected_gene_indices <- match(protein_coding_genes_ensemblid, dds_genes_ensemblid)
  } else {
    selected_gene_indices <- dds_genes_ensemblid
  }
  selected_gene_indices <- selected_gene_indices[!is.na(selected_gene_indices)]
  dds <- dds[selected_gene_indices, ]
  
  #--- Normalize counts ---#
  dds <- DESeq2::estimateSizeFactors(dds)
  normalized_counts <- DESeq2::counts(dds, normalized=T)
  
  #--- Filter low-expressed genes ---#
  mask <- 0
  n_conditions <- dim(normalized_counts)[2]/n_replicates
  for (cond in 1:n_conditions){
    mask <- mask + as.integer(rowSums(normalized_counts[,cond:(cond+n_replicates-1)] >= 10) > 1)
  }
  mask         <- mask > n_conditions/2.
  keep         <- mask
  dds_filtered <- dds[keep,]
  
  #--- Selected gene lists ---#
  filtered_gene_indices         <- rownames(dds_filtered)
  filtered_coding_genes_ids     <- intersect(filtered_gene_indices, protein_coding_genes_ensemblid)
  filtered_non_coding_genes_ids <- setdiff(filtered_gene_indices, protein_coding_genes_ensemblid)
  
  #--- Print counts ---#
  message(sprintf("Total filtered genes: %d", length(filtered_gene_indices)))
  message(sprintf("Protein-coding genes: %d", length(filtered_coding_genes_ids)))
  message(sprintf("Non-coding genes: %d", length(filtered_non_coding_genes_ids)))
}


######################################
#         PCA for All Samples        #
######################################

message("Performing PCA on normalized data...")

{
  #--- PCA ---#
  data_mat           <- SummarizedExperiment::assay(dds_filtered)
  rlog_data          <- DESeq2::rlog(dds_filtered, blind=T)
  rlog_data_mat      <- SummarizedExperiment::assay(rlog_data)
  mask               <- rowVars((rlog_data_mat), useNames=TRUE)>0
  pca                <- prcomp(t(rlog_data_mat[mask,]), scale.=TRUE)
  variance_explained <- round((pca$sdev^2 / sum(pca$sdev^2))*100, digits=0)
  
  #--- Prepare plot dataframe ---#
  df_pca     <- pca$x
  conditions <- gsub("\\.", "_", interaction(metadata$baseline, metadata$treatment))
  conditions <- ifelse(grepl("_2024", metadata$name), paste0(conditions, "_2024"), conditions)
  
  {
    # To customize baseline colors:
    #baseline_colors_dict <- list(baseline_level1=color_cerise, 
    #                             baseline_level2=color_yellow, 
    #                             baseline_level3=color_palblue)
    #baseline_colors <- get_colors(metadata$baseline, baseline_colors_dict)
    #color_map <- setNames(baseline_colors, metadata$baseline)
    
    #--- Plot PCA ---#
    df <- data.frame(x=df_pca[, "PC1"], y=df_pca[, "PC2"], 
                     baseline=metadata$baseline, treatment=metadata$treatment)
    p <- ggplot(df) +
      geom_point(aes(x=x, y=y, color=baseline, shape=treatment), 
                 size=3.5, alpha=0.65) +
      #scale_color_manual(values=color_map) +
      theme_classic() +  
      labs(color="Baseline", shape="Treatment",
           x=sprintf("PC1: %%%.0f variance", variance_explained[1]), 
           y=sprintf("PC2: %%%.0f variance", variance_explained[2])) + 
      ggtitle(sprintf("line: %s", line_name)) + common_theme
    
    figure_name <- sprintf("PCA_%s.jpeg", line_name)
    figure_path <- file.path(path_figures, figure_name)
    ggsave(figure_path, plot=p, device='jpeg', scale=1, width=10, height=8, 
           units='cm', dpi=500, create.dir=TRUE)
  }
}

message("PCA plot saved.")


######################################
#            DGE (DESeq2)            #
######################################

message("Running Differential Gene Expression analysis...")

#--- Validate design ---#
{
  allowed_factors <- metadata_columns
  design_factors  <- dge_params$design
  
  if (length(design_factors) < 1 || length(design_factors) > 2) {
    stop("design_factors must contain one or two elements.")
  }
  
  if (!all(design_factors %in% allowed_factors)) {
    stop(
      sprintf(
        "Invalid design factor(s): %s. Allowed values are: %s",
        paste(setdiff(design_factors, allowed_factors), collapse = ", "),
        paste(allowed_factors, collapse = ", ")
      )
    )
  }
  
  if (length(design_factors) == 2 &&
      !setequal(design_factors, allowed_factors)) {
    stop("When specifying two design factors, they must be 'baseline' and 'treatment'.")
  }
  message("   >>> Design validated: ", paste(design_factors, collapse = " and "))
}

#--- Set design ---#
{
  if (length(design_factors) == 1) {
    
    design_factor        <- design_factors[1]
    design_name          <- design_factor
    design(dds_filtered) <- as.formula(paste("~", design_factor))
    available_levels     <- levels(factor(colData(dds_filtered)[[design_factor]]))
    
  } else {
    design_name <- "baseline_and_treatment"
    dds_filtered[[design_name]] <- factor(paste0(dds_filtered$baseline, 
                                                 dds_filtered$treatment))
    design(dds_filtered) <- as.formula(paste("~", design_name))
    available_levels     <- levels(dds_filtered[[design_name]])
  }
  message("   >>> Single factor levels: ", paste(available_levels, collapse = ", "))
}

#--- Validate contrasts ---#
{
  if (!all(sapply(contrasts, length) == 3)) {
    stop("Each contrast must be a character vector of length 3: c(design_name, level_to_compare, base_level).")
  }
  
  for (contrast in contrasts) {
    level_to_compare <- contrast[2]
    base_level       <- contrast[3]
    
    missing_levels <- setdiff(c(level_to_compare, base_level), available_levels)
    if (length(missing_levels) > 0) {
      stop(sprintf("Invalid contrast levels: %s. Available levels are: %s",
                   paste(missing_levels, collapse = ", "),
                   paste(available_levels, collapse = ", ")))
    }
  }
  message("   >>> All contrasts validated successfully.")
}

#--- Run DGE analysis ---#
{
  results_dict <- list()
  dds_filtered <- DESeq(dds_filtered)

  for (contrast in contrasts) {
    
    contrast_design  <- contrast[1]
    level_to_compare <- contrast[2]
    base_level       <- contrast[3]
    
    if (contrast_design != design_name) {
      stop(
        sprintf("Contrast design '%s' does not match active design '%s'.", contrast_design, design_name))
    }
    
    if (level_to_compare %in% available_levels &&
        base_level %in% available_levels) {
      
      contrast_name <- sprintf("%svs%s", level_to_compare, base_level)
      res           <- results(dds_filtered, contrast=contrast)
      res_shrunk    <- lfcShrink(dds_filtered, res=res, type="ashr")
      results_dict[[contrast_name]] <- res_shrunk
    }
  }
}

message("DGE analysis complete. Results stored in results_dict.")


######################################
#           LFC vs Gene Length       #
######################################

message("Generating LFC vs gene length plots...")

{
  #--- Prepare gene widths ---#
  gene_df                <- gtf_df[gtf_df$type == "gene", ]
  genes_ID               <- gene_df$gene_id
  genes_ID_noversion     <- stringr::str_replace(genes_ID, pattern=".[0-9]+$", replacement="")
  rownames(gene_df)      <- genes_ID_noversion
  coding_gene_widths     <- gene_df[filtered_coding_genes_ids, ]$width
  non_coding_gene_widths <- gene_df[filtered_non_coding_genes_ids, ]$width
  
  #--- Loop through contrasts ---#
  for (name in names(results_dict)){
    plot1 <- plot_lfc_vs_length(ensembl_ids=filtered_coding_genes_ids,
                                gene_len=coding_gene_widths,
                                comparison=name, 
                                results_dict=results_dict, 
                                line_name=line_name,
                                gene_type="prontein-coding")
    
    noncoding_gene_type <- "non-coding"
    if (gene_type=="lncRNA_ProteinCoding"){
      noncoding_gene_type <- "long non-coding"
    }
    
    plot2 <- plot_lfc_vs_length(ensembl_ids=filtered_non_coding_genes_ids,
                                gene_len=non_coding_gene_widths,
                                comparison=name, 
                                results_dict=results_dict, 
                                line_name=line_name,
                                gene_type=noncoding_gene_type)
    
    #--- Combine plots ---#
    combined_plot <- plot1 + plot2 + plot_layout(ncol=2) +
                     plot_annotation(title=sprintf("%s - %s", line_name, name))
    parameters_description <- sprintf("line#%s", line_name)
    figure_name <- sprintf("GeneslfcVSlength__contrast#%s__%s.jpeg", name, parameters_description)
    figure_path <- sprintf("%s%s/%s", path_figures, name, figure_name)
    ggsave(figure_path, plot=combined_plot, device='jpeg', scale=1, width=15, 
           height=6, units='cm', dpi=500, create.dir=TRUE)
  }
}

message("LFC vs gene length plots completed.")


######################################
#     GO GSEA (ClusterProfiler)      #
######################################

message("Running GSEA analysis using clusterProfiler...")

#--- Read parameters ---#
{
  ranking_list   <- gsea_params$ranking_list
  exponent_list  <- gsea_params$exponent_list
  cutbypval_list <- gsea_params$cutbypval_list
  keyType        <- gsea_params$keyType
  nPermSimple    <- gsea_params$nPermSimple
  minGSSize      <- gsea_params$minGSSize
  maxGSSize      <- gsea_params$maxGSSize
  pvalueCutoff   <- gsea_params$pvalueCutoff
  padjTh         <- gsea_params$padjTh
  pAdjustMethod  <- gsea_params$pAdjustMethod
  ont            <- gsea_params$ont
  seed           <- gsea_params$seed
  eps            <- gsea_params$eps
  OrgDb          <- org.Hs.eg.db
  string_parameters <- sprintf('nPermSimple#%d_minGSSize#%d_maxGSSize#%d_pAdjustMethod#%s', 
                               nPermSimple, minGSSize, maxGSSize, pAdjustMethod)
  selected_go_terms <- gsea_params$selected_go_terms
}

#--- Print parameters ---#
{
  print(sprintf("Ranking stat list: %s", ranking_list))
  print(sprintf("NES exponent list: %s", exponent_list))
  print(sprintf("Cut by p-val list: %s", cutbypval_list))
  print(sprintf("P-val adj method:  %s", pAdjustMethod))
  print(sprintf("Min gene set size: %s", minGSSize))
  print(sprintf("Max gene set size: %s", maxGSSize))
  print(sprintf("N. permutations:   %s", nPermSimple)) 
}

#---  Run gseGO analyses  ---#
{
  GOgse_results <- list()
  
  #--- Loop over p-value cutoff ---#
  for (cutbypval in cutbypval_list){
    cutbypval_str <- sprintf('cutbypval#%s',cutbypval)
    print(cutbypval_str)
    GOgse_results[[cutbypval_str]] <- list()
    
    #--- Loop over exponent ---#
    for (exponent in exponent_list){
      exponent_str <- sprintf('exponent#%d',exponent)
      print(exponent_str)
      GOgse_results[[cutbypval_str]][[exponent_str]] <- list()
      
      #--- Loop over ranking list ---#
      for (ranking in ranking_list){
        print(ranking)
        GOgse_results[[cutbypval_str]][[exponent_str]][[ranking]] <- list()
        
        #--- Loop over contrasts ---#
        for (name in names(results_dict)){
          print(name)
          
          #--- Prepare ranking stat ---#
          res_local  <- results_dict[[name]]
          gene_names <- rownames(res_local)
          if (ranking=='ranking#lfc'){
            stat_table <- res_local$log2FoldChange
          }
          else if (ranking=='ranking#SlfcLogP'){
            pvalue_table      <- res_local$pvalue
            mask_pvalue_table <- (pvalue_table>0)&!is.na(pvalue_table)
            min_pvalue_table  <- min(pvalue_table[mask_pvalue_table])
            pvalue_table[pvalue_table==0] <- min_pvalue_table / 10.
            stat_table <- sign(res_local$log2FoldChange)*(-log10(pvalue_table))
          }
          else {
            stop(paste("Invalid ranking option:", ranking))
          }
          
          names(stat_table) <- gene_names
          
          #--- Mask by adjusted p-value ---#
          if (cutbypval){
            padj_values <- res_local$padj
            mask_padj   <- padj_values<padjTh
            stat_table  <- stat_table[mask_padj]
          }
          
          #--- Remove NaNs / infinite values ---#
          mask_notnull <- (!is.na(stat_table)) & is.finite(stat_table)
          stat_table   <- stat_table[mask_notnull]
          
          #--- Add tiny noise & sort ---#
          stat_table <- stat_table + rnorm(length(stat_table), mean=0, sd=1e-6)
          stat_table <- sort(stat_table, decreasing=TRUE)
          
          #--- Run GSEA ---#
          gsego <- gseGO(stat_table,
                         keyType=keyType, 
                         nPermSimple=nPermSimple, 
                         minGSSize=minGSSize,
                         maxGSSize=maxGSSize, 
                         pvalueCutoff=pvalueCutoff,
                         exponent=exponent,
                         pAdjustMethod=pAdjustMethod,
                         OrgDb=OrgDb,
                         ont=ont,
                         seed=seed,
                         eps=eps)
          
          #--- Save result ---#
          GOgse_results[[cutbypval_str]][[exponent_str]][[ranking]][[name]] <- gsego
        }
      }
    }
  }
}

message("GSEA analysis completed.")


######################################
#     Running Enrich. Score Plots    #
######################################

message("Generating running enrichment score plots...")

{
  print("----------- Producing running enrichment score plots -----------")
  
  for (cutbypval in cutbypval_list){
    cutbypval_str <- sprintf('cutbypval#%s',cutbypval)
    
    for (exponent in exponent_list){
      exponent_str <- sprintf('exponent#%d',exponent)
      
      for (ranking in ranking_list){
        
        for (name in names(results_dict)){
          
          for (pathway_name in selected_go_terms) {
            # Running enrichment plot
            gsego_local <- GOgse_results[[cutbypval_str]][[exponent_str]][[ranking]][[name]]
            GO_term_idx <- match(c(pathway_name), gsego_local$Description)
            if (!is.na(GO_term_idx)) {
              print(paste(pathway_name, "- enrichment score:", gsego_local[GO_term_idx,]$enrichmentScore,
                          "- NES:", gsego_local[GO_term_idx,]$NES))
              GO_term_ID        <- gsego_local[GO_term_idx,]$ID
              GO_term_NES       <- gsego_local[GO_term_idx,]$NES
              print(sprintf('Size set: %d', length(gsego_local[[GO_term_ID]])))
              title             <- sprintf("%s - %s - %s, NES=%.2f, exp:%d", line_name, name, pathway_name, GO_term_NES, exponent)
              runningscore_plot <- gseaplot(gsego_local, geneSetID=GO_term_idx, by="runningScore", title=pathway_name) + 
                #coord_cartesian(ylim=c(-0.4, 0.4)) +
                ggtitle(title)
              parameters_description <- sprintf("line#%s__%s__%s__%s", line_name, ranking, exponent_str, cutbypval_str)
              pathway_name_formatted <- gsub(" ", "", pathway_name)
              figure_name            <- sprintf("RunningScoreGOterm#%s__contrast#%s__%s.jpeg", pathway_name_formatted, name, parameters_description)
              figure_dir             <- sprintf("%s%s/running_enrichment_score/GOterms/", path_figures, name)
              dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
              figure_path <- sprintf("%s/%s", figure_dir, figure_name)
              ggsave(figure_path, plot=runningscore_plot, device='jpeg', scale=1, 
                     width=18, height=9, units='cm', dpi=500, create.dir=TRUE)
            }
          }
        }
      }
    }
  } 
}

message("Running enrichment score plots saved.")


######################################
#        GO terms categories         #
######################################

message("Grouping and exporting GO terms by category...")

#--- Group GO terms into categories ---#
{
  keyword_pathway_dict <- list()
  for (cutbypval in cutbypval_list){
    cutbypval_str <- sprintf('cutbypval#%s',cutbypval)
    keyword_pathway_dict[[cutbypval_str]] <- list()
    
    for (exponent in exponent_list){
      exponent_str <- sprintf('exponent#%d',exponent)
      keyword_pathway_dict[[cutbypval_str]][[exponent_str]] <- list()
      
      for (ranking in ranking_list){
        keyword_pathway_dict[[cutbypval_str]][[exponent_str]][[ranking]] <- list()
        
        for (name in names(results_dict)){
          gsego_local     <- GOgse_results[[cutbypval_str]][[exponent_str]][[ranking]][[name]]
          GO_descriptions <- gsego_local$Description
          keyword_pathway_dict[[cutbypval_str]][[exponent_str]][[ranking]][[name]] <- find_matching_pathways(GO_descriptions, go_category_keywords_dict)
        }
      }
    }
  }
}

#--- Export GO terms in each category ---#
{
  for (cutbypval in cutbypval_list){
    cutbypval_str <- sprintf('cutbypval#%s',cutbypval)
    
    for (exponent in exponent_list){
      exponent_str <- sprintf('exponent#%d',exponent)
      
      for (ranking in ranking_list){
        
        for (name in names(results_dict)){
          gsego_local           <- GOgse_results[[cutbypval_str]][[exponent_str]][[ranking]][[name]]
          all_pathways          <- gsego_local$Description
          category_pathway_dict <- keyword_pathway_dict[[cutbypval_str]][[exponent_str]][[ranking]][[name]]
          
          #--- Output file ---#
          parameters_description <- sprintf("line#%s__%s__%s__%s", line_name, ranking, exponent_str, cutbypval_str)
          file_name              <- sprintf("GOtermsCategories__contrast#%s__%s.txt", name, parameters_description)
          file_path              <- sprintf("%s%s/%s%s", path_results, name, "GOterms/Categories/", file_name)
          dir.create(dirname(file_path), recursive=TRUE, showWarnings=FALSE)
          file_conn <- file(file_path, "w")
          gsego_local_copy           <- data.frame(gsego_local[,c("NES", "p.adjust", "Description")])
          NES_ranks                  <- length(gsego_local_copy$NES) + 1 - rank(abs(gsego_local_copy$NES))
          gsego_local_copy$rank      <- NES_ranks
          gsego_local_copy           <- gsego_local_copy[, c("rank", "NES", "p.adjust", "Description")]
          rownames(gsego_local_copy) <- gsego_local_copy$Description
          
          #--- Write categories ---#
          for (keyword in names(category_pathway_dict)) {
            pathways              <- category_pathway_dict[[keyword]]
            pathways_table        <- gsego_local_copy[pathways,]
            pathways_table_sorted <-pathways_table[order(abs(pathways_table$NES), decreasing=TRUE),]
            substring      <- paste(rep("-", times=50), collapse="")
            string         <- sprintf("%s %s %s", substring, keyword, substring)
            start_position <- (nchar(string) - 50) / 2 + 1
            string         <- sprintf("\n\n %s \n\n", substr(string, start_position, start_position + 49))
            cat(string, file=file_path, append=TRUE)
            pathways_table_sorted$NES      <- sprintf("%.2f", pathways_table_sorted$NES)
            pathways_table_sorted$p.adjust <- sprintf("%.2e", pathways_table_sorted$p.adjust)
            write.table(pathways_table_sorted, file=file_path, sep="\t\t", append=TRUE, quote=FALSE, row.names=FALSE)
          }
          close(file_conn)
        }
      }
    }
  }
}

message("GO terms categories exported.")


######################################
#      Plot grouped GO terms NES     #
######################################

message("Plotting grouped GO terms NES...")

{
  for (cutbypval in cutbypval_list){
    cutbypval_str <- sprintf('cutbypval#%s',cutbypval)
    
    for (exponent in exponent_list){
      exponent_str <- sprintf('exponent#%d',exponent)
      
      for (ranking in ranking_list){
        
        for (name in names(results_dict)){
          title               <- sprintf("%s - %s, exp:%d", line_name, name, exponent)
          keyword_pathway_set <- keyword_pathway_dict[[cutbypval_str]][[exponent_str]][[ranking]][[name]]
          gsego_local         <- GOgse_results[[cutbypval_str]][[exponent_str]][[ranking]][[name]]
          triples             <- find_xysize_NES_categories(keyword_pathway_set, gsego_local)
          x_range  <- c(-6.5, 6.5)
          x_breaks <- c(-6, -3, 0, 3, 6)
          if(ranking=='ranking#SlfcLogP'){
            x_range  <- c(-9, 9)
            x_breaks <- c(-8, -4, 0, 4, 8)
          }
          if(exponent==1){
            x_range  <- c(-3, 3)
            x_breaks <- c(-3, -1.5, 0, 1.5, 3)
          }
          categories_plot        <-  plot_GO_categories(triples, keyword_pathway_set, x_range=x_range, x_breaks=x_breaks, title=title)
          parameters_description <- sprintf("line#%s__%s__%s__%s", line_name, ranking, exponent_str, cutbypval_str)
          figure_name            <- sprintf("CategoriesGOterms__contrast#%s__%s.jpeg", name, parameters_description)
          file_path              <- sprintf("%s%s/GOterms_categories/%s", path_figures, name, figure_name)
          dir.create(dirname(file_path), recursive=TRUE, showWarnings=FALSE)
          ggsave(file_path, plot=categories_plot, device='jpeg', scale=1, width=14, height=13, units='cm', dpi=500, create.dir=TRUE)
        }
      }
    }
  } 
}

message("GO term category plots saved.")
message("Analysis completed")
