###############################
#           Colors            #
###############################

{
  color_lavander="#DBB2D1"
  color_lred='#FF616D'
  color_bittersweet="#F96B57"
  color_whitelavander="#FFEBF0"
  color_folly="#ff3366"
  color_palepurple="#FEEBFF"
  color_mauve="#FBC2FF"
  color_fuchsia="#F231FF"
  color_darkfuchsia="#AB00B8"
  color_cerise="#DF00A3"
  color_violet="#8D00B8"
  color_merblue="#44427B"
  color_zaffre="#0021A3"
  color_palblue="#1741FF"
  color_lblue='#3C8DAD'
  color_skyblue="#05C1FF"
  color_paleazure='#85e0ff'
  color_paleazurelight="#C2F0FF"
  color_turquoise="#55D6BE"
  color_nonphotoblue="#92d5e6"
  color_tiffany="#9CE8E0"
  color_mustard="#ffd35c"
  color_yellow="#ffc31f"
  color_orange="#f5a300"
  color_sand="#EEC7AA"
  color_linen="#F8E9DD"
  color_magenta="#FF00FF"
  color_gold="#FFD700"
}


###############################
#       Plot parameters       #
###############################

common_theme <- theme(
  plot.title=element_text(size=11, color="black"),  # Plot title font size
  axis.title.x=element_text(size=11, color="black"),  # X-axis title font size
  axis.title.y=element_text(size=11, color="black"),  # Y-axis title font size
  axis.text.x=element_text(size=9, color="black"),  # X-axis text font size
  axis.text.y=element_text(size=9, color="black"),  # Y-axis text font size
  legend.title=element_text(size=11, color="black"),  # Legend title font size
  legend.text=element_text(size=9, color="black"),  # Legend text font size
  axis.line=element_line(linewidth=0.3),  # Axis line width
  axis.ticks=element_line(linewidth=0.3)  # Ticks line width and color
  
)


###############################
#      GO term keywords       #
###############################

{
  go_category_keywords_dict <- list(RNA_processing_splicing=list(c("RNA", "processing"),
                                                              c("RNA", "metabol"),
                                                              c("splicing")),
                                    DNA_replication_repair=list(c("DNA", "replication"), 
                                                                c("DNA", "repair"), 
                                                                c("DNA", "process"), 
                                                                c("helicase")), 
                                    Transcription=list(c("transcription")),
                                    Translation=list(c("translation")),
                                    Cell_cycle=list(c("cell cycle"), 
                                                    c("mitosis"), 
                                                    c("meiosis")),
                                    Nucleus_chromosome=list(c("nucleo"), 
                                                            c("nuclea"),
                                                            c("nucleu"),
                                                            c("chromosom")),
                                    Ribosome=list(c("ribo")),
                                    Development=list(c("development")), 
                                    Membrane_transport=list(c("membrane"), 
                                                            c("transport")),
                                    Mitochondria=list(c("membrane")),
                                    Cytoskeleton=list(c("cytockelet"),
                                                      c("microtubul")),
                                    Extracellular_matrix=list(c("extracel", "matrix")))
  
  go_category_names <- c("RNA processing/splicing", 
                         "DNA replication/repair", 
                         "Transcription", 
                         "Translation",
                         "Cell cycle", 
                         "Nucleus/chromosome", 
                         "Ribosome",
                         "Development", 
                         "Membrane/transport", 
                         "Mitochondria", 
                         "Cytoskeleton", 
                         "Extracellular matrix")
}


###############################
#      Metadata columns       #
###############################

metadata_columns <- c("baseline", "treatment")


###############################
#          Functions          #
###############################

{
  #----- find_matching_pathways -----#
  find_matching_pathways <- function(GO_descriptions, go_category_keywords_dict) {
    keyword_pathway_dict <- list()
    unique_pathways_set <- list()
    # Iterate through each pathway_category in go_category_keywords_dict
    for (pathway_category in names(go_category_keywords_dict)) {
      string_lists <- go_category_keywords_dict[[pathway_category]]
      matching_pathways <- list()
      # Iterate through each string_list for the current pathway_category
      for (string_list in string_lists) {
        current_matching_pathways <- list()
        # Iterate through each pathway in GO_descriptions
        for (pathway in GO_descriptions) {
          # Check if all strings in the current string_list are present in the pathway
          if (all(sapply(string_list, function(string) grepl(string, pathway)))) {
            current_matching_pathways <- c(current_matching_pathways, list(pathway))
          }
        }
        # Add the current matching pathways to the overall list
        matching_pathways <- c(matching_pathways, current_matching_pathways)
      }
      # Combine the lists of matching pathways into a single list
      combined_pathways <- if(length(matching_pathways) > 0) do.call(c, matching_pathways) else list()
      # Remove duplicates from the combined list
      unique_pathways <- unique(combined_pathways)
      # Remove pathways that are already in the set
      unique_pathways <- unique_pathways[!unique_pathways %in% unique_pathways_set]
      unique_pathways <- if(length(unique_pathways) > 0) unique_pathways else character(0) 
      # Update the set of unique pathways
      unique_pathways_set <- union(unique_pathways_set, unique_pathways)
      # Store the unique list of pathways for the current pathway_category in the new dictionary
      keyword_pathway_dict[[pathway_category]] <- unique_pathways
    }
    return(keyword_pathway_dict)
  }
  
  
  #----- find_xysize_NES_categories -----#
  find_xysize_NES_categories <- function(keyword_pathway_set, gsego_results) {
    points <- list()
    y_val <- 1
    y <- numeric(0)
    x <- numeric(0)
    size <- numeric(0)
    for (category in names(keyword_pathway_set)){
      for (pathway in keyword_pathway_set[[category]]) {
        index <- which(pathway==gsego_results$Description)
        x <- c(x, gsego_results$NES[[index]])
        size <- c(size, -log10(gsego_results$p.adjust[[index]]))
        y <- c(y, y_val)
      }
      y_val <- y_val + 1
    }
    # Your list of triples (x, y, size)
    triples <- data.frame(x=x,
                          y=y,
                          size=size)
    return(triples)
  }
  
  
  #----- plot_GO_categories(triples) -----#
  plot_GO_categories <- function(triples, keyword_pathway_set, x_range=c(-6.5, 6.5), 
                                 x_breaks=c(-6, -3, 0, 3, 6), title='NES GO by categories'){
    y_range <- c(1, length(keyword_pathway_set))
    y_breaks <- seq_along(names(keyword_pathway_set))
    y_labels <- go_category_names
    categories_plot <- ggplot(triples, aes(x=x, y=y, size=size)) +
      geom_point(shape=1, color=color_folly) +  # Use shape=1 for circular markers
      scale_size_continuous(range=c(1, 8), name='-log(padj)') +  # Adjust the range of marker sizes
      scale_y_continuous(breaks=y_breaks, labels=y_labels) +  # Specify breaks and labels for the y-axis
      scale_x_continuous(breaks=x_breaks) + 
      labs(x="Enrichment (NES)", y=NULL) +
      theme_classic() +  # Use a minimal theme
      theme(axis.text.x=element_text(angle=0, hjust=1, size=11, color="black"),  # Rotate x-axis labels
            panel.grid.major.y=element_line(color="gray80"),  # Add horizontal grid lines
            axis.line.x=element_line(color="black", linewidth=0.5),
            axis.ticks.y=element_blank(),
            axis.line.y=element_blank(), 
            axis.text.y=element_text(size=11, color="black"), 
            axis.title.x=element_text(size=11, color="black")) +
      geom_vline(xintercept=0, color="black", size=0.5) +
      coord_cartesian(xlim=x_range, ylim=y_range) +
      ggtitle(title)
    print(categories_plot)
    return(categories_plot)
  }
  
  
  #----- plot_lfc_vs_length -----#
  plot_lfc_vs_length <- function(ensembl_ids, gene_len, comparison, results_dict, 
                                 line_name, gene_type="", y_range=c(-8, 8), 
                                 y_breaks=c(-8, -4, 0, 4, 8), method="spearman") {
    comparison_results <- results_dict[[comparison]]
    
    # Ensure gene names (rownames) match ensembl_ids
    selected_indices <- match(ensembl_ids, rownames(comparison_results))
    selected_indices <- selected_indices[!is.na(selected_indices)]  # Remove NA matches
    
    # Subset only relevant genes
    lfc <- comparison_results$log2FoldChange[selected_indices]
    padj <- comparison_results$padj[selected_indices]

    # Create data frame and categorize genes
    lfc_length_data <- data.frame(length=gene_len, lfc=lfc, padj=padj)
    lfc_length_data <- na.omit(lfc_length_data)  # Remove rows with NA values

    # Categorize genes based on the conditions
    lfc_length_data$category <- "no diff."
    lfc_length_data$category[(lfc_length_data$padj < 1e-3) & (lfc_length_data$lfc > 0.75)] <- "up"
    lfc_length_data$category[(lfc_length_data$padj < 1e-3) & (lfc_length_data$lfc < -0.75)] <- "down"
    lfc_length_data$category <- factor(lfc_length_data$category, levels=c("no diff.", "down", "up"))
    lfc_length_data <- lfc_length_data[order(lfc_length_data$category), ]
    
    # Compute correlation
    lfc_length_data$log10_length <- log10(lfc_length_data$length)
    cor_data <- lfc_length_data#[lfc_length_data$category %in% c("up", "down"), ]
    if (nrow(cor_data) >= 100) {
      correlation <- cor(cor_data$lfc, cor_data$log10_length, method=method)
      cor_label <- sprintf("r=%.2f", correlation)
    } else {
      cor_label <- "r=NA"
    }
    
    # Define custom colors
    custom_colors <- c("up"=color_fuchsia, "down"=color_palblue, "no diff."="gray90")
    
    # Create the scatter plot
    lfc_length_scatterplot <- ggplot(lfc_length_data, aes(x=length, y=lfc, fill=category)) +
      geom_point(shape=21, size=0.8, alpha=1., stroke=0.) +  
      scale_fill_manual(values=custom_colors, breaks=c("up", "down")) + 
      theme_classic() +
      common_theme +
      scale_x_log10(labels=scales::trans_format("log10", scales::math_format(10^.x))) +
      coord_cartesian(ylim=y_range) +
      scale_y_continuous(breaks=y_breaks) +
      annotate("text", x=Inf, y=Inf, label=cor_label, hjust=1.1, vjust=1.5, size=3.5) +
      labs(x="length (bp)", y="lfc", title=gene_type)
    
    # Print and return the plot
    print(lfc_length_scatterplot)
    return(lfc_length_scatterplot)
  }
  
  
  #----- rename_duplicates -----#
  rename_duplicates <- function(strings) {
    # Create a named vector to store counts
    counts <- table(strings)
    # Loop through the strings
    for (i in 1:length(strings)) {
      # Get the current string
      current_string <- strings[i]
      # Check if the current string has duplicates
      if (sum(strings==current_string)>1) {
        # If it has duplicates, rename it with appended index
        counter <- 1
        for (j in i:length(strings)) {
          if (strings[j]==current_string){
            strings[j] <- sprintf("%s#%d", strings[j], counter)
            counter <- counter + 1
          }
        }
      }
    }
    return(strings)
  }
  
  
  #----- get_colors -----#
  ## Function to create vector of colors ##
  get_colors <- function(levels, colors_dict) {
    # Match levels to colors in dictionary using the mapping function
    colors <- sapply(levels, function(level) {
      color_key <- gsub("\\.", "_", level)
      colors_dict[[color_key]]
    })
    return(colors)
  }
}