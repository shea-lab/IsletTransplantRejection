# Gene conversion functions ---- 
# These functions are called out in the first central function of the package
convMensID.Msym <- function(counts, mouse_mart = NULL) { # Convert mouse IDs to symbols
  genes <- rownames(counts) # counts df must have genes as rownames, not the first column
  require("biomaRt") # Necessary package
  if (is.null(mouse_mart)) { # Default mart
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE,
                     host="https://dec2021.archive.ensembl.org/")
  } else {
    mouse <- mouse_mart # Mart specified by user
  }
  new_genes <- getLDS(attributesL = c("mgi_symbol"), filters = "ensembl_gene_id", 
                      values = genes, mart = mouse, attributes = c("ensembl_gene_id"), 
                      martL = mouse, uniqueRows=TRUE)
  colnames(new_genes) <- c("old", "new")
  new_genes <- new_genes[!duplicated(new_genes[1]),]
  new_genes <- new_genes[!duplicated(new_genes[2]),]
  rownames(new_genes) <- new_genes$old
  counts_conv <- merge(new_genes, counts, by = "row.names", all.x=T)
  rownames(counts_conv) <- counts_conv$new
  counts_conv <- counts_conv[,-1:-3]
  return(counts_conv)
}
convHensID.Hsym <- function(counts, human_mart = NULL) { # Convert human IDs to symbols
  genes <- rownames(counts) # counts df must have genes as rownames, not the first column
  require("biomaRt") # Necessary package
  if (is.null(human_mart)) { # Default mart
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE,
                     host="https://dec2021.archive.ensembl.org/")
  } else {
    human <- human_mart # Mart specified by user
  }
  new_genes <- getLDS(attributesL = c("hgnc_symbol"), filters = "ensembl_gene_id", 
                      values = genes, mart = human, attributes = c("ensembl_gene_id"), 
                      martL = human, uniqueRows=TRUE)
  colnames(new_genes) <- c("old", "new")
  new_genes <- new_genes[!duplicated(new_genes[1]),]
  new_genes <- new_genes[!duplicated(new_genes[2]),]
  rownames(new_genes) <- new_genes$old
  counts_conv <- merge(new_genes, counts, by = "row.names", all.x=T)
  rownames(counts_conv) <- counts_conv$new
  counts_conv <- counts_conv[,-1:-3]
  return(counts_conv)
}

# Central functions ----

# In the first function, counts are filtered as specified.
flexiDEG.function1 <- function(counts, # Counts df, genes=rows & samples=cols, 
                               # gene names can be row names or the first col
                               metadata, # Metadata with first cell in the col reading "Samples"
                               # also needs a col with first cell reading "Group"
                               # and a col with first cell reading "Batch" if going to batch correct
                               convert_genes = FALSE, # Convert gene IDs to symbols 
                               exclude_riken = TRUE, # Remove riken genes prior to analyses
                               exclude_pseudo = TRUE, # Remove psuedo genes prior to analyses 
                               batches = FALSE, # Analyze multiple data sets
                               quality = FALSE, # Report on data quality
                               variance = TRUE, # Filter out genes that change little
                               use_pseudobulk = FALSE # New argument to specify pseudobulk usage
) {
  # Check inputs ---- 
  sample_names <- metadata$Samples # Extract sample names 
  sample_num <- length(sample_names) # Extract number of samples 
  groups <- metadata$Group # Extract sample groups 
  message <- "Please indicate if the genes are rownames (1) or the first column (2) of the counts df: "
  user_input <- readline(prompt = message) # Check for consistency between metadata & counts
  if (user_input == 1) {
    if (ncol(counts) != nrow(metadata)) {stop("Number of samples in counts & metadata are inconsistent")}
    if (!all(sample_names %in% colnames(counts))) {stop("Sample names in counts & metadata are inconsistent")} 
    genes <- make.names(rownames(counts), unique = T) # Cant start rows w/ numbers, add a letter to genes starting w/ numbers
  } else if (user_input == 2) {
    if (ncol(counts) - 1 != nrow(metadata)){stop("Number of samples in counts & metadata are inconsistent")}
    if (!all(sample_names %in% colnames(counts[, -c(1)]))) {stop("Sample names in counts & metadata are inconsistent")} 
    genes <- make.names(counts[, 1], unique = T) # Cant start rows w/ numbers, add a letter to genes starting w/ numbers
    counts <- counts[, -1] # Remove first column of genes to add it back as row names instead 
    colnames(counts) <- c(sample_names) # Add sample names as column names
    rownames(counts) <- genes # Add gene names as row names
  } else {stop("Genes must be listed as rownames or first column in the counts")
  }
  if (sum(duplicated(genes)) != 0) { # Check for duplicate genes
    genes2 <- genes[!duplicated(genes)] # Genes which are not duplicated
    counts <- counts[genes2, ] # Remove duplicated genes
    removed <- length(genes) - length(genes2) # Number of removed genes
    print(paste(removed, " duplicate genes removed")) # Display number 
    genes <- rownames(counts) # Redefine genes without duplicates
  } 
  
  # Convert genes ---- 
  if (convert_genes == TRUE) { # Only converts ensembl IDs to symbols for mouse or human
    if (substring(genes[1], 1, 3) != "ENS") { # Check for genes as ensembl IDs
      stop("Only genes as ensembl IDs can be converted by this function.")
    } else {
      if (substring(genes[1], 4, 7) == "MUSG") { # If they are mouse
        counts <- convMensID.Msym(counts) # Function to convert
        dt <- "mouse"
      } else if (substring(genes[1], 4, 4) == "G") { # If they are human
        counts <- convHensID.Hsym(counts) # Function to convert
        dt <- "human"
      } else { 
        stop("Genes are not mouse or human; convert manually and then input into function.")
      }}
    genes <- rownames(counts) # Redefine genes once converted
  } else { # Check if genes are human or mouse, even though we won't be converting
    letters_before_numbers <- gsub("[^A-Za-z]+.*$", "", genes[1:10]) # Find letters before the first number
    all_caps <- all(grepl("^[A-Z]+$", letters_before_numbers)) # Check if all characters are uppercase
    if (all_caps == TRUE) {
      dt <- "human"
    } else {
      dt <- "mouse"
    }}
  
  # Exclude genes ---- 
  if (exclude_riken == TRUE | exclude_pseudo == TRUE) { 
    if (substring(genes[1], 1, 3) == "ENS") { # Check for genes as ensembl IDs
      stop("Riken genes cannot be excluded when given as ensembl IDs")
    }}
  if (exclude_riken == TRUE) { 
    if (dt == "mouse") {
      genes2 <- genes[grep("Rik", genes, invert = T)] # Find Riken genes
      removed <- length(genes) - length(genes2) # Number of removed genes
      print(paste(removed, "Riken genes removed")) # Display number 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Riken genes
    } else if (dt == "human") {
      genes2 <- genes[grep("Rik", genes, invert = T)] # Find Riken genes
      removed <- length(genes) - length(genes2) # Number of removed genes
      print(paste(removed, "RIKEN genes removed")) # Display number 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Riken genes
    }}
  if (exclude_pseudo == TRUE) { 
    if (dt == "mouse") {
      genes2 <- genes[grep("Ps", genes, invert = T)] # Find Pseudo genes
      genes2 <- genes2[grep("Gm", genes2, invert = T)] # Find Pseudo genes
      removed <- length(genes) - length(genes2)  # Number of removed genes
      print(paste(removed, "Pseudo genes (starting with Ps or Gm) removed")) 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Pseudo genes
    } else if (dt == "human") {
      genes2 <- genes[grep("PS", genes, invert = T)] # Find Pseudo genes
      genes2 <- genes2[grep("GM", genes2, invert = T)] # Find Pseudo genes
      removed <- length(genes) - length(genes2)  # Number of removed genes
      print(paste(removed, "Pseudo genes (starting with PS or GM) removed")) 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Pseudo genes
    }}
  
  if (use_pseudobulk== FALSE){
    # Remove unexpressed genes ---- 
    counts[is.na(counts)] <- 0 # Set any missing values equal to zero
  }
  else{
    # Check fraction of missing values in counts ----
    missing_fraction <- mean(is.na(counts)) # Calculate fraction of missing values
    print(paste("Fraction of missing values in counts data:", round(missing_fraction, 4)))
    
    # Function to filter genes based on a missing value threshold
    filter_genes <- function(df, max_missing_threshold = 0.2) {
      keep_genes <- apply(df, 1, function(x) mean(is.na(x)) <= max_missing_threshold)
      df_filtered <- df[keep_genes, ]
      return(df_filtered)
    }
    message <- "Please maximum fraction of missing values allowed: "
    max_missing_threshold <- as.numeric(readline(prompt = message))
    counts  <- filter_genes(counts , max_missing_threshold = 0.2)
    # Impute missing values with median across samples
    median_impute <- function(df) {
      df[] <- lapply(df, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
      return(df)
    }
    counts <- median_impute(counts)
  }
  counts2 <- counts[rowSums(counts) > 0,] # Remove genes unexpressed across all samples
  removed <- nrow(counts) - nrow(counts2) # Number of removed genes
  print(paste(removed, "unexpressed genes removed")) # Display number
  counts <- counts2 # Update counts df
  counts2 <- counts[rowSums(counts == 0) <= sample_num*(0.85),] # Unexpressed in 85%+ of samples
  removed <- nrow(counts) - nrow(counts2) # Number of removed genes
  print(paste(removed, "mostly zero genes removed")) # Display number
  counts <- counts2 # Update counts df
  
  # Batch correction ---- 
  if (batches == TRUE) { # Batch correction requires three things:
    # 1) counts df w/ RAW data from all batches combined, rows=genes cols=samples
    # 2) vector defining sample batches
    # 3) vector defining sample groups
    Batch <- metadata$Batch 
    coldata <- data.frame(cbind(groups, Batch)) # Combine groups and batch 
    row.names(coldata) <- sample_names # DESeq2 uses raw counts; rows:genes & cols:samples
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Batch+groups)
    paste(nrow(dds), " genes input for batch correction", sep="")
    suppressMessages(dds <- DESeq(dds)) # Analyze previously specified variables
    vsd <- vst(dds, blind = FALSE) # Estimates dispersion trend & stabilizes transformation
    mat <- assay(vsd)
    mat <- limma::removeBatchEffect(mat, vsd$Batch)
    assay(vsd) <- mat
    counts <- as.data.frame(assay(vsd)) # Summarized assay value, Update counts df
  }
  
  #ps_counts <- log(counts + 1, base = 2) # Pseudo-normalize
  ps_counts<- counts
  # All filters are per sample per gene
  message <- "Please enter expression cutoff for low filter (required): "
  filter_low <- as.numeric(readline(prompt = message))
  if (is.na(filter_low)) {stop("Invalid input. Please enter a numeric value.")} 
  keep1 <- rowSums(ps_counts) > (filter_low*sample_num) # Low count filter
  ps_filt <- ps_counts[keep1,] # Apply filter to pseudonorm counts
  removed <- nrow(ps_counts) - nrow(ps_filt) # Number of filtered genes
  print(paste(removed, "genes removed due to expression below the low filter"))
  message <- "Please enter expression cutoff for group filter (0 for no group filter): "
  filter_group <- as.numeric(readline(prompt = message))
  if (is.na(filter_group)) {stop("Invalid input. Please enter a numeric value.")
  } else if (filter_group != 0) { 
    keep2 <- rowSums(ps_filt) >= (filter_group*sample_num) # Group filter, optional
    ps_filt2 <- ps_filt[keep2,] # Apply filter to pseudonorm counts
    removed <- nrow(ps_filt) - nrow(ps_filt2) # Number of filtered genes
    print(paste(removed, "genes removed due to expression below the low group filter"))
    ps_filt <- ps_filt2 # Update counts df
  }
  message <- "Please enter expression cutoff for high filter (0 for no high filter): " 
  filter_high <- as.numeric(readline(prompt = message))
  if (is.na(filter_high)) {stop("Invalid input. Please enter a numeric value.")
  } else if (filter_high != 0) { 
    keep3 <- rowSums(ps_filt) < (filter_high*sample_num) # High count filter, optional
    ps_filt2 <- ps_filt[keep3,] # Apply filter to pseudonorm counts
    removed <- nrow(ps_filt) - nrow(ps_filt2) # Number of filtered genes
    print(paste(removed, "genes removed due to expression above the high filter"))
    ps_filt <- ps_filt2 # Update counts df
  }
  genes_filt <- rownames(ps_filt) # Update genes vector for downstream 
  counts_filt <- counts[genes_filt,] # Update counts df for downstream 
  
  # Sample quality ---- 
  if (quality == TRUE) { 
    suppressMessages(ps_LF <- melt(ps_counts, variable.name = "Samples", value.name = "Count")) # Long form
    ps_LF <- merge(ps_LF, metadata, by = "Samples") # Merge w/ metadata to get the group for each sample
    ps_LF <- ps_LF[!is.na(ps_LF$Count), ] # Remove any missing values
    suppressMessages(psf_LF <- melt(ps_filt, variable.name = "Samples", value.name = "Count")) # Long form
    psf_LF <- merge(psf_LF, metadata, by = "Samples") # Merge w/ metadata to get the group for each sample
    psf_LF <- psf_LF[!is.na(psf_LF$Count), ] # Remove any missing values
    if (sample_num <= 20) { # If 20 samples or less, print to 1 page
      par(mfrow=c(4,ceiling(sample_num/4))) # All histograms on a page w/ 4 rows
    } else {par(mfrow=c(3,5))} # Histograms on multiple pages w/ 3 rows & 5 cols
    par(mar = c(4,1,2,1), oma = c(0.5,1,0,0)) # Adjust margins to better use space
    sapply(1:sample_num, function(x) hist(counts[,x], ylim=c(0,1000), breaks=100,
                                          xlab=sample_names[x], ylab = "", main = ""))
    mtext("Counts, 100 bins", side = 3, line = -2, outer = T) # Add title
    sapply(1:sample_num, function(x) hist(ps_counts[,x], ylim=c(0,1000), breaks=100,
                                          xlab=sample_names[x], ylab = "", main = ""))
    mtext("Log Normalized Counts, 100 bins", side = 3, line = -2, outer = T) # Add title
    sapply(1:sample_num, function(x) hist(ps_filt[,x], ylim=c(0,1000), breaks=100,
                                          xlab=sample_names[x], ylab = "", main = ""))
    mtext("Log Normalized Counts Filt, 100 bins", side = 3, line = -2, outer = T) # Add title
    # Box & Density Plots, how samples compare in mean, variance & any outliers
    p1 <- ggplot(ps_LF, aes(x = Samples, y = Count)) + 
      geom_boxplot() + xlab("") +  ylab("Log Normalized Counts") +
      scale_x_discrete(breaks = sample_names, labels = sample_names) + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
      ggtitle("Counts Distribution")
    grid.draw(p1) # Plot distribution
    p2 <- ggplot(ps_LF, aes(x = Count, fill = Group, group = Group)) + 
      geom_density(alpha = 0.4, color = "black") + 
      theme(legend.position = "top") + 
      xlab("Log Normalized Counts") + 
      ggtitle("Combined Density Distribution of Counts")
    grid.draw(p2) # Plot distribution
    # Box & Density Plots, how samples compare in mean, variance & any outliers
    p1 <- ggplot(psf_LF, aes(x = Samples, y = Count)) + 
      geom_boxplot() + xlab("") +  ylab("Log Normalized Counts") +
      scale_x_discrete(breaks = sample_names, labels = sample_names) + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
      ggtitle("Counts Distribution")
    grid.draw(p1) # Plot distribution
    p2 <- ggplot(psf_LF, aes(x = Count, fill = Group, group = Group)) + 
      geom_density(alpha = 0.4, color = "black") + 
      theme(legend.position = "top") + 
      xlab("Log Normalized Counts") + 
      ggtitle("Combined Density Distribution of Counts")
    grid.draw(p2) # Plot distribution
  }
  
  # Variance filtering ---- 
  if (variance == TRUE) { 
    # Filter by variance across all samples to exclude genes that don't change much
    mn_cut <- 0.2 # Mean cutoff
    var_cut <- 0.01 # Variance cutoff
    df <- cbind.data.frame("x" = rowMeans(counts_filt), # Calculate Mean & Var/Mean
                           "y" = rowVars(counts_filt) / rowMeans(counts_filt)) 
    g_var_all <- rownames(subset(df, x > mn_cut & y > var_cut)) # Genes to keep
    # Plot variance and mean
    par(mfrow = c(1, 1)) # Set the layout to a single plot in a 1x1 grid
    with(df, plot(x, y, main = paste("All Var vs Mean, ", length(g_var_all), 
                                     " genes w/ Mean > ", mn_cut, " & Var > ", var_cut, sep = ""), 
                  pch = 20, cex = 1, xlab = "Mean across all samples", 
                  ylab = "Var/Mean across All samples"))
    with(subset(df, x > mn_cut & y > var_cut), 
         points(x, y, pch = 20, col = "red3", cex = 1)) # Genes above cutoff are red
    abline(v = mn_cut, col = "red3", lty = 2, lwd = 1.5) # Line for mean cutoff
    abline(h = var_cut, col = "red3", lty = 2, lwd = 1.5) # Line for variance cutoff
    removed <- nrow(counts_filt) - length(g_var_all) # Number of filtered genes
    print(paste(removed, "genes removed due to expression below the variance filter"))
    counts_filt <- counts_filt[g_var_all,]  # Update counts_filt 
    #   if (length(unique(groups)) == 2) { # Additional variance by group 
    #     num_groups <- length(unique(groups)) # Determine the number of groups 
    #     # ++++ This needs to be finished still
    #     
    #     removed <- nrow(counts_filt) - length(g_var_group) # Number of filtered genes 
    #     print(paste(removed, " genes removed by the pairwise group variance filter")) 
    #     counts_filt <- counts_filt[g_var_group,]  # Update counts_filt 
    # }
  }
  return(counts_filt) # Final function 1 output
}

# In the second function, we perform ANOVA by sample groups.
flexiDEG.function2 <- function(counts_filt, # df w/ genes as rownames, Function 1 output
                               metadata, sig_cutoff = NULL) {
  if (ncol(counts_filt) > nrow(counts_filt)) {counts_filt <- t(counts_filt)} # Transpose if cols=genes 
  groups <- metadata$Group # Extract sample groups 
  groups <- factor(groups, levels = unique(groups), labels = unique(groups)) # Groups as factors
  gene_pvals <- apply(counts_filt, 1, function(gene_row) { # Perform ANOVA for each gene
    aov_result <- aov(gene_row ~ groups, data = counts_filt)
    return(summary(aov_result)[[1]]$`Pr(>F)`[1])
  })
  gene_pvals <- data.frame(gene = row.names(counts_filt), p_val = gene_pvals) # Genes & p-values
  if (is.null(sig_cutoff)) {sig_cutoff <- 0.05} # If empty, set p value cutoff to 0.05
  gene_sig <- gene_pvals[gene_pvals$p_val < sig_cutoff, ] # Filter results by p value
  counts_anova <- counts_filt[gene_sig$gene,]
  # Print the message indicating the number of groups and significant genes/pathways after ANOVA
  if (sum(grepl("[_ ]", rownames(counts_filt))) >= 5) {
    # If there are 2 or more rownames in counts_anova with "_" or space, consider them as pathways
    print(paste("ANOVA performed for", length(unique(groups)), "groups with", nrow(counts_anova), "significant pathways"))
  } else {
    # Otherwise, consider them as genes
    print(paste("ANOVA performed for", length(unique(groups)), "groups with", nrow(counts_anova), "significant genes"))
  }
  return(counts_anova) # Final function 2 output
}

# In the third function, we perform t-tests between pairwise sample groups.
flexiDEG.function3 <- function(counts_filt, metadata, heatmap = FALSE,
                               fdr_cutoff = NULL, logfc_cutoff = NULL) { 
  groups <- metadata$Group # Extract sample groups 
  if ((length(unique(groups)) != 2) & !("Simple" %in% colnames(metadata))) { 
    stop("Samples must be in 2 groups for comparison") # Needs pairwise comparison 
  } else if ((length(unique(groups)) != 2) & ("Simple" %in% colnames(metadata))) {
    simple <- metadata$Simple 
    if (length(unique(simple)) != 2) {
      stop("Samples must be in 2 groups for comparison") # Needs pairwise comparison 
    } 
    extracted_info <- metadata[metadata$Simple %in% simple, ]
    sample_ids <- extracted_info$Samples # Extract SampleID values
    extracted_data <- counts_filt[, colnames(counts_filt) %in% sample_ids] # colnames = sample IDs
    group <- metadata$Simple[metadata$Samples %in% sample_ids] # Prepare data for analysis
    group_data <- extracted_data[, group %in% simple]
    t_group0 <- group_data[, group %in% unique(simple)[1]] # Extract two groups
    t_group1 <- group_data[, group %in% unique(simple)[2]]
  } else if (length(unique(groups)) == 2) {
    extracted_info <- metadata[metadata$Group %in% groups, ]
    sample_ids <- extracted_info$Samples # Extract SampleID values
    extracted_data <- counts_filt[, colnames(counts_filt) %in% sample_ids] # colnames = sample IDs
    group <- metadata$Group[metadata$Samples %in% sample_ids] # Prepare data for analysis
    group_data <- extracted_data[, group %in% groups]
    t_group0 <- group_data[, group %in% unique(groups)[1]] # Extract two groups
    t_group1 <- group_data[, group %in% unique(groups)[2]]
  }
  log_fc <- rowMeans(t_group1) - rowMeans(t_group0) # Calculate log-fold change
  t_test <- t.test(t_group1, t_group0) # Perform t-test & extract p-values
  pvalues <- t_test$p.value
  fdr <- p.adjust(pvalues, method = "fdr") # Multiple testing via Benjamini-Hochberg
  results <- data.frame(gene = row.names(group_data), log_fc, fdr) # df with genes, log-fold change, & FDR
  if (is.null(fdr_cutoff)) {fdr_cutoff <- 0.05} # Set cutoff if the user hasn't
  if (is.null(logfc_cutoff)) {logfc_cutoff <- 1.5} # Set cutoff if the user hasn't
  filt_results <- results[results$fdr < fdr_cutoff & results$log_fc > abs(logfc_cutoff),] # Filter by FDR & FC
  filt_data <- extracted_data[row.names(extracted_data) %in% filt_results$gene, ] # Filtered expression
  if (heatmap == TRUE) { # Create heatmap 
    groups <- metadata$Group[metadata$Samples %in% colnames(filt_data)] # Group info
    num_groups <- length(unique(groups)) # Number of unique groups
    group_colors <- setNames(brewer.pal(num_groups, "Set1"), unique(group)) # Named palette for groups
    colSide <- group_colors[groups] # Matrix of colors for the samples
    heatmap.2(as.matrix(filt_data), scale = "row", key = TRUE, key.ylab = NA, keysize = 1.0, 
              col = colorRampPalette(c("navy", "white", "firebrick3"))(100), trace = "none",
              xlab = "Samples", ylab = "Genes", margins = c(13, 13), key.title = NA, 
              ColSideColors = colSide, dendrogram = "column")
  }
  # Print the message indicating the number of significant genes/pathways after pairwise t-tests
  if (sum(grepl("[_ ]", rownames(counts_filt))) >= 5) {
    # If there are 5 or more rownames in counts_filt with "_" or space, consider them as pathways
    print(paste("Pairwise t-tests resulted in", nrow(filt_data), "significant pathways"))
  } else {
    # Otherwise, consider them as genes
    print(paste("Pairwise t-tests resulted in", nrow(filt_data), "significant genes"))
  }
  return(filt_data) # Final function 3 output
}

# In the fourth function, we run an Elastic Net regression with validation selected by group 
# size. There are three versions of the Elastic Net validation, and function 4 calls one.
flexiDEG.function4 <- function(counts_filt, metadata, validation_option = NULL) { 
  opt <- validation_option # Just to shorten this object name
  if (is.null(opt)) {opt <- 0} # Set to zero if not user specified
  groups <- as.integer(factor(metadata$Group, labels = 1:length(unique(metadata$Group))))
  groups <- as.data.frame(table(groups))
  min_group <- min(groups$Freq) # Number of samples in smallest group
  if (min_group < 3) { # CV method chosen based on group sizes & user input
    stop("Samples per group must be at least 3 for Elastic Net regression.")
  } else {
    if ((min_group == 3 && (opt != 0 || opt != 1))) {
      print("Samples per group must be >3 for Elastic Net regression with partitioning. 
              Performing simple Elastic Net regression instead.")
      genes <- flexiDEG.ENsimple(counts_filt, metadata)
      return(genes)
    } 
    if ((min_group >= 3 && min_group < 5 && opt !=2 && opt != 3 ) || opt == 1) {
      print("Performing simple Elastic Net regression.")
      genes <- flexiDEG.ENsimple(counts_filt, metadata)
      return(genes)
    } 
    if ((min_group >= 5 && min_group < 10 && opt != 1 && opt != 3) || opt == 2) {
      print("Performing Elastic Net regression with partitioning.")
      genes <- flexiDEG.ENpartition(counts_filt, metadata)
      return(genes)
    } 
    if ((min_group >= 10 && opt != 1 && opt!= 2 )|| opt == 3) {
      print("Performing Elastic Net regression with partitioning and bootstrapping.")
      genes <- flexiDEG.ENbootstrap(counts_filt, metadata)
      return(genes)
    }}}

flexiDEG.ENsimple <- function(counts_filt, metadata) { # No cross-validation, n=3 sample size
  start_t <- Sys.time() # Record when function begins processing data
  num_groups <- length(unique(metadata$Group)) # Groups for Elastic Net 
  groups <- as.integer(factor(metadata$Group, labels = 0:(num_groups - 1))) # Groups as factors 
  levels(groups) <- 0:(num_groups - 1) 
  groups <- factor(groups) 
  dataset <- as.data.frame(cbind(as.data.frame(groups), t(counts_filt))) # Add factor as first col
  dataset <- dataset[complete.cases(dataset), ] # For cases w/out missed items
  set.seed(123)
  train_samples <- dataset$groups %>% createDataPartition(p = 1, list = F) # Training set not split
  train_data <- dataset[train_samples, ]
  x <- model.matrix(groups~., train_data)[ ,-1] # Predictor variables
  y <- train_data$groups # Outcome variable
  foldid <- sample(1:length(y), size = length(y), replace = TRUE)
  num_alphas <- 21 # 21 alphas
  # Fit EN & calculate AUC for each alpha
  if (num_groups < 2) {
    stop("Samples in one group; Elastic Net regression cannot be performed")
  } else if (num_groups == 2) {
    print("Two sample groups; binomial Elastic Net regression will be performed")
    fit_evaluate_EN <- function(alphas, x, y, foldid, num_groups) { # Fit EN & calculate AUC for each alpha
      # Perform cross-validated Elastic Net regression for binary classification
      cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alphas, foldid = foldid, family = "binomial"))
      # Determine the lambda value that minimizes cross-validation error
      best.lambda <- cvfit$lambda.min
      # Fit the final Elastic Net model using the best lambda value obtained from cross-validation
      fit <- suppressMessages(glmnet(x, y, alpha = alphas, lambda = best.lambda, family = "binomial"))
      my_lambda <- cvfit$lambda.min
      Coef <- coef(fit, s = my_lambda)
      Index <- c()
      # Extracting non-zero indices
      non_zero_indices <- which(Coef != 0)
      # Exclude the intercept (which is usually the first coefficient)
      Index <- non_zero_indices[non_zero_indices != 1]
      return(list(ENGenes = Index)) # Return EN Genes for this alpha
    } 
  } else {
    print("More than two sample groups; multinomial Elastic Net regression will be performed")
    fit_evaluate_EN <- function(alphas, x, y, foldid, num_groups) { # Fit EN & calculate AUC for each alpha
      cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alphas, foldid = foldid, family = "multinomial"))
      my_lambda <- cvfit$lambda.min
      fit <- suppressMessages(glmnet(x, y, alpha = alphas, lambda = my_lambda, family = "multinomial"))
      Coef <- coef(fit, s = my_lambda)
      Index <- c()
      for (group in 1:num_groups) {
        Index <- c(Index, Coef[[as.character(group)]]@i[-1])
      }
      return(list(ENGenes = Index)) # Return EN Genes for this alpha
    }
  }
  ENgenes <- vector("list", length = num_alphas) # List for EN genes of each alpha
  for (i in 1:num_alphas) {
    pb <- txtProgressBar(min = 0, max = num_alphas, initial = i, style = 3)
    setTxtProgressBar(pb, i)
    alphas <- (i - 1) * 0.05 
    ENgenes[[i]] <- fit_evaluate_EN(alphas, x, y, foldid, num_groups)
  }
  ENgenes <- rev(unlist(ENgenes, recursive = FALSE))
  names(ENgenes) <- seq(1, 0, -0.05)
  total_t <- difftime(Sys.time(), start_t, units = "secs") # Calculate elapsed time
  if (total_t >= 3600) {t_unit <- "hours" # Determine appropriate time unit
  total_t <- total_t / 3600
  } else if (total_t >= 60) {t_unit <- "minutes"
  total_t <- total_t / 60
  } else {t_unit <- "seconds"}
  cat("Total time taken:", total_t, t_unit, "\n") # Print time
  return(ENgenes) # Final output for simple EN
}

flexiDEG.ENpartition <- function(counts_filt, metadata) {
  #No separate validation
  start_time <- Sys.time()
  bc_best <- counts_filt
  grouped<- metadata$Group
  samples <- length(grouped) 
  sample_name<- metadata$Samples #Sample names
  sample_num <- length(sample_name) # Number of samples
  user_input <- "Sample"
  batch <- c(rep(user_input,sample_num))
  n=length(unique(grouped)) 
  # Combine lines to convert 'grouped' into a factor with levels 0 to (n-1)
  groups <- factor(as.integer(factor(grouped, labels = 1:n)) - 1, levels = 0:(n-1))
  # Create a data frame with the factor variable
  xfactors <- cbind.data.frame(groups)
  dataset <- t(bc_best) # Load the data
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  
  # Create an empty vector to store the AUC values for each alpha
  num_alphas <- 21  # (0 to 1 with interval of 0.05, so 21 alpha values)
  
  # Pre-assign colors for different alpha values using a custom color palette
  colors <- c("#1F78B4", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", 
              "#6A3D9A", "#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6",
              "#FFD92F", "#FFFFB3", "#B15928", "#CCCCCC", "#000000",
              "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
              "#a65628")  # Add more colors if needed
  
  
  # Function to fit a model and calculate AUC for a given alpha value
  fit_and_evaluate_model <- function(alpha_val, x, y, foldid,num_groups,test.data) {
    tryCatch({
      # Stop execution if all samples are categorized into one group
      if (num_groups < 2) {
        stop("Samples are all categorized as one group; an Elastic Net regression cannot be performed")
      } else if (num_groups == 2) {
        # Perform cross-validated Elastic Net regression for binary classification
        cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alpha_val, foldid = foldid, family = "binomial"))
        # Determine the lambda value that minimizes cross-validation error
        best.lambda <- cvfit$lambda.min
        # Fit the final Elastic Net model using the best lambda value obtained from cross-validation
        fit <- suppressMessages(glmnet(x, y, alpha = alpha_val, lambda = best.lambda, family = "binomial"))
      }else{
        # Perform cross-validated Elastic Net regression for multinomial classification
        cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alpha_val, foldid = foldid, family = "multinomial"))
        # Determine the lambda value that minimizes cross-validation error
        best.lambda <- cvfit$lambda.min
        # Fit the final Elastic Net model using the best lambda value obtained from cross-validation
        fit <- suppressMessages(glmnet(x, y, alpha = alpha_val, lambda = best.lambda, family = "multinomial"))
      }
      # Generate predictions on test data
      x_test <- model.matrix(groups ~ ., test.data)[,-1]  # Predictor variables for test data
      probabilities <-(predict(fit, newx = x_test, type = "response", s = best.lambda))
      # Create a binary matrix for the true labels of the test data
      labels <- as.factor(test.data$groups)
      # Flatten the probabilities matrix to match the dimensions of the labels
      prob_flat <- probabilities
      if (num_groups == 2) {
        # Convert the probability matrix to a data frame for binomial classification
        prob_flat <- as.data.frame(prob_flat) #Binomial
        # Calculate the probability of the second class for binary classification
        prob_flat$s2<- 1 - prob_flat$s1
        # Rename columns to represent class labels
        colnames(prob_flat) <- c("0", "1")
        # Convert the probability data frame back to a matrix
        prob_flat<-as.matrix(prob_flat)
      }else{
        # Extract the probability matrix for multinomial classification
        prob_flat <- as.matrix(prob_flat[, , 1]) #Multinomial
      }
      
      roc_curves <- list()
      for (j in 0:(num_groups-1)) {
        # Check if there are samples with the current class label
        if (sum(labels == j) > 0) {
          # Compute ROC curve for the current class label
          suppressMessages(roc_curves[[j+1]] <- roc(labels == j, prob_flat[, j+1], levels = c(FALSE, TRUE)))
        }
      }
      
      # Calculate microaveraged TPR and FPR
      micro_tpr <- 0
      micro_fpr <- 0
      for (j in 1:length(roc_curves)) {
        if (!is.null(roc_curves[[j]])) {
          # Add sensitivities (TPR) and 1-specificities (FPR) for each class
          micro_tpr <- micro_tpr + roc_curves[[j]]$sensitivities
          micro_fpr <- micro_fpr + (1 - roc_curves[[j]]$specificities)
        }
      }
      # Calculate average TPR and FPR across all classes
      micro_tpr <- micro_tpr / length(roc_curves)
      micro_fpr <- micro_fpr / length(roc_curves)
      # Calculate AUC for the microaveraged ROC curve using the trapezoidal rule
      micro_auc <- -1*sum(diff(micro_fpr) * micro_tpr[-1] + diff(micro_fpr) * (micro_tpr[-1] - micro_tpr[-length(micro_tpr)]) / 2)
      # Extract coefficients and indices for feature selection
      Coef <- coef(fit, s = best.lambda)
      Index <- c()
      for (group in 1:num_groups) {
        if (num_groups == 2) {
          # For binary classification, extract indices of selected features
          Index <- c(Index, Coef@i[-1])
        }else{
          # For multinomial classification, extract indices of selected features for each class
          Index <- c(Index, Coef[[as.character(group - 1)]]@i[-1])
        }
      }
      
      # Return AUC, microaveraged TPR, and microaveraged FPR, EN Genes for this alpha
      return(list(auc = micro_auc, tpr = micro_tpr, fpr = micro_fpr, roc_curves = roc_curves,ENGenes=Index))
    }, error=function(e){
      cat("Skipping Iteration due to insufficient number of samples per group in Cross Validation", conditionMessage(e), "\n")
      return(NULL)  # Return NULL or an appropriate value to indicate the issue
    })
  }
  
  # Initialize lists to store results for each alpha 
  auc_values_list <- vector("list", length = num_alphas)
  micro_tpr_values_list <- vector("list", length = num_alphas)
  micro_fpr_values_list <- vector("list", length = num_alphas)
  # List to store ENgenes for top 3 alpha with partitioning and feature selection
  ENgenes_list <- vector("list", length = 3)
  # Prompt the user to enter a string
  message <- paste("Please enter the number of iterations for Partitioning Code (Minimum Suggested:", 1000, "): ", sep = " ")
  num_partiton_samples <- readline(prompt = message)
  num_partiton_samples<- as.numeric(num_partiton_samples)
  # Print message indicating the start of running iterations to select top 3 alphas
  cat("Running Iterations to Select 3 Alphas\n")
  # Select the top 3 alphas
  for (i in 1:num_alphas) {
    #Progress bars
    pb <- txtProgressBar(min = 0, max = num_alphas, initial = i, style = 3)
    setTxtProgressBar(pb, i)
    # Calculate the current alpha value
    alpha_val <- (i - 1) * 0.05
    # Perform partioning
    num_iters <-num_partiton_samples/10;   # Number of iteration samples for initial selection
    auc <- vector("numeric", length = num_iters)
    micro_tpr <- vector("list", length = num_iters)
    micro_fpr <- vector("list", length = num_iters)
    
    for (j in 1:num_iters) {
      # Generate a partitioned sample from your data
      training.samples <- dataset$groups %>% createDataPartition(p = 0.7, list = F) # Split 70% of data into training
      train.data  <- dataset[training.samples, ]
      test.data <- dataset[-training.samples, ]
      # Generate the design matrix (predictor variables) for the training data
      x <- model.matrix(groups~., train.data)[,-1] 
      # Extract the outcome variable (response variable) for the training data
      y <- train.data$groups # Outcome variable
      num_groups <- length(unique(y)) 
      # Generate fold indices for cross-validation
      foldid <- sample(1:length(y), size = length(y), replace = TRUE)
      # Fit and evaluate the model on the partitioned sample
      result  <- fit_and_evaluate_model(alpha_val,x, y, foldid,num_groups,test.data)
      if(is.null(result)){
        j<-j-1 # Skip the iteration if error throws up
      }else{
        auc[j]<-result$auc
        micro_tpr[[j]] <- result$tpr
        micro_fpr[[j]] <- result$fpr
      }
    }
    # Store the results for this alpha in the lists
    auc_values_list[[i]] <- auc
    # Convert the list of micro TPR values into a matrix and store it in the corresponding list for the current alpha
    micro_tpr_values_list[[i]] <- do.call(rbind, micro_tpr)
    # Convert the list of micro FPR values into a matrix and store it in the corresponding list for the current alpha
    micro_fpr_values_list[[i]] <- do.call(rbind, micro_fpr)
  }
  
  # Convert the lists to matrices
  auc_values <- sapply(auc_values_list, mean)
  # Initialize an empty list to store the column-wise averages
  micro_fpr_values <- list()
  # Loop through each matrix in the micro_fpr_values_list
  for (matrix in micro_fpr_values_list) {
    # Calculate the column-wise averages for the current matrix
    column_averages <- apply(matrix, 2, mean)
    # Append the column-wise averages to the column_averages_list
    micro_fpr_values  <- c(micro_fpr_values , list(column_averages))
  }
  
  # Initialize an empty list to store the column-wise averages
  micro_tpr_values <- list()
  # Loop through each matrix in the micro_fpr_values_list
  for (matrix in micro_tpr_values_list) {
    # Calculate the column-wise averages for the current matrix
    column_averages <- apply(matrix, 2, mean)
    # Append the column-wise averages to the column_averages_list
    micro_tpr_values  <- c(micro_tpr_values , list(column_averages))
  }
  # Combine the column-wise averages into matrices for micro TPR and FPR value
  micro_tpr_values <- do.call(cbind, micro_tpr_values)
  micro_fpr_values <- do.call(cbind, micro_fpr_values)
  # Identify the top 3 alphas with the highest AUC values
  top_indices <- order(auc_values, decreasing = TRUE)[1:3]
  par(mfrow = c(1, 1))
  # Plot the microaveraged ROC curve for the alpha 
  plot(NA, NA, type = "n", col = "blue", lty = 1, lwd = 2, 
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity", 
       main = "Microaveraged ROC for Top 3 Alphas")
  # Add grids to the plot
  grid()
  # Add a straight line joining (0,0) to (1,1)
  abline(a = 0, b = 1, col = "red", lty = 2)
  # Plot the microaveraged ROC curves for the top 5 alphas
  for (i in top_indices) {
    alpha_val <- (i - 1) * 0.05
    lines(micro_fpr_values[, i], micro_tpr_values[, i],
          col = colors[i], lty = 1, lwd = 2)
  }
  top_aucs<-auc_values[top_indices]
  # Add a legend to show the colors and AUC values for each alpha
  legend_text <- paste("Alpha =", (top_indices - 1) * 0.05, "AUC =", round(top_aucs, 3))
  # Add a legend to show the colors and AUC values for each alpha
  legend("bottomright", legend = legend_text, col = colors[top_indices],
         lty = 1, lwd = 2, cex = 0.8)
  
  ## Now run gene selection for these three top alphas
  k<-1
  # Print message indicating the start of running iterations for top 3 alphas
  cat("\n")
  cat("Running Iterations To Select Genes/Pathways for Top 3 Alphas\n")
  for (i in top_indices) {
    pb <- txtProgressBar(min = 1, max = 3, initial = i, style = 3) # Progress bar
    setTxtProgressBar(pb, k)
    alpha_val <- (i - 1) * 0.05
    ENgenes_list_partition <- vector("list", length = num_partiton_samples)
    for (j in 1:num_partiton_samples) {
      # Perform partionioning multiple times
      training.samples <- dataset$groups %>% createDataPartition(p = 0.7, list = F)
      # Extract training data based on the generated training samples
      train.data  <- dataset[training.samples, ]
      # Extract testing data based on the remaining samples
      test.data <- dataset[-training.samples, ]
      # Generate design matrix (predictor variables) for the training data
      x <- model.matrix(groups~., train.data)[,-1] # Predictor variables
      # Extract outcome variable (response variable) for the training data
      y <- train.data$groups # Outcome variable
      # Calculate the number of unique groups/classes in the outcome variable
      num_groups <- length(unique(y))
      # Generate fold indices for cross-validation
      foldid <- sample(1:length(y), size = length(y), replace = TRUE)
      # Fit and evaluate the model on the[ bootstrap partitioned samples
      result  <- fit_and_evaluate_model(alpha_val,x, y, foldid,num_groups,test.data)
      if(is.null(result)){
        j<-j-1
      }else{
        ENgenes_list_partition[[j]]<-result$ENGenes
      }
    }
    # Function to find elements that appear at least 80% of the time
    find_common_elements <- function(lists,alpha_val) {
      # Combine all the lists into one long vector and keep only the unqiue 
      lists<-ENgenes_list_partition
      # Unlist the combined vector and retain only the unique elements
      all_elements <- unlist(lapply(lists, unique))
      # Count the frequency of each element
      element_counts <- table(all_elements)
      # Convert the table to a data frame for easier writing to CSV
      element_counts_df <- as.data.frame(element_counts)
      # Extract the names of the elements for further processing
      elements_name <- (element_counts_df$all_elements)
      # Update the element names based on rownames of bc_best using elements_name
      element_counts_df$all_elements<- rownames(bc_best)[elements_name]
      # Sort the dataframe in descending order based on the 'Freq' colum
      sorted_df <- element_counts_df[order(-element_counts_df$Freq), ]
      # Check if sorted_df$all_elements has 5 or more elements containing "_" or a space
      if (sum(grepl("[_ ]", sorted_df$all_elements)) >= 5) {
        # If condition is met, set final_name as "PathwayFrequency_"
        final_name <- "PathwayFrequency_"
      } else {
        # Otherwise, set final_name as "GeneFrequency_"
        final_name <- "GeneFrequency_"
      }
      # Create the file name based on the value of alpha_val and final_name
      file_name <- paste(final_name, "Alpha", format(alpha_val, nsmall = 2, digits = 2), ".csv", sep = "")
      # Write the data frame to a CSV file named "element_counts.csv"
      write.csv( sorted_df , file = file_name, row.names = FALSE)
      # Find the elements that appear at least 80% of the time
      common_elements <- as.integer(names(element_counts[element_counts >= (0.8*length(lists) )]))
      # Sort the elements based on their frequency in descending order
      sorted_counts <- sort(element_counts, decreasing = TRUE)
      # Get the top 10 elements and their counts for the barplor
      top_25_elements <- rownames(bc_best)[as.numeric(rownames(sorted_counts[1:25]))]
      top_25_counts <- sorted_counts[1:25]
      # Set the margins for the plot
      par(mar = c(5, 40, 4, 2))  # c(bottom, left, top, right) - Adjust the values as needed
      # Plot the bar graph
      barplot(top_25_counts, horiz = TRUE, las = 1, names.arg = top_25_elements, main = paste("Top 25 Common Elements for Alpha:", alpha_val),
              xlab = "Frequency", col = "steelblue",cex.names = 0.7)
      return(common_elements)
    }
    ENgenes_list[[k]] <-find_common_elements(ENgenes_list_partition,alpha_val)# Get elements showing up 80% of the time
    names(ENgenes_list)[k] <- as.character(alpha_val)
    k<-k+1
  }
  #dev.off()
  end_time <- Sys.time()
  # Calculate the time difference
  total_time <- difftime(end_time, start_time, units = "secs")
  # Define conversion factors
  secs_per_min <- 60
  secs_per_hour <- 60 * 60
  
  # Determine the appropriate time unit
  if (total_time >= secs_per_hour) {
    time_unit <- "hours"
    total_time_value <- total_time / secs_per_hour
  } else if (total_time >= secs_per_min) {
    time_unit <- "minutes"
    total_time_value <- total_time / secs_per_min
  } else {
    time_unit <- "seconds"
    total_time_value <- total_time
  }
  
  # Print the result
  cat("\n")
  cat("Total time taken:", total_time_value, time_unit, "\n")
  return(ENgenes_list)
}

flexiDEG.ENbootstrap <- function(counts_filt, metadata) { # n >= 10 sample size
  start_t <- Sys.time()
  bc_best <- counts_filt
  grouped<- metadata$Group
  samples <- length(grouped) 
  
  # Extract Sample names and number
  sample_name<- metadata$Samples #Sample names
  sample_num <- length(sample_name) # Number of samples
  user_input<-"Samples"
  batch <- c(rep(user_input,sample_num))
  n=length(unique(grouped)) 
  # Combine lines to convert 'grouped' into a factor with levels 0 to (n-1)
  groups <- factor(as.integer(factor(grouped, labels = 1:n)) - 1, levels = 0:(n-1))
  # Create a data frame with the factor variable
  xfactors <- cbind.data.frame(groups)
  # Transpose the data matrix bc_best to create the dataset
  dataset <- t(bc_best) # Load the data
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groups %>% createDataPartition(p = 0.7, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groups~., train.data)[,-1] # Predictor variables
  y <- train.data$groups # Outcome variable
  num_groups <- length(unique(y)) 
  foldid <- sample(1:length(y), size = length(y), replace = TRUE)
  # Create an empty vector to store the AUC values for each alpha
  num_alphas <- 21  # (0 to 1 with interval of 0.05, so 21 alpha values)
  # Pre-assign colors for different alpha values using a custom color palette
  colors <- c("#1F78B4", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", 
              "#6A3D9A", "#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6",
              "#FFD92F", "#FFFFB3", "#B15928", "#CCCCCC", "#000000",
              "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
              "#a65628")  # Add more colors if needed
  
  
  # Function to fit a model and calculate AUC for a given alpha value
  fit_and_evaluate_model <- function(alpha_val, x, y, foldid,num_groups,test.data) {
    tryCatch({
      if (num_groups < 2) {
        stop("Samples are all categorized as one group; an Elastic Net regression cannot be performed")
      } else if (num_groups == 2) {
        # Fit an Elastic Net regression model for binomial family
        cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alpha_val, foldid = foldid, family = "binomial"))
        # Extract the value of lambda corresponding to the minimum cross-validated error
        best.lambda <- cvfit$lambda.min
        # Fit the final model using the selected lambda value
        fit <- suppressMessages(glmnet(x, y, alpha = alpha_val, lambda = best.lambda, family = "binomial"))
      }else{
        # Fit an Elastic Net regression model for multinomial family
        cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alpha_val, foldid = foldid, family = "multinomial"))
        # Extract the value of lambda corresponding to the minimum cross-validated error
        best.lambda <- cvfit$lambda.min
        # Fit the final model using the selected lambda value
        fit <- suppressMessages(glmnet(x, y, alpha = alpha_val, lambda = best.lambda, family = "multinomial"))
      }
      
      # Generate predictions on test data
      x_test <- model.matrix(groups ~ ., test.data)[,-1]  # Predictor variables for test data
      probabilities <-(predict(fit, newx = x_test, type = "response", s = best.lambda))
      # Create a binary matrix for the true labels of the test data
      labels <- as.factor(test.data$groups)
      # Flatten the probabilities matrix to match the dimensions of the labels
      if (num_groups == 2) {
        # Convert probabilities to a data frame for binomial case
        probabilities <- as.data.frame(probabilities) #Binomial
        # Calculate the second class probability for binomial case
        probabilities$s2<- 1 - probabilities$s1
        # Rename the columns 
        colnames(probabilities) <- c("0", "1")
        # Convert the probabilities data frame to a matrix
        probabilities<-as.matrix(probabilities)
      }else{
        # Convert probabilities directly to a matrix for multinomial case
        probabilities <- as.matrix(probabilities[, , 1]) #Multinomial
      }
      
      roc_curves <- list()
      for (j in 0:(num_groups-1)) {
        # Check if there are samples belonging to the current class
        if (sum(labels == j) > 0) {
          #calculate ROC curve for the current class
          suppressMessages(roc_curves[[j+1]] <- pROC::roc(labels == j, probabilities[, j+1],levels = c(FALSE, TRUE)))
        }
      }
      # Calculate microaveraged TPR and FPR
      micro_tpr <- 0
      micro_fpr <- 0
      for (j in 1:length(roc_curves)) {
        if (!is.null(roc_curves[[j]])) {
          # Calculate micro-averaged true positive rate and false positive rate
          micro_tpr <- micro_tpr + roc_curves[[j]]$sensitivities
          micro_fpr <- micro_fpr + (1 - roc_curves[[j]]$specificities)
        }
      }
      
      # Compute the average true positive rate and false positive rate
      micro_tpr <- micro_tpr / length(roc_curves)
      micro_fpr <- micro_fpr / length(roc_curves)
      # Calculate AUC for the microaveraged ROC curve using the trapezoidal rule
      micro_auc <- -1*sum(diff(micro_fpr) * micro_tpr[-1] + diff(micro_fpr) * (micro_tpr[-1] - micro_tpr[-length(micro_tpr)]) / 2)
      # Extract coefficients from the fitted model for the selected lambda
      Coef <- coef(fit, s = best.lambda)
      Index <- c()
      for (group in 1:num_groups) {
        if (num_groups == 2) {
          # Extract the index of selected variables for binomial case
          Index <- c(Index, Coef@i[-1])
        }else{
          # Extract the index of selected variables for multinomial case
          Index <- c(Index, Coef[[as.character(group - 1)]]@i[-1])
        }
      }
      
      # Return AUC, microaveraged TPR, and microaveraged FPR, EN Genes for this alpha
      return(list(auc = micro_auc, tpr = micro_tpr, fpr = micro_fpr, roc_curves = roc_curves,ENGenes=Index))
    }, error=function(e){
      cat("Skipping Iteration due to insufficient number of samples per group in Cross Validation", conditionMessage(e), "\n")
      return(NULL)  # Return NULL or an appropriate value to indicate the issue
    })
  }
  # Initialize lists to store results for each alpha and bootstrapped sample
  auc_values_list <- vector("list", length = num_alphas)
  micro_tpr_values_list <- vector("list", length = num_alphas)
  micro_fpr_values_list <- vector("list", length = num_alphas)
  # List to store ENgenes for top 3 alpha with bootstrapping and feature selection
  ENgenes_list <- vector("list", length = 3)
  # Prompt the user to enter a string
  message <- paste("Please enter the number of bootstrapped samples (Minimum Suggested:", 1000, "): ", sep = " ")
  num_bootstrap_samples <- readline(prompt = message)
  num_bootstrap_samples<- as.numeric(num_bootstrap_samples)
  cat("Running Iterations To select Top 3 Alphas\n")
  # Select the top 3 alphas
  for (i in 1:num_alphas) {
    # Start the progress bar
    pb <- txtProgressBar(min = 0, max = num_alphas, initial = i, style = 3)
    setTxtProgressBar(pb, i)
    #Calculate the alpha value
    alpha_val <- (i - 1) * 0.05
    # Perform bootstrapping with bagging
    num_iters <-num_bootstrap_samples/10   # Number of bootstrap samples  for initial alpha selection
    auc_bootstrap <- vector("numeric", length = num_iters)
    micro_tpr_bootstrap <- vector("list", length = num_iters)
    micro_fpr_bootstrap <- vector("list", length = num_iters)
    for (j in 1:num_iters) {
      # Generate a bootstrap sample from your data
      min_samples_per_group <- 4  # Adjust this value as needed
      # Initialize a variable to keep track of whether the criterion is met
      criterion_met <- FALSE
      while (!criterion_met) {
        # Perform random sampling with replacement
        boot_indices <- sample(length(y), replace = TRUE)
        # Check the number of samples in each group
        samples_per_group <- table(y[boot_indices])  # Assuming y is your vector of group labels
        # Check if the criterion is met for all groups
        if (all(samples_per_group >= min_samples_per_group)) {
          criterion_met <- TRUE
        }
      }
      # Fit and evaluate the model on the bootstrap sample
      result  <- fit_and_evaluate_model(alpha_val,x[boot_indices, ], y[boot_indices], foldid,num_groups,test.data)
      if(is.null(result)){
        j<-j-1
      }else{
        auc_bootstrap[j]<-result$auc
        micro_tpr_bootstrap[[j]] <- result$tpr
        micro_fpr_bootstrap[[j]] <- result$fpr
      }
      
    }
    ## Calculate the mean AUC across bootstrap samples for this alpha
    # Store the results for this alpha in the lists
    auc_values_list[[i]] <- auc_bootstrap
    micro_tpr_values_list[[i]] <- do.call(rbind, micro_tpr_bootstrap)
    micro_fpr_values_list[[i]] <- do.call(rbind, micro_fpr_bootstrap)
  }
  # Convert the lists to matrices
  auc_values <- sapply(auc_values_list, mean)
  # Initialize an empty list to store the column-wise averages
  micro_fpr_values <- list()
  # Loop through each matrix in the micro_fpr_values_list
  for (matrix in micro_fpr_values_list) {
    # Calculate the column-wise averages for the current matrix
    column_averages <- apply(matrix, 2, mean)
    # Append the column-wise averages to the column_averages_list
    micro_fpr_values  <- c(micro_fpr_values , list(column_averages))
  }
  # Initialize an empty list to store the column-wise averages
  micro_tpr_values <- list()
  # Loop through each matrix in the micro_fpr_values_list
  for (matrix in micro_tpr_values_list) {
    # Calculate the column-wise averages for the current matrix
    column_averages <- apply(matrix, 2, mean)
    # Append the column-wise averages to the column_averages_list
    micro_tpr_values  <- c(micro_tpr_values , list(column_averages))
  }
  micro_tpr_values <- do.call(cbind, micro_tpr_values)
  micro_fpr_values <- do.call(cbind, micro_fpr_values)
  # Identify the top 3 alphas with the highest AUC values
  top_indices <- order(auc_values, decreasing = TRUE)[1:3]
  par(mfrow = c(1, 1))
  # Plot the microaveraged ROC curve for alpha=1 to set the proper axis and labels
  plot(NA, NA, type = "n", col = "blue", lty = 1, lwd = 2, 
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity", 
       main = "Microaveraged ROC for Top 3 Alphas",margins=c(1,3))
  # Add grids to the plot
  grid()
  # Add a straight line joining (0,0) to (1,1)
  abline(a = 0, b = 1, col = "red", lty = 2)
  # Plot the microaveraged ROC curves for the top 5 alphas
  for (i in top_indices) {
    alpha_val <- (i - 1) * 0.05
    lines(micro_fpr_values[, i], micro_tpr_values[, i],
          col = colors[i], lty = 1, lwd = 2)
  }
  top_aucs<-auc_values[top_indices]
  # Add a legend to show the colors and AUC values for each alpha
  legend_text <- paste("Alpha =", (top_indices - 1) * 0.05, "AUC =", round(top_aucs, 3))
  # Add a legend to show the colors and AUC values for each alpha
  legend("bottomright", legend = legend_text, col = colors[top_indices],
         lty = 1, lwd = 2, cex = 0.8)
  ## Now run gene selection for these three top alphas
  k<-1
  cat("\n")
  cat("Running Iterations for Top 3 Alphas\n")
  for (i in top_indices) {
    pb <- txtProgressBar(min = 1, max = 3, initial = i, style = 3)
    setTxtProgressBar(pb, k)
    alpha_val <- (i - 1) * 0.05
    # Perform bootstrapping with bagging
    ENgenes_list_partition <- vector("list", length = num_bootstrap_samples)
    for (j in 1:num_bootstrap_samples) {
      min_samples_per_group <- 4  # Adjust this value as needed
      # Generate a bootstrap sample from the training data
      criterion_met <- FALSE
      while (!criterion_met) {
        # Perform random sampling with replacement
        boot_indices <- sample(length(y), replace = TRUE)
        # Check the number of samples in each group
        samples_per_group <- table(y[boot_indices])  #  y is vector of group labels
        # Check if the criterion is met for all groups
        if (all(samples_per_group >= min_samples_per_group)) {
          criterion_met <- TRUE
        }
      }
      # Fit and evaluate the model on the bootstrap sample
      result  <- fit_and_evaluate_model(alpha_val,x[boot_indices, ], y[boot_indices], foldid,num_groups,test.data)
      #result  <- fit_and_evaluate_model(alpha_val,x, y, foldid,num_groups,test.data)
      if(is.null(result)){
        j<-j-1 #Skip iteration if error comes up
      }else{
        ENgenes_list_partition[[j]]<-result$ENGenes
      }
    }
    # Function to find elements that appear at least 80% of the time
    find_common_elements <- function(lists,alpha_val) {
      # Combine all the lists into one long vector and keep only the unqiue 
      lists<-ENgenes_list_partition
      all_elements <- unlist(lapply(lists, unique))
      # Count the frequency of each element
      element_counts <- table(all_elements)
      # Convert the table to a data frame for easier writing to CSV
      element_counts_df <- as.data.frame(element_counts)
      elements_name <- (element_counts_df$all_elements)
      element_counts_df$all_elements<- rownames(bc_best)[elements_name]
      # Sort the dataframe in descending order based on the 'Freq' column
      sorted_df <- element_counts_df[order(-element_counts_df$Freq), ]
      # Check if sorted_df$all_elements has 5 or more elements containing "_" or a space
      if (sum(grepl("[_ ]", sorted_df$all_elements)) >= 5) {
        # If condition is met, set final_name as "PathwayFrequency_"
        final_name <- "PathwayFrequency_"
      } else {
        # Otherwise, set final_name as "GeneFrequency_"
        final_name <- "GeneFrequency_"
      }
      # Create the file name based on the value of alpha_val and final_name
      file_name <- paste(final_name, "Alpha", format(alpha_val, nsmall = 2, digits = 2), ".csv", sep = "")
      # Write the data frame to a CSV file named "element_counts.csv"
      write.csv( sorted_df , file = file_name, row.names = FALSE)
      # Find the elements that appear at least 80% of the time
      common_elements <- as.integer(names(element_counts[element_counts >= (0.8*length(lists) )]))
      # Sort the elements based on their frequency in descending order
      sorted_counts <- sort(element_counts, decreasing = TRUE)
      # Get the top 10 elements and their counts
      top_25_elements <- rownames(bc_best)[as.numeric(rownames(sorted_counts[1:25]))]
      top_25_counts <- sorted_counts[1:25]
      # Set the margins for the plot
      par(mar = c(5, 30, 4, 2))  # c(bottom, left, top, right) - Adjust the values as needed
      # Plot the bar graph
      barplot(top_25_counts, horiz = TRUE, las = 1, names.arg = top_25_elements, main = paste("Top 25 Common Elements for Alpha:", alpha_val),
              xlab = "Number of Times", col = "steelblue", cex.names = 0.5)
      return(common_elements)
    }
    
    ENgenes_list[[k]] <-find_common_elements(ENgenes_list_partition,alpha_val)
    names(ENgenes_list)[k] <- as.character(alpha_val)
    k<-k+1
  }
  total_t <- difftime(Sys.time(), start_t, units = "secs") # Calculate elapsed time
  if (total_t >= 3600) {t_unit <- "hours" # Determine appropriate time unit
  total_t <- total_t / 3600
  } else if (total_t >= 60) {t_unit <- "minutes"
  total_t <- total_t / 60
  } else {t_unit <- "seconds"}
  cat("\n")
  cat("Total time taken:", total_t, t_unit, "\n") # Print time
  return(ENgenes_list)
}

# Smaller functions ----

# This function is for the special case where function 2, 3, or 4 was performed on 
# GSVA data instead of gene expression to identify recurring genes from a pathway series.
flexiDEG.sharedgenes <- function(pathwaylist, counts_gsva, cutoff = 1, plot = TRUE) {
  common <- unlist(pathwaylist[rownames(counts_gsva)])
  common <- as.data.frame(table(common))
  common <- common[order(common$Freq, decreasing = T),] # Gene frequency
  common_cut <- common[common$Freq > cutoff,] # Genes with > X pathway occurrences
  if (plot == TRUE) {
    par(mar=c(3,3,3,3))# Set margins for plots
    par(mfrow=c(2,2))  # Set up a 2x2 grid for plots
    hist(common$Freq, xlab="Number of Pathways", ylab="Frequency", col="skyblue",
         main="Frequency of All Genes in Pathways")
    boxplot(common$Freq, horizontal=TRUE, xlab="Frequency", col="skyblue")
    hist(common_cut$Freq, xlab="Number of Pathways", ylab="Frequency", col="salmon", 
         main="Frequency of Cutoff Passed Genes in Pathways")
    boxplot(common_cut$Freq, horizontal=TRUE, xlab="Frequency", col="salmon")
    par(mfrow = c(1, 1)) # Reset plotting configuration
  }
  common_cut <- as.vector(common_cut[,1])
  return(common_cut)
}

flexiDEG.colors <- function(metadata) {
  color_palette <- rainbow(length(unique(metadata$Group))) # Generate group colors
  cohort_colors <- setNames(color_palette, unique(metadata$Group)) # Map colors
  colSide <- vector("character", length(metadata$Samples)) # Empty vector for colors
  for (x in 1:length(metadata$Samples)) {cohort_name <- metadata$Group[x] # Name colors
  if (cohort_name %in% names(cohort_colors)) {
    colSide[x] <- cohort_colors[cohort_name]
  } else {colSide[x] <- "yellow"}  # Color for unknown group
  }
  names(colSide) <- metadata$Group
  return(colSide)
}

flexiDEG.doublevolc <- function(counts, blank) {           # ++++ Just started working on this function
  rawpval <- apply(counts, 1, ttest_all2, 
                   grp1 = blank, grp2 = blank)
  lg_counts <- log(counts + 1, base=2)
  control <- apply(lg_counts[, blank], 1, mean) # Means in reference group
  test <- apply(lg_counts[, blank], 1, mean) # Means in test group
  logFC <- test - control # Difference between means
  results1 <- as.data.frame(cbind(logFC, rawpval))
  results1$genes <- rownames(results1)
  # results1 <- merge(results1, ensID, by = "row.names", all.x = F) # ++++ ensID object is bad
  rownames(results1) <- results1[,5]
  results1 <- results1[ , -c(1,6)]
  colnames(results1) <- c("logFC", "rawpval", "IDs", "genes") # ++++ Get rid of IDs for 1, 2, and both
  
  rawpval <- apply(counts, 1, ttest_all2, 
                   grp1 = blank, grp2 = blank)
  lg_counts <- log(counts + 1, base=2)
  control <- apply(lg_counts[, blank], 1, mean) # Means in reference group
  test <- apply(lg_counts[, blank], 1, mean) # Means in test group
  logFC <- test - control # Difference between means
  results2 <- as.data.frame(cbind(logFC, rawpval))
  results2$genes <- rownames(results2)
  # results2 <- merge(results2, ensID, by = "row.names", all.x = F) # ++++ ensID object is bad
  rownames(results2) <- results2[,5]
  results2 <- results2[ , -c(1,6)]
  colnames(results2) <- c("logFC", "rawpval", "IDs", "genes")
  
  results_both <- cbind(results1, results2)
  results_both <- results_both[ , -c(3,4)]
  results_both <- as.data.frame(results_both)
  colnames(results_both) <- c("logFC_S", "rawpval_S", "logFC_H", "rawpval_H", "IDs", "genes")
  results_both$dist <- round((abs(results_both$logFC_S)) + (abs(results_both$logFC_H)), 2)
  VolcanoTop <- results_both[g_cutoff1F, ]
  bothVolcanoTops <- subset(VolcanoTop, abs(logFC_S)>0.2 & abs(logFC_H)>0.25 & (abs(logFC_S) + abs(logFC_H))>0.4)
  
  par(mfrow=c(1,1), mar=c(4,4,2,1), cex=1.0, cex.main=1.4, cex.axis=1.0, cex.lab=1.2)
  ggplot() + # When saving this plot, W:600, H:450
    theme_classic() + xlim(-0.7, 0.7) + ylim(-0.7, 0.7) + # Format axes
    xlab = "STx Log[2] Fold Change" + ylab = "HTx Log[2] Fold Change" + 
    geom_hline(yintercept = 0, color = "black") + geom_vline(xintercept = 0, color = "black") + 
    geom_point(data = results_both, aes(x = logFC_S, y = logFC_H, color = dist), shape = 1, size = 0.8) + 
    scale_color_gradient(low = "grey90", high = "grey40") +
    geom_point(data = VolcanoTop, aes(x = logFC_S, y = logFC_H), color = "grey30") +
    geom_point(data = bothVolcanoTops, aes(x = logFC_S, y = logFC_H), color = "black", fill = "deeppink1", 
               size = 3, shape = 21, stroke = 1.4) + 
    # geom_text_repel(data = results_both[g_allomapnames,], aes(x = logFC_S, y = logFC_H, label = genes), 
    # color = "blue", min.segment.length = 0) +
    geom_text_repel(data = bothVolcanoTops, aes(x = logFC_S, y = logFC_H, label = genes), color = "deeppink4", 
                    min.segment.length = 0, max.overlaps = 30)
  par(mar = c(4,4,4,4), oma = c(1,1,1,1)) # Reset margins
  g_bothVolcanoTops <- bothVolcanoTops[,c(5,6)]
  counts_volc <- counts[g_bothVolcanoTops,]
  return(counts_volc)
}

flexiDEG.ENplots <- function(data, ENlist, colors, colors2) {
  require("ggbiplot")
  require("ggplot2")
  plots <- lapply(1:20, function(alpha) {
    g_EN <- unique(rownames(data)[unlist(ENlist$ENgenes[alpha])])
    counts_EN <- na.omit(data[g_EN, ])
    if (length(counts_EN[,1]) <= 1) {
      paste("Alpha =", names(ENlist$ENgenes[alpha]), "not used")
    } else {
      ggbiplot(prcomp(t(counts_EN), scale.=T), ellipse=T, groups=names(colors), var.axes=F, 
               var.scale=1, circle=T) + theme_classic() + 
        ggtitle(paste("a =", names(ENlist$ENgenes[alpha]))) + 
        geom_point(size=3, color=colors) + scale_color_manual(name="Group", values=colors2) +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none")
    }})
  return(plots)
}

ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE, 
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = muted('red'))
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  # Label the variable axes
  if(var.axes) {
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar, y = yvar, 
                    angle = angle, hjust = hjust), 
                color = 'darkred', size = varname.size)
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}