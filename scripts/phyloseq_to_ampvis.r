is_ampvis = 
  function(data) {
    if (!inherits(data, "ampvis2")) {
      stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis2 functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)", call. = FALSE)
    }
    invisible(TRUE)
  }

normaliseTo100 =
  function(data) {
    ### Data must be in ampvis2 format
    is_ampvis2(data)
    
    if (!abundAreCounts(data)) {
      warning("The data has already been normalised. Setting normalise = TRUE (the default) will normalise the data again and the relative abundance information about the original data of which the provided data is a subset will be lost.", call. = FALSE)
    }
    # normalise each sample to sample totals, skip samples with 0 sum to avoid NaN's
    tmp <- data$abund[, which(colSums(data$abund) != 0), drop = FALSE]
    if (nrow(tmp) == 1L) {
      # apply returns a vector and drops rownames if only 1 row, therefore set to 100 instead
      tmp[1L, ] <- 100L
    } else if (nrow(tmp) > 1L) {
      tmp <- as.data.frame(apply(tmp, 2, function(x) {
        x / sum(x) * 100
      }))
    }
    data$abund[, which(colSums(data$abund) != 0)] <- tmp
    attributes(data)$normalised <- TRUE
    return(data)
  }

phyloseq_to_ampvis2 = 
  function(physeq) {
    #check object for class
    if(!any(class(physeq) %in% "phyloseq"))
      stop("physeq object must be of class \"phyloseq\"", call. = FALSE)
    
    #ampvis2 requires taxonomy and abundance table, phyloseq checks for the latter
    if(is.null(physeq@tax_table))
      stop("No taxonomy found in the phyloseq object and is required for ampvis2", call. = FALSE)
    
    #OTUs must be in rows, not columns
    if(phyloseq::taxa_are_rows(physeq))
      abund <- as.data.frame(phyloseq::otu_table(physeq)@.Data)
    else
      abund <- as.data.frame(t(phyloseq::otu_table(physeq)@.Data))
    
    #tax_table is assumed to have OTUs in rows too
    tax <- phyloseq::tax_table(physeq)@.Data
    
    #merge by rownames (OTUs)
    otutable <- merge(
      abund,
      tax,
      by = 0,
      all.x = TRUE,
      all.y = FALSE,
      sort = FALSE
    )
    colnames(otutable)[1] <- "OTU"
    
    #extract sample_data (metadata)
    if(!is.null(physeq@sam_data)) {
      metadata <- data.frame(
        phyloseq::sample_data(physeq),
        row.names = phyloseq::sample_names(physeq), 
        stringsAsFactors = FALSE, 
        check.names = FALSE
      )
      
      #check if any columns match exactly with rownames
      #if none matched assume row names are sample identifiers
      samplesCol <- unlist(lapply(metadata, function(x) {
        identical(x, rownames(metadata))}))
      
      if(any(samplesCol)) {
        #error if a column matched and it's not the first
        if(!samplesCol[[1]])
          stop("Sample ID's must be in the first column in the sample metadata, please reorder", call. = FALSE)
      } else {
        #assume rownames are sample identifiers, merge at the end with name "SampleID"
        if(any(colnames(metadata) %in% "SampleID"))
          stop("A column in the sample metadata is already named \"SampleID\" but does not seem to contain sample ID's", call. = FALSE)
        metadata$SampleID <- rownames(metadata)
        
        #reorder columns so SampleID is the first
        metadata <- metadata[, c(which(colnames(metadata) %in% "SampleID"), 1:(ncol(metadata)-1L)), drop = FALSE]
      }
    } else
      metadata <- NULL
    
    #extract phylogenetic tree, assumed to be of class "phylo"
    if(!is.null(physeq@phy_tree)) {
      tree <- phyloseq::phy_tree(physeq)
    } else
      tree <- NULL
    
    #extract OTU DNA sequences, assumed to be of class "XStringSet"
    if(!is.null(physeq@refseq)) {
      #convert XStringSet to DNAbin using a temporary file (easiest)
      fastaTempFile <- tempfile(pattern = "ampvis2_", fileext = ".fa")
      Biostrings::writeXStringSet(physeq@refseq, filepath = fastaTempFile)
    } else
      fastaTempFile <- NULL
    
    #load as normally with amp_load
    ampvis2::amp_load(
      otutable = otutable,
      metadata = metadata,
      tree = tree,
      fasta = fastaTempFile
    )
  }