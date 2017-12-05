#' @title Normalize protein abundances by \code{NSAF}
#'
#' @description normalizes proteome data by a \code{Normalized Spectral Abundance Factor}
#' using the protein length as specified by the \code{length} parameter in the
#' \code{\link[phyloseq]{tax_table}} of a phyloseq object
#' The input \code{\link[phyloseq]{otu_table}} should contain proteomic count or intensity values from a mass spec.
#'
#' @param ps A \code{phyloseq} object that contains a \code{\link[phyloseq]{tax_table}}.
#' The \code{tax_table} must contain protein length as specified with the \code{length} parameter.
#'
#' @param length A column in the \code{tax_table} containing the protein length.

#' @return This function returns the input \code{phyloseq} object with a the  \code{otu_table()} normalized
#' by \code{NSAF}.
#' The Normalized Spectral Abundance Factor was
#' proposed by \href{http://pubs.acs.org/doi/pdf/10.1021/pr060161n}{Zybailov \emph{et al.} (2006)}.
#' NSAF for a protein `k` is the number of spectral counts (`SpC`, the total number of MS/MS spectra)
#' identifying a protein, `k`, divided by the protein’s length (`L`),
#' divided by the sum of `SpC`/`L` for all `N` proteins in the experiment.
#'
#' @examples
#' ps <- add_nsaf(ps)

add_nsaf=function(ps, length){
  if(ps@otu_table@taxa_are_rows == TRUE){
    mat <- (otu_table(ps))
  }else{
    mat <- t((otu_table(ps)))
  }
  prot_len <- as.numeric(tax_table(ps)[,length])
  mat <- mat/prot_len
  mat <- mat/rowSums(mat)
  otu_table(ps) <- mat
  return(ps)
}



#' @title Merge Protein Replicates
#'
#' @description This function merges replicate samples from proteomic data
#'
#' @param ps A \code{phyloseq} object that contains a variable in \code{\link[phyloseq]{sample_data}}
#' denoting which samples are replicates of one another.
#'
#' @param repID A column in the \code{sample_data} that identifies which samples are replicates of each other and should be merged.
#' All replicates should have the same value in this column.
#'
#' @param minAbund (optional) A column in the \code{sample_data} specifying which samples are replicates to be merged.
#' This is the minimum value a protein in a set of replicates has to have for it to be considered "real"
#' Default value = 10.
#'
#' @param maxZeros A column in the \code{sample_data} specifying which samples are replicates to be merged.
#' Default value = 0.5
#'
#' @return This function returns the input \code{phyloseq} object with the sample replicates merged.
#'
#' @examples
#' ps <- merge_protein_replicates(ps, repID="Subject_ID", minAbund=10, maxZeros = 0.5)

merge_protein_replicates=function(ps, repID, minAbund=10, maxZeros = 0.5){

  # Enforce inclusion of non-optional arguments
  if(missing(ps)){print("phyloseq object not specified")}
  if(missing(repID)){print("specify which sample_data() column contains the replicate names")}

  # set function:
  evalProts=function(x){
    if(max(x) <= minAbund){
      x <- 0
    }else{
      percZeros <- 1 - sum(x != 0)/length(x)
      ifelse(percZeros <= maxZeros, x <- mean(x[which(x != 0)]), x <- 0)
    }
  }

  df <- psmelt(ps)

  df <- aggregate(as.formula(paste0("Abundance ~ OTU+",repID)), df, evalProts) # summary stats
  df <- df[with(df, order(df[,1], df[,3])), ] # sort df

  df <- dcast(df, as.formula(paste0("OTU ~", repID)), value.var= "Abundance") # re-cast to count matrix
  rownames(df) <- df$OTU; df$OTU <- NULL
  df <- df[match(taxa_names(ps), rownames(df)),] # re-order matrix

  # make a new, merged phyloseq object
  ps_i <- suppressWarnings(merge_samples(ps, repID))
  # replace the merged OTU table with the one you just made
  otu_table(ps_i) <- otu_table(df, taxa_are_rows = TRUE)
  tryCatch(
    ps_i <- prune_taxa(taxa_sums(ps_i) > 0, ps_i), # remove taxa that have zero abundance
    warning=function(w) {print("The settings you have specified do not return any proteins
                               with non-zero abundance. Try relaxing your settings.")},
    error = function(e) {print("The settings you have specified do not return any proteins
                               with non-zero abundance. Try relaxing your settings.")}
  )
  print(paste0(nsamples(ps_i), " samples returned (out of ", nsamples(ps), ")
               and ", ntaxa(ps_i), " proteins returned (out of ",ntaxa(ps), ")
               after merging on ", repID , " with minAbund = ", minAbund, "
               and maxZeros = ", maxZeros,"."))

  return(ps_i)
}

