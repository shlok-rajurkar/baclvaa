#This function takes case control data as an input and outputs a nested R data frame including in this order
#Each given Loci
#Each population given
#The amino acid, codon, and allele output from the BIGDAWG package
#If applicable, gives a "builtout" object which shows the given allele OR compared 
#To the given ORs of the associated silent mutations codon positions in that given population and locus
#Inputs must be Case control Data
#First col must contain subject ID
#Second col must contain phenotype
#Columns after that must contain HLA loci for both chromosomes
#This would be in form of DRB1_1 and DRB1_2
#Check HLAtoools::HLAgazeteer$nuc to see valid HLA loci

codon_associationATF <- function(rawData) {
  # Pre-compute constants outside loops
  codon_list <- list(
    "A" = A_Codons,
    "B" = B_Codons, 
    "C" = C_Codons,
    "DPA1" = DPA1_Codons,
    "DPB1" = DPB1_Codons, 
    "DQA1" = DQA1_Codons, 
    "DQB1" = DQB1_Codons,
    "DRB1" = DRB1_Codons
  )
  
  valid_loci <- HLAtools::HLAgazetteer$nuc
  
  # Pre-allocate list with proper size
  full_ana <- vector("list", length(valid_loci))
  names(full_ana) <- valid_loci
  
  # Process each dataset
  for (i in seq_along(rawData)) {
    current_data <- rawData[[i]]
    current_data_name <- names(rawData)[i]
    
    # Extract all locus columns at once using vectorized operations
    odd_cols <- seq(3, ncol(current_data), by = 2)
    locus_names <- vapply(colnames(current_data)[odd_cols], 
                          function(x) strsplit(x, "_", fixed = TRUE)[[1]][1], 
                          character(1))
    
    # Validate all loci at once
    invalid_loci <- !locus_names %in% valid_loci
    if (any(invalid_loci)) {
      stop(paste("Invalid loci:", paste(locus_names[invalid_loci], collapse = ", ")))
    }
    # Process each locus
    for (j in seq_along(locus_names)) {
      locus_name <- locus_names[j]
      
      # Check if locus has both cases and controls
      locus_check <- add_locus_edited(
        locus = locus_name, 
        codon = codon_list[[locus_name]], 
        data = current_data, 
        phenotype_exist = FALSE
      )
      
      if (!locus_check[[1]]) {
        warning(paste("Locus", locus_name, "in", current_data_name, 
                      "does not contain both cases and controls. Skipping."))
        next
      }
      
      # Initialize locus list if needed
      if (is.null(full_ana[[locus_name]])) {
        full_ana[[locus_name]] <- list()
      }
      
      # Process locus data
      full_ana[[locus_name]][[current_data_name]] <- 
        process_locus_ATF_edited(locus_check[[2]], locus_name, codon_list[[locus_name]])
    }
  }
  
  # Remove empty elements efficiently
  full_ana[lengths(full_ana) > 0]
}

# In process_locus_ATF_edited, pass codon_data through:
process_locus_ATF_edited <- function(data_list, locus_name, codon_data) {
  locus <- vector("list", length(data_list))
  names(locus) <- names(data_list)
  
  needs_amino <- locus_name %in% c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
  
  for (i in seq_along(data_list)) {
    set <- data_list[[i]]
    name <- names(data_list)[i]
    
    raw <- BIGDAWG(Data = set, HLA = FALSE, Return = TRUE, Run.Tests = "L")
    sig <- get_sig(raw)
    
    locus[[name]] <- list(
      "Raw Analysis" = raw,
      "Significant positions" = sig
    )
    
    if (i == 1 && needs_amino) {
      raw_a <- BIGDAWG(Data = set, HLA = TRUE, Return = TRUE, Run.Tests = "A")
      sig_a <- get_sig_a(raw_a)
      
      locus[["Amino"]] <- list(
        "Raw Analysis" = raw_a,
        "Significant positions" = sig_a
      )
    }
    
    if (i == 2) {
      locus <- process_builtout_analysis(locus, data_list, codon_data)  # pass codon_data
    }
  }
  
  locus
}

# Extract builtout processing to separate function for clarity
# Update process_builtout_analysis to accept and use codon_data:
process_builtout_analysis <- function(locus, data_list, codon_data) {
  
  sig_codons <- locus[["Codon"]]$`Significant positions`
  sig_allele <- locus[["Allele"]]$`Significant positions`
  same_aminos_pos <- same_aminos(sig_codons)
  
  if (length(same_aminos_pos) > 0 && nrow(sig_allele) > 0 && nrow(sig_codons) > 0) {
    
    builtout <- sig_allele_amino_pos(
      sig_allele   = sig_allele,
      sig_codon_OR = sig_codons,
      amino_pos    = same_aminos_pos,
      loci_codon   = codon_data
    )
    
    if (!is.null(builtout)) {
      locus[["Builtout"]] <- builtout
    }
    
    locus[["Builtout"]] <- builtout
  }
  
  locus
}

add_locus_edited <- function(locus, codon, data, phenotype_exist) {
  # Find locus columns once for given locus
  locus_cols <- which(grepl(paste0("^", locus), colnames(data)))
  
  if (length(locus_cols) < 2) {
    stop(paste("Insufficient columns found for locus:", locus))
  }
  
  # Process both conversions
  bigdawg_locus <- bigdawg_converter(
    codon, locus_cols[1], locus_cols[2], data, stringent = FALSE
  )
  
  #gives table with subjects, phenotypes, and locus columns
  
  bigdawg_locus_allele <- phenotype_bigdawg(
    Rawdata = data, 
    col1 = locus_cols[1], 
    col2 = locus_cols[2], 
    name1 = paste0(locus, "_1"), 
    name2 = paste0(locus, "_2"), 
    phenotype = phenotype_exist
  )
  
  locus_list <- list("Allele" = bigdawg_locus_allele, "Codon" = bigdawg_locus)
  list(both_phenotypes(locus_list), locus_list)
}

bigdawg_analysis_locus <- function(locus, locus_name) {
  process_locus_ATF_edited(locus, locus_name)
}

both_phenotypes <- function(data) {
  # checks for atleast one case and control in dataset
  allele_data <- data[[1]]
  case <- FALSE
  control <- FALSE
  for (i in 1:nrow(allele_data)) {
    if (allele_data[i, 2] == 0 && !is.na(allele_data[i, 3])){
      control <- TRUE
    }
    if (allele_data[i, 2] == 1 && !is.na(allele_data[i, 3])){
      case <- TRUE
    }
  }
  case && control
}




