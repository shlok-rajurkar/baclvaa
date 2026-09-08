#Input HLAalignments$codon$loci
#Outputs a similar alignment where all exons and 
#indels are removed. Leaving only each codon position
#with the three nucleotides in that position
codon_converter <- function(A){
  nuc <- A[, 5:ncol(A)]
  #Data set with only the nucleotides of the given Locus
  names <- A[1:4]
  #Data set with only the names of the alleles within the same locus
  A_cols <- colnames(nuc)
  #All column names of the differently nucleotides
  #This will include the separate codons under which each nucleotide belongs
  codon_set <- data.frame(stringsAsFactors = FALSE)
  progress <- 0
  #Iterate through each row in the alignment of the given locus
  for (j in 1:nrow(A)) {
    #Update progress to user
    if (round(j/nrow(A)*100) != progress) {
      progress <- round(j/nrow(A)*100)
      print(paste0(progress, "%"))
    }
    #Set a temporary data set. This will act as one row of the total data set
    temp_set <- data.frame(matrix(nrow = 1, ncol = 0), stringsAsFactors = FALSE)
    #Iterate through each column in a given row
    for (i in 1:length(A_cols)) {

      #If i == 1, automatically setting that value as the base column(The codon I am encoding for)
      if (i == 1) {
        base_col <- A_cols[i]
        temp_name <- nuc[j, i] #Represents the current codon
        
        #Checking for either if a column name is the start of a new codon and if that codon is
        #only one more than the previous, or checking if that new codon is 1 (since 1 is)
        #two more than the previous as there is no 0
      } else if (!is.na(as.numeric(A_cols[i]))) {
        if ((!grepl("\\.", A_cols[i]) && as.numeric(A_cols[i]) == as.numeric(base_col) + 1) | (!grepl("\\.", A_cols[i]) && as.numeric(A_cols[i]) - 1 == 0)) { 
          #Since a new codon is started, add the previous codon to the row
          temp_set[[base_col]] <- temp_name
          #Start process of creating a new codon
          base_col <- A_cols[i]
          temp_name <- nuc[j, i]
        } else if (!grepl("\\.", nuc[j,i]) && !grepl("\\|", nuc[j,i])){
          temp_name <- paste0(temp_name, nuc[j, i]) #Add nucleotide to codon
        }
          #Check if value contains a nucleotide or if it is an insertion/unknown
      } 
      
    }
    codon_set <- rbind(codon_set, temp_set) #Add new row to the original set of codons
    
  }
  Full_Dataset_A <- cbind(names, codon_set)  #Add names back to the set
  
  
  return(Full_Dataset_A)
}

bigdawg_converter <- function(Full_Dataset_A, col1, col2, Rawdata, stringent = TRUE) {
  # takes in a table of subjects and their allele labels to output codons
  phenotype <- Rawdata[2]
  study_locus <- Rawdata[c(1, col1:col2)] # Modify data set to only include our given locus
  study_proc <- data.frame(matrix(nrow = 1, ncol = 0))
  
  for (i in 1:nrow(study_locus)) { # Iterate through each row of data
    person_1 <- study_locus[i, ]
    allele_1 <- person_1[[2]] # Initialize the allele of the first chromosome
    allele_2 <- person_1[[3]] # Initialize the allele of the second chromosome
    if (stringent == TRUE) {
      if (is.na(allele_1) | is.na(allele_2)) {
        next # If the data is NA, skip
      }
    }
    
    
    # If allele name includes the generic version of a silent mutation (01) it's removed
    # 02:01:01 -> 02:01
    if (length(strsplit(allele_1, ":")[[1]]) > 2 && substring(allele_1, nchar(allele_1)) == "1") {
      positions <- gregexpr(":", allele_1)[[1]]
      allele_1 <- substring(allele_1, 1, positions[2] - 1)
    }
    
    # If allele name is missing the first 0, it is added 
    # 7:01:02 -> 07:01:02
    if(!is.na(allele_1)){
      if (nchar(strsplit(allele_1, ":")[[1]][1]) < 2) {
        allele_1 <- paste0("0", allele_1)
      }
    }
    # Repeat for allele_2
    if (length(strsplit(allele_2, ":")[[1]]) > 2 && substring(allele_2, nchar(allele_2)) == "1") {
      positions <- gregexpr(":", allele_2)[[1]]
      allele_2 <- substring(allele_2, 1, positions[2] - 1)
    }
    if(!is.na(allele_2)){
      if (nchar(strsplit(allele_2, ":")[[1]][1]) < 2) {
        allele_2 <- paste0("0", allele_2)
      }
    }
    # Check to see if current allele name is found within the alignment
    if (nrow(Full_Dataset_A[Full_Dataset_A$allele == allele_1, ]) != 1) {
      # If not, add an ":01" to encode for the silent mutation 
      sequence_1 <- Full_Dataset_A[grepl(paste0(allele_1, ":01"), Full_Dataset_A$allele), ] # Use grepl to get all cases
      sequence_1 <- sequence_1[1, ] # Pull from the first one
    } else {
      # If current allele name is found within alignment
      sequence_1 <- Full_Dataset_A[Full_Dataset_A$allele == allele_1, ]
    }
    
    # Repeat for allele_2
    if (nrow(Full_Dataset_A[Full_Dataset_A$allele == allele_2, ]) != 1) {
      sequence_2 <- Full_Dataset_A[grepl(paste0(allele_2, ":01"), Full_Dataset_A$allele), ]
      sequence_2 <- sequence_2[1, ]
    } else {
      sequence_2 <- Full_Dataset_A[Full_Dataset_A$allele == allele_2, ]
    }
    if (stringent == TRUE) {
      if (any(is.na(sequence_1)) | any(is.na(sequence_2))){
        next
      }
    }
    
    study_temp <- data.frame(matrix(nrow = 1, ncol = 0))
    
    # Iterate through sequences adding each sequence together with the corresponding chromosome at the end
    for (j in 1:ncol(sequence_1)) {
      # Adds codon to table according to codon number 
      study_temp[paste0(colnames(sequence_1)[j], "_1")] <- sequence_1[j]
      study_temp[paste0(colnames(sequence_2)[j], "_2")] <- sequence_2[j]
    }
    
    study_temp <- cbind(person_1[1], study_temp) # Add subject ID to the sequence
    study_proc <- rbind(study_proc, study_temp)  # Add new row to the data set
  }
  
  subject_ID <- study_proc[1]
  study_cod <- study_proc[10:ncol(study_proc)] # Create new data set with only codons
  
  print(nrow(subject_ID))
  print(nrow(phenotype))
  print(nrow(study_cod))
  
  # Create new data set according to BIGDAWG input needs
  study_bigdawg <- cbind(subject_ID, phenotype, study_cod)
  
  return(study_bigdawg)
}

bigdawg_allele_converter <- function(Full_Dataset_A, col1, col2, Rawdata, stringent = TRUE) {
  
  study_locus <- Rawdata[c(1, col1:col2)] # Modify data set to only include our given locus
  study_proc <- data.frame(matrix(nrow = 1, ncol = 0))
  
  for (i in 2:nrow(study_locus)) { # Iterate through each row of data
    person_1 <- study_locus[i, ]
    allele_1 <- person_1[[2]] # Initialize the allele of the first chromosome
    allele_2 <- person_1[[3]] # Initialize the allele of the second chromosome
    if (stringent == TRUE) {
      if (is.na(allele_1) | is.na(allele_2)) {
        next # If the data is NA, skip
      }
    } else {
      if (is.na(allele_1) && is.na(allele_2)) {
        next # If the data is NA, skip
      }
    }
    
    
    # If allele name includes the generic version of a silent mutation (01) it's removed
    # 02:01:01 -> 02:01
    if (length(strsplit(allele_1, ":")[[1]]) > 2 && substring(allele_1, nchar(allele_1)) == "1") {
      positions <- gregexpr(":", allele_1)[[1]]
      allele_1 <- substring(allele_1, 1, positions[2] - 1)
    }
    
    # If allele name is missing the first 0, it is added 
    # 7:01:02 -> 07:01:02
    if(!is.na(allele_1)){
      if (nchar(strsplit(allele_1, ":")[[1]][1]) < 2) {
        allele_1 <- paste0("0", allele_1)
      }
    }
    # Repeat for allele_2
    if (length(strsplit(allele_2, ":")[[1]]) > 2 && substring(allele_2, nchar(allele_2)) == "1") {
      positions <- gregexpr(":", allele_2)[[1]]
      allele_2 <- substring(allele_2, 1, positions[2] - 1)
    }
    if(!is.na(allele_2)){
      if (nchar(strsplit(allele_2, ":")[[1]][1]) < 2) {
        allele_2 <- paste0("0", allele_2)
      }
    }
    # Check to see if current allele name is found within the alignment
    if (nrow(Full_Dataset_A[Full_Dataset_A$allele == allele_1, ]) != 1) {
      # If not, add an ":01" to encode for the silent mutation 
      sequence_1 <- Full_Dataset_A[grepl(paste0(allele_1, ":01"), Full_Dataset_A$allele), ] # Use grepl to get all cases
      sequence_1 <- sequence_1[1, ] # Pull from the first one
    } else {
      # If current allele name is found within alignment
      sequence_1 <- Full_Dataset_A[Full_Dataset_A$allele == allele_1, ]
    }
    
    # Repeat for allele_2
    if (nrow(Full_Dataset_A[Full_Dataset_A$allele == allele_2, ]) != 1) {
      sequence_2 <- Full_Dataset_A[grepl(paste0(allele_2, ":01"), Full_Dataset_A$allele), ]
      sequence_2 <- sequence_2[1, ]
    } else {
      sequence_2 <- Full_Dataset_A[Full_Dataset_A$allele == allele_2, ]
    }
    if (stringent == TRUE) {
      if (any(is.na(sequence_1)) | any(is.na(sequence_2))){
        next
      }
    }
    
    study_temp <- data.frame(matrix(nrow = 1, ncol = 0))
    
    # Iterate through sequences adding each sequence together with the corresponding chromosome at the end
    for (j in 1:ncol(sequence_1)) {
      study_temp[paste0(colnames(sequence_1)[j], "_1")] <- sequence_1[j]
      study_temp[paste0(colnames(sequence_2)[j], "_2")] <- sequence_2[j]
    }
    
    study_temp <- cbind(person_1[1], study_temp) # Add subject ID to the sequence
    study_proc <- rbind(study_proc, study_temp)  # Add new row to the data set
  }
  
  
  alleles <- study_proc[, 4:5]
  
  subject_ID <- study_proc[1]
  phenotype <- data.frame(matrix(nrow = 1, ncol = 0))
  study_cod <- study_proc[10:ncol(study_proc)] # Create new data set with only codons
  
  for (val in subject_ID[[1]]) {
    phen_temp <- data.frame(matrix(nrow = 1, ncol = 0))
    # Check if row is a case or a control
    if (substring(val, 3, 3) == "C") {
      phen_temp[["phenotype"]] <- 0
    } else {
      phen_temp[["phenotype"]] <- 1
    }
    phenotype <- rbind(phenotype, phen_temp)
  }
  
  # Create new data set according to BIGDAWG input needs
  study_bigdawg <- cbind(subject_ID, phenotype, alleles, study_cod)
  
  return(study_bigdawg)
}


bigdawg_DRBX_converter <- function(DRB3_Codons, DRB4_Codons, DRB5_Codons, Rawdata, col1, col2) {
  study_A <- Rawdata[c(1, col1:col2)]
  
  DRBX_full <- data.frame(matrix(nrow = 1, ncol = 0))
  #Iterate through every row of the raw data
  for (i in 2:nrow(study_A)) {
    person_1 <- study_A[i, ] #Initialize ith row
    if (!is.na(person_1[[2]]) && !is.na(person_1[[3]])) { #Check row has data
      #Check what combination of DRBX person has
      if (substring(person_1[[2]], 1, 1) == "3") {
        if (substring(person_1[[3]], 1, 1) == "3") {
          #Pass data to sequencer which will return the given sequence of the persons two alleles
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "4") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "5") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "0") {
          row <- sequence_DRBX(DRBX_1 = DRB3_Codons, allele_1 = substring(person_1[[2]], 2), name = person_1[1])
        }
        
      } else if (substring(person_1[[2]], 1, 1) == "4") {
        
        if (substring(person_1[[3]], 1, 1) == "3") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "4") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "5") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "0") {
          row <- sequence_DRBX(DRBX_1 = DRB4_Codons, allele_1 = substring(person_1[[2]], 2), name = person_1[1])
        }
        
      } else if (substring(person_1[[2]], 1, 1) == "5") {
        
        if (substring(person_1[[3]], 1, 1) == "3") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "4") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "5") {
          row <- sequence_DRBX(DRB3_Codons, DRB3_Codons, substring(person_1[[2]], 2), substring(person_1[[3]], 2), person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "0") {
          row <- sequence_DRBX(DRBX_1 = DRB5_Codons, allele_1 = substring(person_1[[2]], 2), name = person_1[1])
        }
        
      } else if (substring(person_1[[2]], 1, 1) == "0") {
        
        if (substring(person_1[[3]], 1, 1) == "3") {
          row <- sequence_DRBX(DRBX_2 = DRB3_Codons, allele_2 = substring(person_1[[3]], 2), name = person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "4") {
          row <- sequence_DRBX(DRBX_2 = DRB4_Codons, allele_2 = substring(person_1[[3]], 2), name = person_1[1])
          
        } else if (substring(person_1[[3]], 1, 1) == "5") {
          row <- sequence_DRBX(DRBX_2 = DRB5_Codons, allele_2 = substring(person_1[[3]], 2), name = person_1[1])
        }else if (substring(person_1[[3]], 1, 1) == "0") {
          row <- sequence_DRBX(name = person_1[1])
        }
        
      }
      
      DRBX_full <- rbind(DRBX_full, row)
    }
  }
  
  subject_ID <- DRBX_full[1]
  phenotype <- data.frame(matrix(nrow = 1, ncol = 0))
  study_cod <- DRBX_full[10:ncol(DRBX_full)]
  
  
  for(val in subject_ID[[1]]) {
    phen_temp <- data.frame(matrix(nrow = 1, ncol = 0))
    #Check if row is case or control, and assign them a phenotype
    if (substring(val, 3, 3) == "P") {
      phen_temp[["phenotype"]] <- 0
    } else phen_temp[["phenotype"]] <- 1
    phenotype <- rbind(phenotype, phen_temp)
  }
  #Bind subject ID, phenotype, and codon sequence according to BIGDAWG input needs
  bigdawg_DRBX <- cbind(subject_ID, phenotype, study_cod)
  return(bigdawg_DRBX)
}


sequence_DRBX <- function(DRBX_1 = NA, DRBX_2 = NA, allele_1 = NA, allele_2 = NA, name){
  
  
  
  col_0_names <- c("locus", "allele", "trimmed_allele", "allele_name")
  
  # Generate numbered column names from -29 to 237, skipping 0
  numbered_names <- as.character(c(-29:-1, 1:237))
  
  # Combine them
  col_0_names <- c(col_0_names, numbered_names)
  
  # Create the data frame with 1 row and all values set to "000"
  zeros <- data.frame(matrix("000", nrow = 1, ncol = length(col_0_names)), stringsAsFactors = FALSE)
  
  # Assign the column names
  colnames(zeros) <- col_0_names
  
  #Go through same process as seen above
  if (!is.na(allele_1)){
    if( length(strsplit(allele_1, ":")[[1]]) > 2 && substring(allele_1, nchar(allele_1)) == "1"){
      positions <- gregexpr(":", allele_1)[[1]]
      allele_1 <- substring(allele_1, 1, positions[2] - 1)
    }
    if(nchar(strsplit(allele_1, ":")[[1]][1]) < 2) {
      allele_1 <- paste0("0", allele_1)
    }
  }
  
  if (!is.na(allele_2)){
    if( length(strsplit(allele_2, ":")[[1]]) > 2 && substring(allele_2, nchar(allele_2)) == "1"){
      positions <- gregexpr(":", allele_2)[[1]]
      allele_2 <- substring(allele_2, 1, positions[2] - 1)
    }
    if(nchar(strsplit(allele_2, ":")[[1]][1]) < 2) {
      allele_2 <- paste0("0", allele_2)
    }
  }
  
  
  if (!is.na(allele_1)){
    if (nrow(DRBX_1[DRBX_1$allele == allele_1, ]) < 1) {
      sequence_1 <- DRBX_1[grepl(paste0(allele_1, ":01"), DRBX_1$allele), ]
      sequence_1 <- sequence_1[1, ]
    } else {
      sequence_1 <- DRBX_1[DRBX_1$allele == allele_1, ]
    }
  }
  if (!is.na(allele_2)) {
    if (nrow(DRBX_2[DRBX_2$allele == allele_2, ]) < 1) {
      sequence_2 <- DRBX_2[grepl(paste0(allele_2, ":01"), DRBX_2$allele), ]
      sequence_2 <- sequence_2[1, ]
    } else {
      sequence_2 <- DRBX_2[DRBX_2$allele == allele_2, ]
    }
  }
  
  if(is.na(allele_1)){
    sequence_1 <- zeros
  }
  if(is.na(allele_2)){
    sequence_2 <- zeros
  }
  
  
  study_temp <- data.frame(matrix(nrow = 1, ncol = 0))
  for (j in 1:ncol(sequence_1)) {
    study_temp[paste0(colnames(sequence_1)[j], "_1")] <- sequence_1[j]
    study_temp[paste0(colnames(sequence_2)[j], "_2")] <- sequence_2[j]
  }
  
  study_temp <- cbind(name, study_temp)
  return (study_temp)
  
}

unique_bigdawg <- function(bigdawg_file){
  
  unique_col <- c(1,2)
  for (i in 3:ncol(bigdawg_file)) {
    if(length(unique(bigdawg_file[,i])) > 1) {
      unique_col <- append(unique_col, i, after = length(unique_col))
    }
  }
  unique_col
  
  for (val in unique_col){
    if (val %% 2 != 0){
      if (!((val + 1) %in% unique_col)) {
        pos <- which(unique_col == val)
        
        if (length(pos) > 0) {
          insert_pos <- pos - 1  # insert before
          unique_col <- append(unique_col, val + 1, after = insert_pos)
        }
      }
    }
    if (val %% 2 == 0){
      if (!((val - 1) %in% unique_col)) {
        pos <- which(unique_col == val)
        
        if (length(pos) > 0) {
          insert_pos <- pos - 1  # insert before
          unique_col <- append(unique_col, val - 1, after = insert_pos)
        }
      }
    }
  }
  return (bigdawg_file[, unique_col ])
}

codon_to_amino_acid <- function(codon) {
  codon <- toupper(codon)
  codon_table <- list(
    "TTT" = "Phe", "TTC" = "Phe",
    "TTA" = "Leu", "TTG" = "Leu",
    "CTT" = "Leu", "CTC" = "Leu", "CTA" = "Leu", "CTG" = "Leu",
    "ATT" = "Ile", "ATC" = "Ile", "ATA" = "Ile",
    "ATG" = "Met",
    "GTT" = "Val", "GTC" = "Val", "GTA" = "Val", "GTG" = "Val",
    "TCT" = "Ser", "TCC" = "Ser", "TCA" = "Ser", "TCG" = "Ser",
    "CCT" = "Pro", "CCC" = "Pro", "CCA" = "Pro", "CCG" = "Pro",
    "ACT" = "Thr", "ACC" = "Thr", "ACA" = "Thr", "ACG" = "Thr",
    "GCT" = "Ala", "GCC" = "Ala", "GCA" = "Ala", "GCG" = "Ala",
    "TAT" = "Tyr", "TAC" = "Tyr",
    "TAA" = "Stop", "TAG" = "Stop",
    "CAT" = "His", "CAC" = "His",
    "CAA" = "Gln", "CAG" = "Gln",
    "AAT" = "Asn", "AAC" = "Asn",
    "AAA" = "Lys", "AAG" = "Lys",
    "GAT" = "Asp", "GAC" = "Asp",
    "GAA" = "Glu", "GAG" = "Glu",
    "TGT" = "Cys", "TGC" = "Cys",
    "TGA" = "Stop",
    "TGG" = "Trp",
    "CGT" = "Arg", "CGC" = "Arg", "CGA" = "Arg", "CGG" = "Arg",
    "AGT" = "Ser", "AGC" = "Ser",
    "AGA" = "Arg", "AGG" = "Arg",
    "GGT" = "Gly", "GGC" = "Gly", "GGA" = "Gly", "GGG" = "Gly"
  )
  
  if (codon %in% names(codon_table)) {
    return(codon_table[[codon]])
  } else {
    return(NA)
  }
}

get_sig <- function(result) {
  OR_SB <- as.data.frame(result$L$Set1$OR)
  sig_SB <- OR_SB[OR_SB$sig == "*", ]
}

get_sig_a <- function(result){
  OR_SB <- as.data.frame(result$A$Set1$OR)
  sig_SB <- OR_SB[OR_SB$sig == "*", ]
}



same_aminos <- function(sig) {
  sig_codons <- sig[, 1:2]
  positions <- c()
  sig_codons$Locus <- as.character(sig_codons$Locus)
  sig_codons$Allele <- as.character(sig_codons$Allele)
  for(i in 1:nrow(sig_codons)){
    if (nrow(sig_codons[sig_codons$Locus == sig_codons$Locus[i], ]) != 2){
      next
    } else {
      same_position <- sig_codons[sig_codons$Locus == sig_codons$Locus[i], ]
      if (!is.na(codon_to_amino_acid(same_position$Allele[[1]])) && !is.na(codon_to_amino_acid(same_position$Allele[[2]]))){
        if (codon_to_amino_acid(same_position$Allele[[1]]) == codon_to_amino_acid(same_position$Allele[[2]])) {
          positions <- c(positions, same_position$Locus[[1]])
        }
      }
    }
  }
  positions <- unique(positions)
  return(positions)
}

process_locus <- function(data_list, locus_name, bigdawg_ana) {
  locus <- list()
  for (i in seq_along(data_list)) {
    set <- data_list[[i]]
    name <- names(data_list)[i]
    raw <- BIGDAWG(Data = set, HLA = FALSE, Return = TRUE, Run.Tests = "L")
    sig <- get_sig(raw)
    same_aminos <- same_aminos(sig)
    current_data <- list(
      "Raw Analysis" = raw,
      "Significant positions" = sig,
      "Same Amino + Sig Positions" = same_aminos
    )
    locus[[name]] <- current_data
  }
  aminos_overlap_temp <- list()
  for (i in seq_along(data_list)) {
    val <- locus[[i]]
    name <- names(locus)[i]
    aminos_overlap_temp[[name]] <- val$`Same Amino + Sig Positions`
  }
  aminos_overlap <- Reduce(intersect, aminos_overlap_temp)
  overlap <- list(
    "Overlapping Aminos + Sig" = aminos_overlap
  )
  locus[["Overlapping"]] <- overlap
  bigdawg_ana[[locus_name]] <- locus
  return(bigdawg_ana)
}

bigdawg_analysis <- function(A = NA, B = NA, C = NA, DPA1 = NA, DPB1 = NA, DQA1 = NA, DQB1 = NA, DRB1 = NA) {
  bigdawg_ana <- list()
  
  if (!missing(A)) {
    if (class(A) == "list" && length(A) > 1) {
      bigdawg_ana <- process_locus(A, "A", bigdawg_ana)
    } else {
      stop("Input A must be a list and have a minimum length of 2")
    }
  }
  
  if (!missing(B)) {
    if (class(B) == "list" && length(B) > 1) {
      bigdawg_ana <- process_locus(B, "B", bigdawg_ana)
    } else {
      stop("Input B must be a list and have a minimum length of 2")
    }
  }
  
  if (!missing(C)) {
    if (class(C) == "list" && length(C) > 1) {
      bigdawg_ana <- process_locus(C, "C", bigdawg_ana)
    } else {
      stop("Input C must be a list and have a minimum length of 2")
    }
  }
  
  if (!missing(DPA1)) {
    if (class(DPA1) == "list" && length(DPA1) > 1) {
      bigdawg_ana <- process_locus(DPA1, "DPA1", bigdawg_ana)
    } else {
      stop("Input DPA1 must be a list and have a minimum length of 2")
    }
  }
  
  if (!missing(DPB1)) {
    if (class(DPB1) == "list" && length(DPB1) > 1) {
      bigdawg_ana <- process_locus(DPB1, "DPB1", bigdawg_ana)
    } else {
      stop("Input DPB1 must be a list and have a minimum length of 2")
    }
  }
  
  if (!missing(DQA1)) {
    if (class(DQA1) == "list" && length(DQA1) > 1) {
      bigdawg_ana <- process_locus(DQA1, "DQA1", bigdawg_ana)
    } else {
      stop("Input DQA1 must be a list and have a minimum length of 2")
    }
  }
  
  if (!missing(DQB1)) {
    if (class(DQB1) == "list" && length(DQB1) > 1) {
      bigdawg_ana <- process_locus(DQB1, "DQB1", bigdawg_ana)
    } else {
      stop("Input DQB1 must be a list and have a minimum length of 2")
    }
  }
  
  if (!missing(DRB1)) {
    if (class(DRB1) == "list" && length(DRB1) > 1) {
      bigdawg_ana <- process_locus(DRB1, "DRB1", bigdawg_ana)
    } else {
      stop("Input DRB1 must be a list and have a minimum length of 2")
    }
  }
  
  return(bigdawg_ana)
}


phenotype_bigdawg <- function(Rawdata, col1, col2, name1, name2, phenotype = TRUE) {
  #takes in a dataset of many subjects and their alleles
  if (phenotype == TRUE){
    subject_ID <- Rawdata[1:nrow(Rawdata), 1]
    phenotype <- data.frame()
    for (val in subject_ID) {
      phen_temp <- data.frame(matrix(nrow = 1, ncol = 0))
      # Check if row is a case or a control
      if (substring(val, 3, 3) == "P") {
        phen_temp[["phenotype"]] <- 0
      } else {
        phen_temp[["phenotype"]] <- 1
      }
      phenotype <- rbind(phenotype, phen_temp)
    }
    
    
    
    dataPhenotype <- cbind(Rawdata[, 1], phenotype, Rawdata[, col1:col2])
    colnames(dataPhenotype) <- c("Name", "Phenotype", name1, name2)
    dataPhenotype <- dataPhenotype[2:nrow(dataPhenotype), ]
    
    print(dataPhenotype)
    
  } else {
    dataPhenotype <- cbind(Rawdata[, 1:2], Rawdata[, col1:col2])
    colnames(dataPhenotype) <- c("Name", "Phenotype", name1, name2)
    dataPhenotype <- dataPhenotype[1:nrow(dataPhenotype), ]
  }
  return(dataPhenotype)
}



same_amino_positions_in_data <- function(rawdata, amino_positions){

  new_data <- data.frame(matrix(nrow = nrow(rawdata)))
  for (pos in amino_positions){
    new_data <- cbind(new_data, rawdata[, startsWith(names(rawdata), paste0(pos, "_"))])
  }
  new_data <- cbind(rawdata[, 1:4], new_data[, 2:ncol(new_data)])
  return(new_data)
}



same_aminos_OR <- function(OR, positions) {
  new_OR <- data.frame(matrix(ncol = 7))
  colnames(new_OR) <- colnames(OR)
  for (pos in positions) {
    new_OR <- rbind(new_OR, OR[OR$Locus == pos, ])
  }
  return (new_OR[2:nrow(new_OR), ])
}

codon_and_allele_OR <- function(same_amino_data, sig_allele_OR, sig_codon_OR) {
  allele_data <- same_amino_data[, 3:ncol(same_amino_data)]
  final_df <- data.frame(matrix(ncol = ncol(allele_data)*2), check.names = FALSE)
  for(i in 1:nrow(allele_data)) {
    row <- allele_data[i, ]
    temp_row <- list()
    evens <- TRUE
    odds <- TRUE
    for (j in 1:ncol(row)) {
      colname <- colnames(row)[j]
      pos_name <- strsplit(colname, "_")[[1]][1]
      val <- row[[j]]
      temp_row[[colname]] <- val
      if (j %% 2 == 0) {
        if (evens == FALSE) {
          temp_row[[paste0(colname, "_OR")]] <- NA
          next
        }
      } else {
        if (odds == FALSE) {
          temp_row[[paste0(colname, "_OR")]] <- NA
          next
        }
      }
      if (j == 1 | j == 2) {
        if (length(strsplit(val, ":")[[1]]) > 3) {
          positions <- gregexpr(":", val)[[1]]
          val <- substring(val, 1, positions[3] - 1)
        }
        if (length(sig_allele_OR$OR[grepl(val, sig_allele_OR$Allele)]) != 0) {
          temp_row[[paste0(colname, "_OR")]] <- as.numeric(sig_allele_OR$OR[grepl(val, sig_allele_OR$Allele)][[1]])
        } else {
          temp_row[[paste0(colname, "_OR")]] <- NA
          if (j == 1) {
            odds <- FALSE
          } else {
            evens <- FALSE
          }
        }
      } else {
        if (length(sig_codon_OR$OR[sig_codon_OR$Allele == val & sig_codon_OR$Locus == pos_name]) != 0) {
          temp_row[[paste0(colname, "_OR")]] <- as.numeric(sig_codon_OR$OR[sig_codon_OR$Allele == val & sig_codon_OR$Locus == pos_name][[1]])
        } else {
          temp_row[[paste0(colname, "_OR")]] <- NA
        }
      }
    }
    for (i in seq_along(temp_row)) {
      if (length(temp_row[[i]]) == 0) {
        temp_row[[i]] <- NA
      }
    }
    temp_row <- as.data.frame(temp_row, check.names = FALSE)

    colnames(final_df) <- colnames(temp_row)

    final_df <- rbind(final_df, temp_row)
  }
  final_df <- final_df[2:nrow(final_df),]
  return(final_df)
}

codon_and_allele_count <- function(same_amino_data, sig_allele_OR, sig_codon_OR, sig_allele_count, sig_codon_count) {
  finaldf <- suppressWarnings({
    allele_data <- same_amino_data[, 3:ncol(same_amino_data)]
    final_df <- data.frame(matrix(ncol = ncol(allele_data)*3), check.names = FALSE)
    sig_allele_count <- as.data.frame(sig_allele_count)
    sig_codon_count <- as.data.frame(sig_codon_count)
    for(i in 1:nrow(allele_data)) {
      row <- allele_data[i, ]
      temp_row <- list()
      evens <- TRUE
      odds <- TRUE
      for (j in 1:ncol(row)) {
        
        colname <- colnames(row)[j]
        val <- row[[j]]
        temp_row[[colname]] <- val
        if (j %% 2 == 0) {
          if (evens == FALSE) {
            temp_row[[paste0(colname, "_group0")]] <- NA
            temp_row[[paste0(colname, "_group1")]] <- NA
            next
          }
        } else {
          if (odds == FALSE) {
            temp_row[[paste0(colname, "_group0")]] <- NA
            temp_row[[paste0(colname, "_group1")]] <- NA
            next
          }
        }
        if (j == 1 | j == 2) {
          if (length(strsplit(val, ":")[[1]]) > 3) {
            positions <- gregexpr(":", val)[[1]]
            val <- substring(val, 1, positions[3] - 1)
          }
          if (length(sig_allele_OR$OR[grepl(val, sig_allele_OR$Allele)]) != 0) {
            temp_row[[paste0(colname, "_group0")]] <- as.numeric(sig_allele_count$Group.0[grepl(val, sig_allele_count$Allele)][[1]])
            temp_row[[paste0(colname, "_group1")]] <- as.numeric(sig_allele_count$Group.1[grepl(val, sig_allele_count$Allele)][[1]])
          } else {
            temp_row[[paste0(colname, "_group0")]] <- NA
            temp_row[[paste0(colname, "_group1")]] <- NA
            if (j == 1) {
              odds <- FALSE
            } else {
              evens <- FALSE
            }
          }
        } else {
          if (length(sig_codon_OR$OR[grepl(val, sig_codon_OR$Allele)]) != 0) {
            temp_row[[paste0(colname, "_group0")]] <- as.numeric(sig_codon_count$Group.0[grepl(val, sig_codon_count$Allele)][[1]])
            temp_row[[paste0(colname, "_group1")]] <- as.numeric(sig_codon_count$Group.1[grepl(val, sig_codon_count$Allele)][[1]])
          } else {
            temp_row[[paste0(colname, "_group0")]] <- NA
            temp_row[[paste0(colname, "_group1")]] <- NA
          }
        }
      }
      temp_row <- as.data.frame(temp_row)
      colnames(final_df) <- colnames(temp_row)
      final_df <- rbind(final_df, temp_row)
    }
    final_df <- final_df[2:nrow(final_df), ]
  })
  return(final_df)
}

sig_codon_plus1 <- function (sig_amino, sig_codon) {
  overlapping_aminos <- data.frame(matrix(ncol = 7))
  colnames(overlapping_aminos) <- colnames (sig_codon)
  sig_amino <- rbind(sig_SDQA1_A[grepl(position_ap, sig_SDQA1_A$Position), 1:7], sig_SDQA1_A[grepl(position_am, sig_SDQA1_A$Position), 1:7])
}

three_field_data <- function (rawdata) {

  for (i in 1:nrow(rawdata)){
    for (j in 3:ncol(rawdata)){
      current_val <- rawdata[i, j]
      current_val <- as.character(current_val)
      if (!is.na(current_val) && length(strsplit(current_val, ":")[[1]]) > 3){
        parts <- strsplit(current_val, ":")[[1]]
        trimmed <- paste(parts[1:3], collapse = ":")
        rawdata[i,j] <- trimmed
      }
    }
  }
  return (rawdata)
}

three_field_data_col <- function (rawdata) {
  for (i in 1:nrow(rawdata)){
    for (j in 2){
      current_val <- rawdata[i, j]
      current_val <- as.character(current_val)
      if (!is.na(current_val) && length(strsplit(current_val, ":")[[1]]) > 3){
        parts <- strsplit(current_val, ":")[[1]]
        trimmed <- paste(parts[1:3], collapse = ":")
        rawdata[i,j] <- trimmed
      }
    }
  }
  return (rawdata)
}

average_dis_OR <- function(orTable) {
  complete_average <- 0.0
  complete_count <- 0
  for (i in 1:nrow(orTable)) {
    current_person <- orTable[i, ]
    allele_1_OR <- orTable[i, 2]
    allele_2_OR <- orTable[i, 4]
    if (is.na(allele_1_OR) && is.na(allele_2_OR)){
      next
    } 
    if (!is.na(allele_1_OR)) {
      codon_1_colnames <- colnames(orTable)[5:length(colnames(orTable))]
      codon_1_OR <- orTable[, 5:length(orTable)][i, grepl("_1_", codon_1_colnames)]
      codon_sum_1 <- sum(as.numeric(codon_1_OR), na.rm = TRUE)
      codon_count_1 <- sum(!is.na(codon_1_OR))
      codon_average_1 <- codon_sum_1 / codon_count_1
      average_diff_1 <- abs(allele_1_OR - codon_average_1)
      if (complete_count == 0) {
        complete_average <- average_diff_1
        complete_count <- complete_count + 1
      } else {
        complete_average <- (complete_average * complete_count + average_diff_1) / (complete_count + 1)
        complete_count <- complete_count + 1
      }
    } 
    if (!is.na(allele_2_OR)) {
      codon_2_colnames <- colnames(orTable)[5:length(colnames(orTable))]
      codon_2_OR <- orTable[, 5:length(orTable)][i, grepl("_2_", codon_2_colnames)]
      codon_sum_2 <- sum(as.numeric(codon_2_OR), na.rm = TRUE)
      codon_count_2 <- sum(!is.na(codon_2_OR))
      codon_average_2 <- codon_sum_2 / codon_count_2
      average_diff_2 <- abs(allele_2_OR - codon_average_2)
      if (complete_count == 0) {
        complete_average <- average_diff_2
        complete_count <- complete_count + 1
      } else {
        complete_average <- (complete_average * complete_count + average_diff_2) / (complete_count + 1)
        complete_count <- complete_count + 1
      }
    }
  }
  return(complete_average)
}


sig_allele_amino_pos <- function(sig_allele, sig_codon_OR, amino_pos, loci_codon) {
  
  all_alleles <- data.frame()
  
  for (i in 1:nrow(sig_allele)) {
    
    temp_row <- data.frame(matrix(nrow = 1))
    row <- sig_allele[i, ]
    
    allele <- row[[2]]
    OR <- row[[3]]
    
    temp_row[["Allele"]] <- allele
    temp_row[["Allele_OR"]] <- OR
    
    subset_rows <- loci_codon[grepl(allele, loci_codon$allele), amino_pos]
    
    if (is.null(nrow(subset_rows)) || nrow(subset_rows) == 0 || ncol(subset_rows) == 0) {
      next
    }
    
    raw_sequence <- subset_rows[1, ]
    
    for (j in 1:ncol(raw_sequence)) {
      
      pos_name <- colnames(raw_sequence)[j]
      pos <- raw_sequence[[j]]
      
      match_vals <- sig_codon_OR$OR[
        sig_codon_OR$Allele == pos & sig_codon_OR$Locus == pos_name
      ]
      
      if (length(match_vals) > 0) {
        pos_OR <- as.numeric(match_vals[1])
      } else {
        pos_OR <- NA
      }
      
      temp_row[[pos_name]] <- pos
      temp_row[[paste0(pos_name, "_OR")]] <- pos_OR
    }
    
    all_alleles <- rbind(all_alleles, temp_row)
  }
  
  if (nrow(all_alleles) == 0) return(NULL)
  return(all_alleles[, 2:ncol(all_alleles)])
}



sig_allele_amino_pos_binned <- function (sig_allele, sig_codon_OR, amino_pos, loci_codon) {
  all_alleles <- data.frame()
  for (i in 1:nrow(sig_allele)) {
    temp_row <- data.frame(matrix(nrow = 1))
    row <- sig_allele[i, ] 
    allele <- row[[2]]
    temp_row[["Allele"]] <- allele
    temp_row[["Allele_Group.0"]] <- row[[3]]
    temp_row[["Allele_Group.1"]] <- row[[4]]
    raw_sequence <- loci_codon[grepl(allele, loci_codon$allele), amino_pos][1, ]
    for (j in 1:ncol(raw_sequence)) {
      pos_name <- colnames(raw_sequence)[[j]]
      pos <- raw_sequence[[j]]
      print(pos)
      if (length(sig_codon_OR$OR[sig_codon_OR$Allele == pos]) != 0) {
        pos_OR <- as.numeric(sig_codon_OR$OR[sig_codon_OR$Allele == pos][[1]])
      } else pos_OR <- NA
      temp_row[[pos_name]] <- pos
      temp_row[[paste0(pos_name, "_OR")]] <- pos_OR
    }
    
    all_alleles <- rbind(all_alleles, temp_row)
  }
  return(all_alleles[, 2:ncol(all_alleles)])
}



