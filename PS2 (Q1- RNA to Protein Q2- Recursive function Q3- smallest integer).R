## Q1 - 6pt
## Every amino acid in a protein is encoded by three nucleotides. 
## Execute the following two lines to get a list of all codons and corresponding amino acids
codons = c('UUU','UUC','UUA','UUG','UCU','UCC','UCA','UCG','UAU','UAC','UAA','UAG','UGU','UGC','UGA','UGG','CUU','CUC','CUA','CUG','CCU','CCC','CCA','CCG','CAU','CAC','CAA','CAG','CGU','CGC','CGA','CGG','AUU','AUC','AUA','AUG','ACU','ACC','ACA','ACG','AAU','AAC','AAA','AAG','AGU','AGC','AGA','AGG','GUU','GUC','GUA','GUG','GCU','GCC','GCA','GCG','GAU','GAC','GAA','GAG','GGU','GGC','GGA','GGG')
amino_acids = c('F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G' )

## Write a function that will take a coding region sequence as input. You can assume the sequence is starting with AUG and is a multiple of three nucleotides. 
## The output should be the corresponding protein sequence. Report up to the first stop codon. 


# Example input/output:
# translate_rna_to_protein("AUGCUGGUGUAGUCGUGGUUAUUCUUU")
# [1] "MLV"

# "MLV*SWLFF" will be acceptable but
#  the protein sequence should ideally end at stop codons.

translate_rna_to_protein = function ( rna_seq) { 
  protein_seq = ""
  separated_rna <- gsub("(.{3})", "\\1 ", rna_seq)
  separated_rna = unlist(strsplit(separated_rna, split=" "))
  
  for (codon in separated_rna){
    amino_acid <- which(codons==codon)
    if (amino_acids[amino_acid]=="*"){
      break
    }
    else{
      protein_seq = paste(protein_seq, amino_acids[amino_acid], sep="")
    }
    }
  
  
  return (protein_seq)
}

## Q2 - 6pt

## Given a positive integer n,

## write a recursive function that calculates of squares of consecutive integers up to n

## series_sum(4) = 4^2 + 3^2 + 2^2 + 1^1 = 30

# Example input/output:
# series_sum(20)
# [1] 2870

series_sum <- function(n){
  if (n==1){
    return(1)
  }
  else{
    return((n^2) + series_sum((n-1)))
  }
  
  
}

## Q3 - 9pt
## Let's define
## F(n) = F(n-1) + 3* F(n-2)
## F(1) = 1
## F(0) = 0
## Write a function that takes any positive integer k as input and 
## prints the value of smallest "n" such that F(n) >= k

# Example input/output:
# smallest_n_finder(20)
# [1] 6

smallest_n_finder = function ( k) { 
  n=0
  recursive_function <- function(n){
    if (n==1){
      return(1)
    }
    if (n==0){
      return(0)
    }
    return(recursive_function(n-1)+3*recursive_function(n-2))
  }
  
  while (k>0){
    fn <- recursive_function(n)
    if (fn<k){
      n=n+1
    }
    if (fn>=k){
      return(n)
      break
    }
  }
  
}


