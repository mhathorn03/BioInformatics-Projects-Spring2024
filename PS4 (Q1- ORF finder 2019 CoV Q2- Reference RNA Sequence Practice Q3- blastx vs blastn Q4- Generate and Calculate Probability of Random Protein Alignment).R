## BCH 339N Problem Set
## Sequence Alignment including BLAST. 


## Q1 - 6pt
## Sequences for the 2019-nCoV are available through NCBI: https://www.ncbi.nlm.nih.gov/genbank/2019-ncov-seqs/Links to an external site.
## Wuhan_nCov19.fasta (see Wuhan_nCov19-1.fasta attached Download attached) contains the complete sequence of one isolate of the 2019 Coronavirus
## An open reading frame (ORF) is a stretch of sequence that starts with a start codon (ATG) 
## and ends with a stop codon (TAG, TAA, TGA) without any other stop codons in between. 
## By enumerating ORFs, we can create a list of candidate proteins encoded by this coronavirus. 
## Write a function that will find all ORFs with at least 59 aa. 
## Your function should output start and stop position of all such ORFs. 

ncov = readLines('Wuhan_nCov19.fasta')
sequence = paste(ncov[2:length(ncov)], collapse= "" )

## We can collapse the multiline FASTA as follows:
ncov_orf_finder <- function(query_sequence, orf_length_threshold = 60) {
  start_codon <- "ATG"
  stop_codon <- c("TAG", "TAA", "TGA")
  ORF_list <- c()
 
  
  find_ORF <- function(ORF_start){
    for (i in seq(ORF_start, nchar(query_sequence), by = 3)){
      codon <- substr(query_sequence, i, i+2)
      if (codon %in% stop_codon){
        if ((i+2 - ORF_start) >= (orf_length_threshold*3)){
          ORF_list <<- c(ORF_list, paste(ORF_start,"-",i+2, sep = ""))
          return(i+3)
        }
        else{
          return(i+3)
        }
      }
      else if (i>=(nchar(query_sequence))-2){
        return(nchar(query_sequence))
      }
    }
  }
  
  
  
  nt=1
  while (nt<(nchar(query_sequence)-2)){
    if (substr(query_sequence, nt, nt+2)==start_codon){
      nt <- find_ORF(nt)
    }
    else{
      nt = nt+3
    }
     
  }
  nt=2
  while (nt<(nchar(query_sequence)-2)){
    if (substr(query_sequence, nt, nt+2)==start_codon){
      nt <- find_ORF(nt)
    }
    else{
      nt = nt+3
    }
    
  }
  nt=3
  while (nt<(nchar(query_sequence)-2)){
    if (substr(query_sequence, nt, nt+2)==start_codon){
      nt <- find_ORF(nt)
    }
    else{
      nt = nt+3
    }
    
  }
  for (ORF in ORF_list){
    print(ORF)
  }
  
  }
  
 

ncov_orf_finder(sequence)



## Example Output: 

## Note that the order doesn't need to be identical to the example output.
# [1] "28734-28955"
# [1] "28284-28577"
# [1] "27894-28259"
# [1] "26523-27191"
# [1] "21936-22199"
# [1] "10215-10400"
# [1] "6156-6350"
# [1] "2958-3206"
# [1] "28274-29533"
# [1] "21536-25384"
# [1] "15461-15667"
# [1] "266-13483"
# [1] "27394-27759"
# [1] "27202-27387"
# [1] "26245-26472"
# [1] "25393-26220"
# [1] "13768-21555"



## Q2 - 6pt
## A potential ORF we identified in Q1 spans the nucleotides "26523-27191"
substr(sequence, 26523, 27191)
## This is one of the proteins encoded by the new 2019 coronavirus genome. 
## Use blast to search for similar sequences.

## For database, select: Reference RNA sequences (refseq_rna). 
## You will notice that the first few hits match the Severe acute respiratory syndrome coronavirus 2.
## What other organism(s) did you find? 

## Hint: Use the "exclude" option. 
## Click on at least one such alignment and explain in your own words the output. Include a description of 
## Score, Expect and Identities and include a screenshot. 


# Filtering for highly similar sequences, I found the SARS coronavirus Tor2, complete genome
# which had a score of 706, an E-value of 0, and identity of 85.85%. I also found the Bat coronavirus which had a
# score of 427, E-value of 2e-115, and identity of 78.56%. Both alignments scores were equal to their maximum
# score but the higher score of the corona virus 2 indicates that the alignments result in higher overall
# similarity. The low e-values show that it's very unlikely a similar sequence would be found by random chance, showing
# that these two sequences probably evolved from each other. Finally, the identity shows the percentage of nucleotides
# that were the same between the two sequences. 


## Q3 - 3pt
## We will use the same sequence for this question
substr(sequence, 26523, 27191)
## Use blastx to search the database. 
## Describe how blastx differs from blastn.
## Given the results page, what do you think is the most likely function of this protein? 


# Blastx translates the nucleotide sequence into protein sequences and compares them to the protein sequences 
# in the database. Blastn compares the inputted nucleotide sequence to a nucleotide
# sequence database. These differences are useful depending on what you're searching for (homologous 
# chromosomes, proteins that may have evolved, etc.). Blastx provided many results but they were all variations of
# the covid membrane glycoprotein. Based on this, the most likely function of this protein is related to 
# the membrane glycoprotein and could help with things such as viral entry, membrane fusion, or interaction with 
# the host cell. 


## Q4 - 6pt

## We talked about the intuition behind assessing statistical significance in our discussion of BLAST. 

## In this problem, we will see another example of how to think about statistical significance. 

## One of the potential ORF we identified in Q1 spans the nucleotides "266-13483"
candidate_orf = substr(sequence, 266, 13483)

## Install and load the Biostrings package.

install.packages("Biostrings")
library (Biostrings)
## We can calculate the nucleotide frequencies with the following code
candidate_orf = substr(sequence, 266, 13483)
nuc_freq = letterFrequency( DNAString(candidate_orf), letters = c("A", "C", "G", "T"), as.prob = T )
nuc_freq

## Similarly, we can count trinucleotide frequencies 
tri_freq_orf = trinucleotideFrequency(DNAString(candidate_orf), as.prob =  T )
tri_freq_orf



## For a random sequence,
## we would expect its trinucleotide frequencies to be a simple multiplication of the individual base frequencies. 
## For example, Freq (ATG) = Freq(A) * Freq(T) * Freq(G) 
freq_ATG_expected = nuc_freq[1] * nuc_freq[4] * nuc_freq[3]
## One way to quantify the difference between expected and observed trinucleotide frequencies is to use the sum of squared of the differences. 
( tri_freq_orf['ATG'] - freq_ATG_expected ) ^2

## The following function quantifies this across all codons. 

diff_trinuc_freq = function (nuc_freq, tri_freq) { 
  expected_freqs = numeric(0)
  for (first in 1:4) { 
    for (second in 1:4) { 
      for (third in 1:4) { 
        expected_freqs = c(expected_freqs, nuc_freq[first] * nuc_freq[second] * nuc_freq[third] ) 
      }
    }  
  }
  return ( sum ( (tri_freq - expected_freqs)^2 ) ) 
}

candidate_orf_dif  = diff_trinuc_freq (nuc_freq,  tri_freq = tri_freq_orf)



## Please generate 1000 random sequences of length 13218  using nuc_freq.

## In other words, each sequence will be a string of length 13218 such that the probability that each nucleotide is (A,C,G,T) is equal to nuc_freq.  

random_sequences <- lapply(1:1000, function(nucleotide){
  paste(sample(c("A", "C", "G", "T"), 13218, replace = TRUE, prob = nuc_freq), collapse = "")})


random_sequences <- unlist(random_sequences)
str(random_sequences)  
# chr [1:1000] "GCGTCCTTTCCAGGTGAGGCATTATTCCTAATATATTGGTTATCCCTTTATCACAATAGTGTCTCATGCCAATTTAAGCTTTCGTGACTCTTAGAAGTATTGACTTATCGG"| __truncated__ ...

## Note that the sequences need not have the exact single nucleotide frequencies as they are randomly generated. 

## For example

letterFrequency( DNAString(random_sequences[2]), letters = c("A", "C", "G", "T"), as.prob = T )        

# A         C         G         T 
# 
# 0.2972462 0.1751400 0.2040399 0.3235739 
# 
# nuc_freq        
# 
# A         C         G         T 
# 
# 0.2988349 0.1763504 0.2009381 0.3238765 




## Plot the diff_trinuc_freq of these 1000 sequences and compare it to candidate_orf_dif. 
## Write a paragraph to interpret the results and how this exercise relates to statistical significant calculation in BLAST. 

random_sequences_trifreq = sapply (random_sequences, function(x){ trinucleotideFrequency (DNAString( x ), as.prob = T ) } )

random_diffs = apply (random_sequences_trifreq, 2, function(x) {diff_trinuc_freq(nuc_freq, x) } )

hist( random_diffs, xlim = c(0, 0.002))
abline ( v = candidate_orf_dif, col = "red")



# The histogram results indicate a statistically significant difference in trinucleotide frequencies of the candidate ORF. 
# The candidate ORF has a much higher random difference value than the random trinucleotide frequencies
# (which were based on random chance). This suggests that the candidate ORF is probably influenced by selective biological 
# pressures, making it difficult to encounter randomly in nature. This corresponds to a low E-value in blast. Low E-value 
# indicates high statistical significance and a low probability that alignment between sequences occurred by chance. Similarily, 
# since there is a large discrepancy between the candidate ORF and the randomly generate trinucleotide frequencies (containing the
# same nucleotide probability), both the histogram and low E-value in BLAST show that the candidate ORF pattern is unlikely to occur 
# by chance, suggesting a biological relevance. 