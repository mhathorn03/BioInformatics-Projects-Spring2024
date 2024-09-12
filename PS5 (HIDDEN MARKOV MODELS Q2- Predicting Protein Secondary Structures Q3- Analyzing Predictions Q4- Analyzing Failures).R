## BCH 339N
## PS5 - Hidden Markov Models

## We have covered the Viterbi and Forward Algorithms in class. 
## These methods are implemented in the R package - HMM

install.packages("HMM")
library(HMM)
?initHMM
?forward
?viterbi

## In this homework, we are going to develop an HMM that can predict secondary structure of proteins.
## Kaggle (https://www.kaggle.com) organizes data-science and machine learning challenges for the community
## They have many challenges related to bioinformatics. 
## We will explore the https://www.kaggle.com/alfrandom/protein-secondary-structure/data challenge. 
## Read the associated introduction to familiarize yourself with this dataset. 

## Q1 -6pt
## What is protein secondary structure? 
## We are going to work with the three state secondary structure encoding (sst3). 
## Explain the biochemical features of the three states that we will model (C, E, and H ).

# The protein secondary structure is the local configuration formed from the peptide backbone.
# Interactions between amino acid residues can cause two different secondary structure configurations, 
# alpha sheets, and beta sheets. 
# C represents coils or regions within the protein that don't have a defined secondary
# structures. The E state represents beta strands, which are zig-zag-like structures within the protein.
# They're formed by H-bonds between adjacent strands. The H state represents alpha helices, which are
# tightly coiled structures formed by H-bonds between every couple Amino acids.


## Q2 -6pt

## We will use the data with the strict quality controls for this homework. Here Download Hereis the file for convenience. 

## You will have to change the next line of code to specify the folder where the file is. 
## First, we will remove all examples with non-standard amino acids 

prot_sec_data = read.csv('./2018-06-06-pdb-intersect-pisces.csv', stringsAsFactors = F)
str(prot_sec_data)
prot_sec_data_std = prot_sec_data[as.logical(prot_sec_data$has_nonstd_aa) == F, ]



## Next, we are going to split our dataset into two parts. 
## The first part should contain 10% of the sequences and will be used for training our HMM
## The remaining 90% will be used as a test set to assess how well our HMM does in predicting secondary structure. 
## To ensure reproducibility of your code, we are going to use the set.seed function. 
?set.seed
set.seed(3)

train  = sample(nrow(prot_sec_data_std), size = floor(nrow(prot_sec_data_std)*.1 ) )
test =  setdiff(1:nrow(prot_sec_data_std), train  )

train_data = prot_sec_data_std[train, ]
test_data = prot_sec_data_std[test, ]

##  We will use the training data to infer the parameters of the HMM model. 
## Specifically, use sst3 and seq columns of the data and determine the 
## transition, emission and initial probabilities.
## Write a function that will take the training data as input. 
## The output should be a list with names
## c(“initial_probs”, “emission_probs”, “transition_probs”)
## Note that emission_probs and transition_probs should be matrices of the form defined in initHMM. 
symbols = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") 
states =  c("C", "E", "H")

get_hmm_probs<- function(train_data, symbols, states){
  # INIT COUNTS
  c_count = 0
  e_count = 0
  h_count = 0
  ttl_count = 0
  
  for (state in unique(train_data$sst3)) {
    if (substr(state, 1, 1) == states[1]){
      c_count = c_count+1
    }
    if (substr(state, 1, 1) == states[2]){
      e_count = e_count+1
    }
    if (substr(state, 1, 1) == states[3]){
      h_count = h_count+1
    }
    ttl_count=ttl_count+1
  }
  init_counts <- c(c_count/ttl_count, e_count/ttl_count, h_count/ttl_count)
  
  
  # EMISSION PROBS
  emissionProb_mat <- matrix(0, nrow=length(states), ncol=length(symbols))
  ttl_count=0
  c_count = 0
  e_count = 0
  h_count = 0
  
  for (protein in unique(train_data$seq)){
    count=0
    amino_acids <- unlist(strsplit(protein, ""))
    current_state <- unlist(strsplit(train_data$sst3[train_data$seq == protein], ""))
    for (aa in amino_acids){
      count = count+1
      ttl_count = ttl_count+1
      col_index <- which(symbols==aa)
      
      row_index <- which(states==current_state[count])
      
      
      if (current_state[count]=="C"){
        c_count = c_count+1
      }
      if (current_state[count]=="E"){
        e_count = e_count+1
      }
      if (current_state[count]=="H"){
        h_count = h_count+1
      }
      
     
      emissionProb_mat[row_index, col_index] <- emissionProb_mat[row_index, col_index] + 1
    }
    
  }
  emissionProb_mat[1, ] <- emissionProb_mat[1, ] / c_count
  emissionProb_mat[2, ] <- emissionProb_mat[2, ] / e_count
  emissionProb_mat[3, ] <- emissionProb_mat[3, ] / h_count
  
# TRANSITION MATRIX
  
  # Need to use indexes: the transition matrix will be C:C over ttl, C:E over ttl, etc
  c_count <- 0
  e_count <- 0
  h_count <- 0
  transProb_mat <- matrix(0, nrow = length(states), ncol = length(states))
  
  for (protein in (train_data$sst3)) {
    transition_states <- unlist(strsplit(protein, ""))
    for (index in 1:(length(transition_states) - 1)) {
      state <- transition_states[index]
      next_state <- transition_states[index + 1]
      
      col_index <- which(states == state)
      row_index <- which(states == next_state)
      
      transProb_mat[row_index, col_index] <- transProb_mat[row_index, col_index] + 1
      
      if (state == "C") {
        c_count <- c_count + 1
      }
      if (state == "E") {
        e_count <- e_count + 1
      }
      if (state == "H") {
        h_count <- h_count + 1
      } 
    }

  }
  transProb_mat[1, ] <- transProb_mat[1, ] / c_count
  transProb_mat[2, ] <- transProb_mat[2, ] / e_count
  transProb_mat[3, ] <- transProb_mat[3, ] / h_count
  
  
 
  
  return(list(initial_probs=init_counts, emission_probs=emissionProb_mat, transition_probs=transProb_mat)) 
}

hmm_values =  get_hmm_probs(train_data, symbols, states)
hmm_values

## The expected output should be as follows. 
## Do not worry if yours is slightly different. 
# $initial_probs
# [1] 1 0 0
# 
# $emission_probs
# [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]
# [1,] 0.06585424 0.01165331 0.07881398 0.05734681 0.03039266 0.11675254 0.03626634
# [2,] 0.06147861 0.01684346 0.03523090 0.04828457 0.05855106 0.04896633 0.02452327
# [3,] 0.11237240 0.01019533 0.05272974 0.09006702 0.04173868 0.03504955 0.02049012
# [,8]       [,9]      [,10]      [,11]      [,12]      [,13]      [,14]
# [1,] 0.03188721 0.05509976 0.06162143 0.02146717 0.05971928 0.07833321 0.03431193
# [2,] 0.09332077 0.04435443 0.10398829 0.02221732 0.02648834 0.02263841 0.02845341
# [3,] 0.05793930 0.06470303 0.11973293 0.02671922 0.03325915 0.02469259 0.04735854
# [,15]      [,16]      [,17]      [,18]      [,19]      [,20]
# [1,] 0.04505597 0.07720446 0.05752448 0.04383315 0.01026327 0.02659880
# [2,] 0.04784344 0.05181368 0.06474705 0.12821078 0.01862806 0.05341782
# [3,] 0.05923237 0.04983277 0.04145271 0.05986647 0.01582762 0.03674048
# 
# $transition_probs
# [,1]        [,2]       [,3]
# [1,] 0.80925703 0.110084193 0.08065877
# [2,] 0.20296365 0.783140502 0.01389585
# [3,] 0.09893198 0.004737097 0.89633093


## Q3 - 6pt
## Look at how hmms are defined in the HMM package
## We will Use the inferred parameters from Q3 to define an HMM. 
?initHMM
sec_struct_hmm = initHMM(States = states, Symbols = symbols,
                         startProbs = hmm_values$initial_probs, 
                         transProbs = hmm_values$transition_probs, 
                         emissionProbs = hmm_values$emission_probs)


## Next, We are going to assess the performance on the test data
## For each example in the test data, use the given sequence (seq column) to predict the most likely path of hidden states. 
## You can use the appropriate function in the hmm package for this step. 
## For each example, compare the predicted most likely hidden state path with  the experimentally identified values (sst3 column). 
## Output a vector named percent_correct containing the percentage of amino acids whose secondary structure was correctly predicted for each example in test data.

## YOUR CODE HERE

percent_correct <- numeric(length(test_data$seq))
index = 0
count=0

for (sequence in test_data$seq){
  
  index = index+1
  aa <- unlist(strsplit(sequence, ""))
  predicted_path <- viterbi(sec_struct_hmm, observation = aa)
  
  actual_path <- unlist(strsplit(test_data$sst3[index], ""))
  
  
  num_similarities <- sum(predicted_path == actual_path)
  percent_correct[index] <- num_similarities/length(actual_path)
 
}





## Plot the distribution of percent_correct 
hist(percent_correct, 40)

## Determine the mean, median and 90th percentile of this distribution 
mean(percent_correct)

median(percent_correct)

quantile(percent_correct, .90)

# Write a sentence about your interpretation of these results.

# The histogram had a normal distribution, the mean of the results was 0.4549, the median was: 0.4448, 
# and the quantile was 0.7079 at 90%. The model is right about 45% of the time. 


## Q4 - 3pt 
## Identify the training examples where we succeed less than 1%. Output the pdb_ids of these examples
## Examine the predicted hidden state for these and the actual secondary structure. 
## You will notice that these test examples have a high percentage of "E" in their secondary structure. 

failures <- c()
index = 0

for (value in percent_correct){
  index = index+1
  if (value < 0.01){
    failures <- c(failures, index)
  }
}
print(test_data$pdb_id[failures])



## Explain why you think our HMM doesn't work for this certain class of proteins?
## Hint: Think about the biology of the "E" state and the definition of a Markov chain

# The HMM doesn't work for this class of protein because the Markov property 
# assumes that the probability of transitioning depends only on the current state, 
# and not the sequences preceding it. Since these proteins have a lot of beta sheets 
# (meaning they depend on the previous motifs and patterns to form H-bonds), The HMM 
# isn't able to accurately predict a beta sheet since the beta sheet depends on the
# amino acids multiple nucleotides away. 

