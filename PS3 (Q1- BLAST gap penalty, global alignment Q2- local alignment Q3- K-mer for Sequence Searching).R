## BCH 339N Problem Set
## Sequence Alignment including BLAST. 

## Q1 - 12 pt
## In class, we have discussed dynamic programming for global alignment of two sequences. 
## Please implement this algorithm. 
## Specifically, the function should take two strings in additon to scores for match, mismatch, gap penalty
## Note that all mismatches will have the same score for this question. 
## We will use a linear gap penalty (no separate penalty for gap open). 
## The output should be the optimal alignment and the associated score

## Expected Output Example: 
## ATT-AGC
## ATTCAGG
## Score: 18 

# alignment_function("ATTAGC", "ATTCAGG", 1, -1, -2)

alignment_function = function (str1, str2, match_score, mismatch_score, gap_penalty) {
  sequence_1 <- unlist(strsplit(str1, ""))
  sequence_2 <- unlist(strsplit(str2, ""))
  
# CREATE MATRIX
  i_axis <- seq(from=0, to = gap_penalty*length(sequence_1), by = gap_penalty)
  j_axis <- seq(from=0, to = gap_penalty*length(sequence_2), by = gap_penalty)
  initial_matrix <- matrix(0, nrow=(length(j_axis)), ncol=(length(i_axis)))
  initial_matrix[1, ] <- i_axis
  initial_matrix[, 1] <- j_axis
  
# CREATE TRACEBACK MATRIX
  traceback_matrix <- matrix(0, nrow=(length(j_axis)), ncol=(length(i_axis)))
  traceback_matrix[1, ]<- "left"
  traceback_matrix[, 1]<- "up"
  traceback_matrix[1,1] <- "done"

# FIND MAX VALUE FOR MATRIX AND ADDING DIRECTION TO TRACEBACK MATRIX
  find_max_in_matrix <- function(j, i){
    if (sequence_1[i-1]==sequence_2[j-1]){
      diagonal <- (initial_matrix[j-1, i-1] + match_score)
    }
    else{
      diagonal <- (initial_matrix[j-1,i-1] + mismatch_score)
    }

    left <- initial_matrix[j, i-1] + gap_penalty
    up <- initial_matrix[j-1, i] + gap_penalty
    
    max_value <- max(diagonal, up, left)
    
    if (max_value == diagonal){
      traceback_matrix[j, i] <<- "diagonal"
    }
    if (max_value == up){
      traceback_matrix[j, i] <<- "up"
    }
    if (max_value == left){
      traceback_matrix[j, i] <<- "left"
    }

    return ( max_value)

  }

# FILLING IN MATRIX
  for (j in 2:length(j_axis)){

    for (i in 2:length(i_axis)){
      initial_matrix[j, i] <- find_max_in_matrix(j, i)

    }
  }
 
# RETURN SEQUENCES
  new_sequence_1 = c()
  new_sequence_2 = c()
  j=length(j_axis)
  i= length(i_axis)
  score <- initial_matrix[j, i]
  
  while (j>=1 || i>=1){
    if (traceback_matrix[j,i]=="diagonal"){
      new_sequence_1 <- c(sequence_1[i-1], new_sequence_1)
      new_sequence_2 <- c(sequence_2[j-1], new_sequence_2)
      j = j-1
      i = i-1
    }
    
    if (traceback_matrix[j,i]=="up"){
      new_sequence_1 <- c("-", new_sequence_1)
      new_sequence_2 <- c(sequence_2[j-1], new_sequence_2)
      j = j-1
    }
    
    if (traceback_matrix[j,i]=="left"){
      new_sequence_1 <- c(sequence_1[i-1], new_sequence_1)
      new_sequence_2 <- c("-", new_sequence_2)
      i = i-1
    }
    
    if (traceback_matrix[j,i]=="done"){
      break
    }
    }
    
  new_sequence_1 <- paste(new_sequence_1, collapse = "")
  new_sequence_2 <- paste(new_sequence_2, collapse = "")
  print(new_sequence_1)
  print(new_sequence_2)
  print(paste("Score:", score))

}



## Q2 - 3pt

## The above algorithm implements a global alignment. Which line in your code would you need to modify if you wanted to implement local alignment?

# I would add an if statement to the function returning max value. If the max value is less than 0, 
# the function would return 0 instead of a negative number. Additionally, I'd set the first row
# and first column equal to 0 instead of a multiple of the gap penalty. 


## Q3 -6pt 
## Enumerate k-mers
## BLAST uses heuristics to speed up the task of searching a query sequence against a large database. 
## For example, given a word size (default 3), BLAST will create a table of all possible short words (k-mers) contained in the query. 
## Write a function that will create this table of all possible words of given size. 
## For example, given  a a word size 3 and query sequence (LRITSLRI): 
## we will have LRI, RIT, ITS, TSL, SLR, LRI => Hence 5 distinct words are possible. 

test_sequence = "LRITSLRI"
test_sequence2 = "LRITSLRIK"

enumerate_words = function(query_seq, word_size = 3) { 
  start=1
  count=0
  distinct_words <- c()
  while (start<= (nchar(query_seq)-word_size+1)){
    substring <- substr(query_seq, start, start+word_size-1)
    start = start+1
    if (!(substring %in% distinct_words)){
      distinct_words <- c(distinct_words, substring)
      count = count+1
    }
  }
  print(distinct_words)
  print(count)
}

enumerate_words (test_sequence) 
## Expected Output
## [1] 5
enumerate_words (test_sequence2) 
## Expected Output
## [1] 6