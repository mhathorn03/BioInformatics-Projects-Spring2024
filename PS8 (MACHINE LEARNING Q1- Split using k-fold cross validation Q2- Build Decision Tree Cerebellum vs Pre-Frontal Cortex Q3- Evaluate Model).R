## BCH 339N
## PS8 - Machine Learning


# In this problem set, we will be using machine learning to create a model
# that can predict the specific region of the brain based on its gene expression
# I preprocessed the RNA expression data from GTEX (https://www.gtexportal.org/home/) for this PS. 
# Specifically, I extracted a small number of samples to speed up the analysis.
# I also selected a subset of genes with high variability across samples. 

library(data.table)
library(tree)

# Download the data using fread
gtex_brain_sample = fread ("https://github.com/cancenik/BCH339N/raw/master/Brain_Gtex_sample_wGeneNamesTarget.csv", data.table =  F)
str(gtex_brain_sample)
# Convert the brain regions into a factor
gtex_brain_sample$Target = as.factor(gtex_brain_sample$Target )
table(gtex_brain_sample$Target)

# Note that the organization of the samples is a little different than the previous homework.
# The column names 'V1' holds the sample names starting with GTEX-...
# The gene names begin with "ENSG". These are the other columns. 
# Finally, the last column of this data.frame is called "Target" similar to the example from class.

## Q1- 6pt 
## The second step in our machine learning framework is splitting the data into training and test set
## We learned about the concept of cross-validation in class. 
## Write a function that will generate training-test set split using k-fold cross-validation.
## Specifically, the input is a data.frame or matrix, and an integer k. 
## The output should be a List of length k, where each element is a vector of indices.
## Each vector element of the output list should include the indices of the samples that will be used for training. 
set.seed(3)
k_fold_cv_generator = function (expr_data, k) { 
  n <- nrow(expr_data)
  fold_size <- n %/% k
  k_vectors <- vector("list", k)
  shuffled_indices <- sample(n)
  
  for (i in 1:k) {
    start_index <- (i - 1) * fold_size + 1
    end_index <- (i * fold_size)
    
    test_indices <- start_index:end_index
    training_indices <- shuffled_indices[-test_indices]
    
    k_vectors[[i]] <- training_indices
   
}

  return (k_vectors)
}

kfold_res<- k_fold_cv_generator(gtex_brain_sample, 5)
str(kfold_res)
# List of 5
# $ : num [1:40] 46 45 37 2 29 50 9 34 32 24 ...

## Q2- 9pt
## Build a decision tree to predict whether a given sample from gtex_brain_sample is from the Cerebellum or the Pre-frontal Cortex.
## Be sure to split the data to test and train set(s). 
## Write a paragraph describing the output of your model. 


for (i in 1:length(kfold_res)){
  train_indices <- unlist(kfold_res[[i]])
  train_data <- gtex_brain_sample[train_indices, ]
  decision_tree <- tree(Target ~ ., data = train_data)
  print(decision_tree)
}
 

# Each decision tree split on one gene. The first decision tree predicted the sample to come from the Frontal Cortex if 
# ENSG00000271895.2 < 1321.5 and from the Cerebellum if it was greater than this. 65% of the samples came from cerebellum, and 
# 35% from Frontal cortex. The second decision tree predicted that samples with gene expression ENSG00000092853.13 < 18.5 
# to belong to the cerebellum, and anything greater than this predicted to come from the frontal cortex. The third tree predicted
# the sample to be from the cerebellum if ENSG00000153904.18 < 8125, and from the frontal cortex if  > 8125. The fourth tree
# predicted the Cerebellum if ENSG00000153904.18 < 8977.5, and frontal cortex if > 8977.5. Finally the fifth tree predicted 
# cerebellum if ENSG00000153904.18 < 9601.5, and frontal cortex if > 9601.5. Each tree had about 65% of samples from the cerebellum
# and 35% from the frontal cortex.


## Q3 - 6pt
## Final step is to evaluate our model
## Generate predictions for your test set using the model you trained. 
## How well did you do? 

true_values <- vector("numeric", length(kfold_res))

for (i in 1:length(kfold_res)){
  train_indices <- unlist(kfold_res[[i]])
  train_data <- gtex_brain_sample[train_indices, ]
  test_data <- gtex_brain_sample[-train_indices, ]
  decision_tree <- tree(Target ~ ., data = train_data)
  tree.pred = predict(decision_tree, test_data, type="class")
  true_values[i] = sum(tree.pred==test_data$Target)/nrow(test_data)
}
print(true_values)


# The accuracy of the predictions is 0.8, 0.8, 0.9, 1, 1, for the five respective data trees. 
# I think this model did a good job predicting the location of the sample



