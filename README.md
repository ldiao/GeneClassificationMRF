# GeneClassificationMRF
Discovery of species specific genes given gene expression data from multiple brain regions and species

We provide an implementation of the simulated field algorithm presented in Celeux et al. 2003, used to solve a Markov random field. In our implementation, we take as input a list of prior probabilities associated with multiple brain regions and 3 species. The 3 species result in k = 5 classes:

1. Nonspecific (same expression in all species)
2. Macaque specific
3. Human specific
4. Chimp specific
5. Differentially expressed in all species

Priors can be defined in any way the user wishes, as long as the the input format is correct. The priors, denoted pp, must be a named list of matrices, where each matrix is the priors for a single brain region. Thus, for k = 5 classes, each matrix must have 5 columns and rows equal to teh number of genes.

A toy example of priors can be found in ex.priors.rda. The code includes a working example for how to estimate the MRF parameters, obtain posterior probabilities for class membership, and the posterior best classifications.

That code is:

run.example <- function() {

  load("ex.priors.rda")
  
  # solve for MRF parameters
  res <- solveMRF(ex.priors)
  
  # get posterior probabilities of class membership
  ex.post <- posterior(res$a, res$b, ex.priors)
  
  # get the posterior best classifications
  postb <- best.mod.list(ex.post)
  
  return(list(res = res, pp = ex.priors, postp = ex.post, postb = postb))
  
}
