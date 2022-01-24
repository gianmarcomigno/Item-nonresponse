The purpose of the work is to understand by using genetics, the behavior of those people who tend not to answer the questions in a 
questionnaire when they have the opportunity to do so (referred as non-responders). 
In this context, non-responders are those who prefer not to answer or to choose the option “I don’t know” in a question when available.

The first step is to develop a missingness phenotype to characterize the non-responders and to identify shared characteristics of individuals driving this bias. Afterwards, heritability has been estimated on the missingness phenotype defined for non-responders and genetic correlation has been calculated with other heritable phenotypes.

In the file "common analysis - PNA & IDK" the data are imported and both the participants and the relevant variables of interest have been subset. Only the participants who answered every question and the questions asked to everybody have been kept in the analysis. Moreover there are some functions which will be useful to plot some results. 

In the file named "Prefer not to answer - PNA" there is the part of the analysis on the items with the "Prefer not to answer" option. After recoding the variable as binary (1 if someone prefer not to answer, 0 otherwise), we run factor analysis with one factor with the idea of extracting off the general latent factor (non-responding behaviour) showing the substructure underneath it. We observed an item redundancy and tried different ways of pruning - aka reducing the dimensionality of the items. Then, we run confirmatory factor analysis with 5 factors, extracted the factor scores and considered the general latent factor as the "prefer not to answer" phenotype. 

Similarly, in the file named "I don't know - IDK" there is the analysis on the items with the "I don't know" option. The only difference is that we didn't pruned in the same way since the distribution of the PNA and IDK in the two analyses is quite different.
