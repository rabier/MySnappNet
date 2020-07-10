# MySnappNet
SnappNet is a new Bayesian method dedicated to phylogenetic network inference

by Rabier, Berry, Glaszmann, Pardi and Scornavacca

ISEM, CIRAD, LIRMM


SnappNet is a new Bayesian method that directly relies on DNA sequences. Our method is implemented in BEAST 2 (Bouckaert et al., 2014) , an improved version of the popular version BEAST 1.x dedicated to Bayesian evolutionary analyses. Our SnappNet package is built on two BEAST packages, Snapp (Bryant et al, 2012), and SpeciesNetwork (Zhang et al., 2017). It incorporates the novel MCMC operators of SpeciesNetwork to move through the network space, and also benefits from operators specific to the mathematical model behind Snapp (e.g. population sizes, mutation rates ...) and extended to handle networks. 
