# MySnappNet
SnappNet is a new Bayesian method dedicated to phylogenetic network inference

******************************************************************************************************************
SnappNet is described in the submitted paper "On the inference of complex phylogenetic networks by Markov Chain Monte-Carlo"

        by Charles-Elie Rabier, Vincent Berry, Jean-Christophe Glaszmann, Fabio Pardi, Scornavacca Celine 


              ISE-M, Univ. Montpellier, CNRS, EPHE, IRD, Montpellier, France

              LIRMM, Univ.Montpellier, CNRS, Montpellier, France

              IMAG, Univ. Montpellier, CNRS, Montpellier, France

              AGAP, CIRAD, Montpellier, France 

***********************************************************************************************************************

SnappNet is a new Bayesian method that directly relies on DNA sequences. Our method is implemented in BEAST 2 (Bouckaert et al., 2014) , an improved version of the popular version BEAST 1.x dedicated to Bayesian evolutionary analyses. Our SnappNet package is built on two BEAST packages, Snapp (Bryant et al, 2012), and SpeciesNetwork (Zhang et al., 2017). It incorporates the novel MCMC operators of SpeciesNetwork to move through the network space, and also benefits from operators specific to the mathematical model behind Snapp (e.g. population sizes, mutation rates ...) and extended to handle networks. 

*************************************************************************************************************************

In the folder workspace-Package-Beast/SnappNet/doc, there are some informations on how to run SnappNet.
In the folder workspace-Package-Beast/SnappNet/example, you will find the rice real data.

SnappNet current version requires Beast 2.6.1

***************************************************************************************************************************
Some extra informations are given at http://charles-elie.rabier.pagesperso-orange.fr/doc/SnappNet.html . Be careful xml informations are linked to another project
called SnappNetForSimSnappNet (also available on github) that we used for analyzing simulated data
