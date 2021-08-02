# SnappNet
SnappNet is a Bayesian method dedicated to phylogenetic network inference, and proposed in 2021

******************************************************************************************************************
SnappNet is described in the paper entitled "On the inference of complex phylogenetic networks by Markov Chain Monte-Carlo"
and published in Plos Computational Biology (2021).

        The authors are

                Charles-Elie Rabier, Vincent Berry, 
                Marnus Stoltz, Joao D Santos, 
                Wensheng Wang, Jean-Christophe Glaszmann, 
                Fabio Pardi, Scornavacca Celine 

        
        from

              ISE-M, Univ. Montpellier, CNRS, EPHE, IRD, Montpellier, France

              LIRMM, Univ.Montpellier, CNRS, Montpellier, France

              IMAG, Univ. Montpellier, CNRS, Montpellier, France

              AGAP, CIRAD, Montpellier, France 
              
              ICS, Chinese Academy of Agricultural Sciences, Beijing, China

***********************************************************************************************************************

SnappNet is a new Bayesian method that directly relies on DNA sequences. Our method is implemented in BEAST 2 (Bouckaert et al., 2014) , an improved version of the popular version BEAST 1.x dedicated to Bayesian evolutionary analyses. Our SnappNet package is built on two BEAST packages, Snapp (Bryant et al, 2012), and SpeciesNetwork (Zhang et al., 2017). It incorporates the novel MCMC operators of SpeciesNetwork to move through the network space, and also benefits from operators specific to the mathematical model behind Snapp (e.g. population sizes, mutation rates ...) and extended to handle networks. 

*************************************************************************************************************************

In the folder workspace-Package-Beast/SnappNet/doc, there are some informations on how to run SnappNet.
In the folder workspace-Package-Beast/SnappNet/example, you will find the rice real data.

SnappNet current version requires Beast 2.6.1

***************************************************************************************************************************
Some extra informations are given at http://charles-elie.rabier.pagesperso-orange.fr/doc/SnappNet.html . Be careful xml informations are linked to another project
called SnappNetForSimSnappNet (also available on github) that we used for analyzing simulated data
