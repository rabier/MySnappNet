// SnappNetPrior  is largely inspired by SnAPPrior.java of Bryant and Bouckaert
//we have kept only the elements required for SnappNet 

//we have a tree as input , just for Beauty ... but we are not using it

//we compute the prior on the coalescent rate on each edge of the network

package snappNetProject.core;

import java.util.List;
import java.util.Random;

//import snap.distribution.ChiSquareNoncentralDist;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree; 

 
@Description("Standard prior for SnAPNet analysis, consisting  " +
        "gamma distribution over the theta values " +
        "(with parameters alpha and beta). " +
        "Thetas are represented by the coalescenceRate parameter where values are theta=2/coalescenceRate")
public class SnappNetPrior extends Distribution {
    public Input<RealParameter> m_pAlpha = new Input<RealParameter>("alpha", "Alpha parameter for the gamma prior on population size (theta) values", Validate.REQUIRED);
    public Input<RealParameter> m_pBeta = new Input<RealParameter>("beta", "Beta parameter for the gamma prior on population size (theta) values", Validate.REQUIRED);
    public Input<RealParameter> m_pCoalescenceRate = new Input<RealParameter>("coalescenceRate", "Populations sizes for the nodes in the tree", Validate.REQUIRED);     
    public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations"); //, Validate.REQUIRED);

    
    @Override
    public double calculateLogP() {
        logP = 0.0;

        double alpha = m_pAlpha.get().getValue();
        double beta = m_pBeta.get().getValue();              
        //Gamma values in tree
        RealParameter coalescenceRate = m_pCoalescenceRate.get();
			
        // we have chosen PRIORCHOICE == 0 in the original code
		//Assume independent gamma distributions for thetas.
			
		//We assume that 2/r has a gamma(alpha,beta) distribution. That means that r has density proportional to
		// 1/(r^2)  * GAMMA(2/r|alpha,beta)
		//which has log (alpha - 1.0)*Math.log(2.0/r) - (beta *(2.0/ r)) - 2*log(r), which in turn simplifies to the expr. below (w/ consts)
		
		//CE in our case, Nodeshave have been replaced by branches
		for (int iEdge = 0; iEdge < coalescenceRate.getDimension(); iEdge++) {
				double r = coalescenceRate.getValue(iEdge);
				logP += -(alpha + 1.0)*Math.log(r) - 2.0* beta / r;								
			}
			 
        return logP;
    } // calculateLogLikelihood

		
	// below is original code from Bryant and Bouckaert, not sure it is useful for SnappNet	
    double heightSum(Node node) {
        if (node.isLeaf()) {
            return 0;
        } else {
            double h = node.getHeight();
            h += heightSum(node.getLeft());
            if (node.getRight() != null) {
            	h += heightSum(node.getRight());
            }
            return h;
        }
    } // heightSum

    //Returns a list of branching times in the tree, sorted in an decreasing sequence. First one is
    //the height of the mrca of the tree.
    List<Double> getSortedHeights(Node node) {
        return null;
    }

    private double p0(double t, double lambda, double mu) {
        return mu*(1-Math.exp(-(lambda-mu)*t))/(lambda - mu*Math.exp(-(lambda-mu)*t));
    }

    private double p1(double t, double lambda, double mu) {
        double denominator = (lambda - mu*Math.exp(-(lambda-mu)*t));
        return (lambda - mu)*(lambda - mu) * Math.exp(-(lambda-mu)*t) / (denominator * denominator);
    }

	@Override public List<String> getArguments() {return null;}
	@Override public List<String> getConditions() {return null;}
	@Override public void sample(State state, Random random) {};
} // class SSSPrior 
