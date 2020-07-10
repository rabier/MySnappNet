package snappNetProject.core;
 
import java.util.List;

import beast.core.util.Log;
import beast.evolution.alignment.TaxonSet;

public class SnappNetLikelihoodCore {

	SiteProbabilityCalculator [] m_siteProbabilityCalculator;
	int numPatterns;
	
	public SnappNetLikelihoodCore(Network speciesNetwork, SnapData data) {    	    	     
		
		numPatterns = data.getPatternCount();
        m_siteProbabilityCalculator = new SiteProbabilityCalculator[numPatterns];
        
        for(int id = 0; id < numPatterns; id++) {
        	m_siteProbabilityCalculator[id]= new SiteProbabilityCalculator(speciesNetwork);
        }
        
	}
	
	
	public double [] computeLogLikelihood(SnapData data, Network speciesNetwork, double u, double v, 
            Double [] coalescenceRate) throws Exception
	{
						
		int numPatterns = data.getPatternCount();
		
		//Temporarily store pattern probabilities... used for numerical checks.
        double [] patternProb = new double[numPatterns];
        List<TaxonSet> taxonSets=data.m_taxonsets.get();
                 
		for(int id = 0; id < numPatterns; id++) {
            int [] dataAtThisSite = data.getPattern(id);
            int [] lineageCounts = data.getPatternLineagCounts(id);	               
            patternProb[id] = m_siteProbabilityCalculator[id].computeSiteLikelihood(dataAtThisSite, taxonSets, lineageCounts, speciesNetwork, u, v, coalescenceRate);            
            m_siteProbabilityCalculator[id]=null;            
        }
	
		return patternProb;
	
		
	}	
		
	
}
