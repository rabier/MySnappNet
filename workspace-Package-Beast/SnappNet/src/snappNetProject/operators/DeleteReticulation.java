package snappNetProject.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import snappNetProject.core.Network;
import snappNetProject.core.NetworkNode;
import snappNetProject.core.SanityChecks;

/**
 * This proposal  deletes a reticulation branch from the species network. If there is no reticulation, or if the branch
 * is connecting two reticulation nodes, this is aborted. The two branches at each connecting point are joined,
 * resulting branches with length l1 and l2 respectively. The gamma prob r is removed.
 * The Jacobian is 1 / (l1 * l2).
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob.
 * Let m be the number of reticulation branches in the current network. The probability of selecting the this branch to
 * remove is 1/m.
 * Let k be the number of branches in the proposed network. The probability of adding this branch is (1/k)(1/k)
 * The Hastings ratio is (1/k)(1/k)(1)(1)(1) / (1/m) = m / k^2.
 *
 * See also AddReticulation.
 *
 * @author Chi Zhang slightly modified by CE Rabier 
 * in order to handle coalescent rates
 */

@Description("Delete a reticulation branch from the species network.")
public class DeleteReticulation extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate", "population sizes");
    
    
    
    // empty constructor to facilitate construction by XML + initAndValidate
    public DeleteReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());               
        
        RealParameter coalescenceRate = m_coalescenceRate.get(this);
        Double [] coalescenceRateValues = coalescenceRate.getValues();        
                       
        final int nHybridNodes = speciesNetwork.getReticulationNodeCount();
        if (nHybridNodes == 0)  // there is no reticulation branch to delete
            return Double.NEGATIVE_INFINITY;
        // number of reticulation branches in the current network
        final int nReticulationBranches = 2 * nHybridNodes;  // m'

        // pick a reticulation branch randomly
        final Integer hybridBranchNr = Randomizer.nextInt(nReticulationBranches) + speciesNetwork.getReticulationOffset();        
        
        final int hybridNodeNr = speciesNetwork.getNodeNumber(hybridBranchNr);
        // branch with hybridBranchNr is connecting hybridNode and parentNode
        NetworkNode hybridNode = speciesNetwork.getNode(hybridNodeNr);
        NetworkNode parentNode = hybridNode.getParentByBranch(hybridBranchNr);
        if (parentNode.isReticulation())  // cannot delete a branch connecting two reticulation nodes
            return Double.NEGATIVE_INFINITY;

        
        // get the parent node and another child node of parentNode
        final Integer pNParentBranchNr = parentNode.gammaBranchNumber;
        NetworkNode pNParentNode = parentNode.getParentByBranch(pNParentBranchNr);
        final Integer pNChildBranchNr;
        if (parentNode.childBranchNumbers.get(0).equals(hybridBranchNr))
            pNChildBranchNr = parentNode.childBranchNumbers.get(1);
        else
            pNChildBranchNr = parentNode.childBranchNumbers.get(0);
        NetworkNode pNChildNode = parentNode.getChildByBranch(pNChildBranchNr);

        // get the child node and another parent node of hybridNode
        final Integer hNChildBranchNr = hybridNode.childBranchNumbers.get(0);
        NetworkNode hNChildNode = hybridNode.getChildByBranch(hNChildBranchNr);
        final Integer hNParentBranchNr;
        if (hybridNode.gammaBranchNumber.equals(hybridBranchNr))
            hNParentBranchNr = hybridNode.gammaBranchNumber + 1;
        else
            hNParentBranchNr = hybridNode.gammaBranchNumber;
        NetworkNode hNParentNode = hybridNode.getParentByBranch(hNParentBranchNr);

        // work out the Jacobian
        final double l1, l2;
        if (parentNode == hNParentNode && hybridNode == pNChildNode) {
            // the two attaching points are on the same branch
            l1 = l2 = pNParentNode.getHeight() - hNChildNode.getHeight();
        } else {
            // the two attaching points are on different branches
            l1 = hNParentNode.getHeight() - hNChildNode.getHeight();
            l2 = pNParentNode.getHeight() - pNChildNode.getHeight();
        }
        double logProposalRatio = - Math.log(l1) - Math.log(l2);

        
        //Handle Coalescent rates at this time , the topology network hasn't changed yet
        Double[] values = new Double[speciesNetwork.getBranchCount()-3]; 	
        
        for (NetworkNode node: speciesNetwork.getLeafNodes() ) {
         //handle leaf edges          	
        	values[node.getNr()]=coalescenceRateValues[node.getNr()];
        }
        	       
        //handle other edges  (based on deleteReticulationBranch in class Network)
        //CE: I handle coalescentRate outside Network like Bryant and Bouckaert
        // that s why i am doing this here !!!
        NetworkNode bifurcNode = parentNode;
        int bifurcNodeNr = bifurcNode.getNr();  
        int bifurcNodeGammaBranchNumber=pNParentBranchNr;
      
        for (int edgeNr=speciesNetwork.getLeafNodeCount(); edgeNr< speciesNetwork.getBranchCount(); edgeNr++) {
        
        	if ( (edgeNr!=bifurcNodeGammaBranchNumber) && (edgeNr!=hybridBranchNr) && (edgeNr!=hNParentBranchNr)) {
        		//edgeNr is a branch which won t disappear in the new topology 
        		if (edgeNr > bifurcNodeNr && edgeNr < hybridBranchNr)
        				values[edgeNr-1]= coalescenceRateValues[edgeNr];
        		else if (edgeNr > hybridBranchNr)
        				values[edgeNr-3]= coalescenceRateValues[edgeNr];
        		else    
        				values[edgeNr]= coalescenceRateValues[edgeNr];
        				;
        	}
    	}
                  	
    	//start editing coalescenceRate
    	coalescenceRate.startEditing(this);    	
    	coalescenceRate.setDimension(values.length);    	
    	
    	for (int i=0; i< values.length; i++) {
    		coalescenceRate.setValue(i, values[i]);    		
    	}    	 
                             
    	// start moving
    	speciesNetwork.startEditing(this);

    	// delete the reticulation branch
    	speciesNetwork.deleteReticulationBranch(hybridBranchNr);

    	SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());
                                 
    	// number of branches in the proposed network
    	final int nBranches = speciesNetwork.getBranchCount();  // k'
    	logProposalRatio += Math.log(nReticulationBranches) - 2 * Math.log(nBranches);

    	return logProposalRatio;
    }
}
