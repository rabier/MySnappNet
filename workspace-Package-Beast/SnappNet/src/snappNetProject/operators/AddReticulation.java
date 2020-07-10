package snappNetProject.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.util.Randomizer;
import snappNetProject.core.Network;
import snappNetProject.core.NetworkNode;
import snappNetProject.core.SanityChecks;

/**
 * This proposal  adds a reticulation event by connecting two existing branches (with length l1 and l2) with a new branch.
 * The same branch can be picked twice (and forms a loop to that branch). The cutting proportion of each picked branch by
 * the connecting point, w1 and w2 ~ Uniform(0,1). Let l11 = l1 * w1, l12 = l1 * (1-w1), l21 = l2 * w2, l22 = l2 * (1-w2)
 * The direction of the new branch is determined by the two connecting points, the higher is speciation node, and the
 * lower is reticulation node. The gamma prob r = w3 ~ Uniform(0,1).
 * The Jacobian is l1 * l2.
 *
 * The AddReticulation and DeleteReticulation are chosen with equal prob. If there is no reticulation in the network,
 * the DeleteReticulation move is aborted.
 * Let k be the number of branches in the current network. The probability of adding this branch is (1/k)(1/k)
 * Let m be the number of reticulation branches in the proposed network. The probability of selecting the same branch to
 * remove is (1/m).
 * The Hastings ratio is (1/m) / [(1/k)(1/k)(g1)(g2)(g3)] = k^2 / m, with g1 = g2 = g3 = 1 (uniform density).
 *
 * See also DeleteReticulation.
 *
 * @author Chi Zhang, slightly modified by CE Rabier 
 * in order to bound the number of reticulations, and to handle coalescent rates
 */
 

@Description("Add a reticulation branch to the species network.")
public class AddReticulation extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public Input<RealParameter> m_coalescenceRate = new Input<RealParameter>("coalescenceRate", "population sizes");
    public Input<Boolean> m_boundReticulationNumber = new Input<Boolean>("boundReticulationNumber", 
			"Check box only if you want to bound the number of reticulations.",
			true);
    public Input<Integer> m_maxReticulationNumber = new Input<Integer>("maxReticulationNumber", 
			"This number will be taken into account only if the box above is checked", 2);
        
    // empty constructor to facilitate construction by XML + initAndValidate
    public AddReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {

        final Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());
        
        final Boolean boundReticulationNumber=m_boundReticulationNumber.get();
        
        final int nHybridNodes=speciesNetwork.getReticulationNodeCount();
        
        if (boundReticulationNumber) {
        	    Integer maxReticulationNumber=m_maxReticulationNumber.get();
        		if (nHybridNodes == maxReticulationNumber)  
                return Double.NEGATIVE_INFINITY;
        }
        
        speciesNetwork.startEditing(this);
                
        // number of branches in the current network
        final int nBranches = speciesNetwork.getBranchCount();  // k
 
        RealParameter coalescenceRate = m_coalescenceRate.get(this);
        Double [] coalescenceRateValues = coalescenceRate.getValues();  
                                 
        // pick two branches randomly, including the root branch
        final Integer pickedBranchNr1 = Randomizer.nextInt(nBranches);
        final Integer pickedBranchNr2 = Randomizer.nextInt(nBranches);  // allow picking the same branch

        // get the nodes associated with each branch
        final int pickedNodeNr1 = speciesNetwork.getNodeNumber(pickedBranchNr1);
        NetworkNode pickedNode1 = speciesNetwork.getNode(pickedNodeNr1);
        final int pickedNodeNr2 = speciesNetwork.getNodeNumber(pickedBranchNr2);
        NetworkNode pickedNode2 = speciesNetwork.getNode(pickedNodeNr2);
        NetworkNode pickedParent1 = pickedNode1.getParentByBranch(pickedBranchNr1);
        NetworkNode pickedParent2 = pickedNode2.getParentByBranch(pickedBranchNr2);

        // propose the attaching position at each branch
        final double l1, l2, l11, l21;
        l1 = pickedParent1.getHeight() - pickedNode1.getHeight();
        l11 = l1 * Randomizer.nextDouble();
        l2 = pickedParent2.getHeight() - pickedNode2.getHeight();
        l21 = l2 * Randomizer.nextDouble();

        double logProposalRatio = Math.log(l1) + Math.log(l2);  // the Jacobian

        // start moving
        speciesNetwork.startEditing(this);

        // create two new nodes
        NetworkNode middleNode1 = new NetworkNode(speciesNetwork);
        NetworkNode middleNode2 = new NetworkNode(speciesNetwork);
        // set height
        middleNode1.setHeight(pickedNode1.getHeight() + l11);
        middleNode2.setHeight(pickedNode2.getHeight() + l21);

        // add a branch joining the two middle nodes (picked branches)
        if (middleNode1.getHeight() < middleNode2.getHeight()) {
            speciesNetwork.addReticulationBranch(middleNode1, middleNode2, pickedBranchNr1, pickedBranchNr2);
            middleNode1.setGammaProb(Randomizer.nextDouble());
        } else {
            speciesNetwork.addReticulationBranch(middleNode2, middleNode1, pickedBranchNr2, pickedBranchNr1);
            middleNode2.setGammaProb(Randomizer.nextDouble());
        }

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // number of reticulation branches in the proposed network
        final int nReticulationBranches = 2 * speciesNetwork.getReticulationNodeCount();  // m
        logProposalRatio += 2 * Math.log(nBranches) - Math.log(nReticulationBranches);
        
                
        //Let us handle coalescent rates                        
        Double[] values = new Double[speciesNetwork.getBranchCount()]; 		
                      
        //be careful : the speciation nodes and their branches have the same numbers except the new speciation node
       for (NetworkNode node: speciesNetwork.getLeafNodes()) {
        	int BranchNumber = node.gammaBranchNumber;  
        	values[BranchNumber] = new Double(coalescenceRateValues[BranchNumber]); 
        }
        
        
     	for (NetworkNode node: speciesNetwork.getSpeciationNodes()) {
    	   	     
    	   	     int BranchNumber = node.gammaBranchNumber;   	   	         	   	     
    	   	     
    	   	     if (BranchNumber!=(speciesNetwork.getReticulationOffset()-1)) {
    	   	    	 //handle an old speciation node 
    	   	    	 //so the branch number did not change between the new and the old network
    	   	    	 values[BranchNumber] = new Double(coalescenceRateValues[BranchNumber]);   
    	   	     }else {
    	   	    	//handle the new branch above the new speciation node
    	   	    	 //we set the same value as coalescentrate as the branch above the speciation child of node 
    	   	    	int branchFirstChild=node.childBranchNumbers.get(0);
    	   	    int branchSecondChild=node.childBranchNumbers.get(1);  	   
    	   	        
    	   	    	NetworkNode node1=node.getChildByBranch(branchFirstChild);
    	   	    	NetworkNode node2=node.getChildByBranch(branchSecondChild);
    	   	    	int idFirstChild=node1.getNr();
    	   	        int idSecondChild=node2.getNr();
    	   	        
    	   	        if (idFirstChild != idSecondChild) {
    	   	        // we picked two different branches for adding the reticulate node and the new speciation node
    	   	        	
    	   	        	if (!(node1.isReticulation() & node2.isReticulation())){
    	   	        	// one child is a speciation, the other one a reticulate node
    	   	        	int branchSpeciation= node1.isSpeciation() ? branchFirstChild: branchSecondChild;
    	   	        	Log.debug.printf("Voila la valeur de la branche Speciation" +  branchSpeciation + "\n");
    	   	        	values[BranchNumber] = new Double(coalescenceRateValues[branchSpeciation]);   	   	    	
    	   	        	//since branchSpeciation has the same number in the old network 
    	   	        	
    	   	        	}else {
    	   	        	// both children are reticulate nodes, yes it can happen !!!    	   	        		
    	   	        		int branchOldRetic = (idFirstChild > idSecondChild) ? branchFirstChild : branchSecondChild;  	   	        		
    	   	        		values[BranchNumber] = new Double(coalescenceRateValues[branchOldRetic-3]);  
    	   	        		//branchOldRetic in the new network had the number branchOldRetic-3 in the old network
    	   	        	}
    	   	        
   	   	        
    	   	        }else {
    	   	        // we picked the same branches for adding the reticulate node and the new speciation node
    	   	        	
    	   	        	NetworkNode Child=node.getChildByBranch(branchFirstChild);    	   	        	
    	   	        	int branchChildOfChild=Child.childBranchNumbers.get(0);    	   	      
    	   	        	
    	   	        	if (Child.getChildByBranch(branchChildOfChild).isReticulation()) {
    	   	        		//child of child is a reticulate node
    	   	        		values[BranchNumber] = new Double(coalescenceRateValues[BranchNumber]); 
    	   	        	}else {
    	   	        	//child of child is a speciation node or a leaf
    	   	        		values[BranchNumber] = new Double(coalescenceRateValues[branchChildOfChild]);   	   	        		
    	   	        	}
    	  	        	
    	   	        	 
    	   	        }
    	   	    	   	   	    	
    	   	     }
    	   	    	   	   	      	   	    	 
    	 }
     	//it is okay since the new speciation node is the last speciation node, ie with the highest number
     	      	
     	
     	//handle now branches above retic nodes 
     	for (NetworkNode node: speciesNetwork.getReticulationNodes()) {
     	
     		int BranchNumber = node.gammaBranchNumber;  
     		     		
     		if (node.getNr()!=speciesNetwork.getReticulationOffset()) {
     			//handle an old retic node
     			values[BranchNumber] = new Double(coalescenceRateValues[BranchNumber-3]);
     			values[BranchNumber+1] = new Double(coalescenceRateValues[BranchNumber-2]);
     		}else {     			     			
     			//handle the new retic node   	
     			int branchChild=node.childBranchNumbers.get(0);
     			//branchChild keep the same number
     			
     			if (!node.getChildByBranch(branchChild).isReticulation()) {
     			    //the children of node is a speciation node
     				//so the brancChild remained the same between the old and the new network
     				values[BranchNumber] = new Double(coalescenceRateValues[branchChild]);
     				values[BranchNumber+1] = new Double(coalescenceRateValues[branchChild]);   			
     			
     			} else {
     				//the children of node is a reticulation
     				// so the branchCHild has changed between the old and the new network
     				values[BranchNumber] = new Double(coalescenceRateValues[branchChild-3]);
     				values[BranchNumber+1] = new Double(coalescenceRateValues[branchChild-3]);    				     				
     			}
     			
     		}
     		  		     		
     	}	
     		   	     
    	
    	//start editing coalescenceRate
    	coalescenceRate.startEditing(this);
    	
    	coalescenceRate.setDimension(values.length);
    	
    	for (int i=0; i< values.length; i++) {
    		coalescenceRate.setValue(i, values[i]);    		
    	}    	
              
        
    	Log.debug.printf("Voici le logProposalRatio" + logProposalRatio +"\n");
        return logProposalRatio;
        
    }
}
