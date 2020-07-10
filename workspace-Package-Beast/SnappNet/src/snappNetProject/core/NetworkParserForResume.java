package snappNetProject.core;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.ClusterTree;
import beast.core.Input.Validate;

/**
 * Parse the network of extended Newick format.
 * @author Huw Ogilvie
 */


@Description("Parse the network of extended Newick format.")
public class NetworkParserForResume extends Network implements StateNodeInitialiser {
    public final Input<Network> networkInput = new Input<>("initial", "Network to initialize.");
    public final Input<Tree> treeInput =
            new Input<>("tree", "Tree initialized from extended newick string.", Validate.REQUIRED);
    public final Input<Boolean> adjustTipHeightsInput =
            new Input<>("adjustTipHeights", "Whether tipHeights shall be adjusted (default is true).", true);

    private List<String> leafOrder;
    private int nextSpeciationNr;
    private int nextReticulationNr;

    public NetworkParserForResume() {
    }

    //CE le met en comment
    public NetworkParserForResume(final Tree tree) {
    		Log.debug.println("Je construit mon Network Parser depuis un arbre");//CE
    		treeInput.setValue(tree, this);
    		initAndValidate();
    }
     
    
    //CE j essaie autre chose
    @Override
    public void initAndValidate() {
    	
    //////////CE tries to put the original version of Chi Zhang ...
    		Log.debug.println("CE Je passe ds le NetworkParser version Chi Zhang\n");
    		final Tree tree = treeInput.get();  	
    	
        final Node treeRoot = tree.getRoot();

        // Step (1) is to initialize the node counts and array
        leafOrder = new ArrayList<>();
        leafNodeCount = 0;
        speciationNodeCount = 0;
        int hybridNodeCount = 0;
        Log.debug.println("Je suis ds Init and validate de Network Parser.java\n");//CE
	
		//CE
		Log.debug.println("Verif du newick avant la construction du reseau\n");
		for (Node n: tree.getNodesAsArray()) {
			Log.debug.println("Je suis le noeud n"+n.getID()+"\n");
			Log.debug.println("Voici ma hauteur"+n.getHeight()+"\n");
		}
		Log.debug.println("Fin de la Verif du newick avant la construction du reseau\n");
		//END CE
	
	
	
        for (Node n: tree.getNodesAsArray()) {
            if (n.getID() != null && n.getID().startsWith("#H")) {
                hybridNodeCount++;
            } else if (n.isLeaf()) {
            	leafOrder.add(n.getID());
                leafNodeCount++;
            } else if (!n.isRoot()) {
                speciationNodeCount++;
            }
        }

        leafOrder.sort(null);

        assert hybridNodeCount % 2 == 0;
        reticulationNodeCount = hybridNodeCount / 2;
        nodeCount = leafNodeCount + speciationNodeCount + reticulationNodeCount + 1;
        nodes = new NetworkNode[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            nodes[i] = new NetworkNode(this); 
        }

        nextSpeciationNr = leafNodeCount;
        nextReticulationNr = leafNodeCount + speciationNodeCount;
	Log.debug.printf("The next speciation number is (in NetworkParser.java)\n" + nextSpeciationNr + "\n");//CE
	Log.debug.printf("The next reticulation number is (in NetworkParser.java)\n" + nextReticulationNr  + "\n");//CE

        // Step (2) is to recursively copy the tree to the network
        rebuildNetwork(treeRoot);
	Log.debug.printf("We have now rebuilt the network (in NetworkParser.java)\n" );//CE


        // Update the cached parents and children for each node
        updateRelationships();

        // Step (3) adjust network tip height to ZERO
        if (adjustTipHeightsInput.get()) {
            // all nodes should be at zero height if no date-trait is available
            for (NetworkNode tip: getLeafNodes()) {
                tip.setHeight(0.0);
            }
        }

        super.initAndValidate();
    	
    	 
    }

    private Integer rebuildNetwork(final Node treeNode) {

	Log.debug.printf("Je suis dans rebuildNetwork (in NetworkParser.java)\n");
        Integer branchNumber;
        NetworkNode newNode;

        final String nodeLabel = treeNode.getID();
        final double nodeHeight = treeNode.getHeight();
        final int matchingNodeNr = getNodeNumber(nodeLabel);
		
	Log.debug.printf("Voici le nodeLabel\n" + nodeLabel + "\n");//CE
	Log.debug.printf("Voici le nodeHeight\n" + nodeHeight + "\n");//CE
	Log.debug.printf("Voici le matchingNodeNr\n" + matchingNodeNr + "\n");//CE

        if (matchingNodeNr < 0) {
            int newNodeNumber;
            double inheritProb = 0.5;

            if (treeNode.isRoot()) {
                newNodeNumber = nodeCount - 1;
            } else if (nodeLabel != null && nodeLabel.startsWith("#H")) {
                if (treeNode.getMetaDataNames().contains("gamma"))
                    inheritProb = (Double) treeNode.getMetaData("gamma");
                newNodeNumber = nextReticulationNr;
                nextReticulationNr++;
            } else if (treeNode.isLeaf()) {
                newNodeNumber = leafOrder.indexOf(nodeLabel);
            } else {
                newNodeNumber = nextSpeciationNr;
                nextSpeciationNr++;
            }
            
            newNode = nodes[newNodeNumber];
            newNode.label = nodeLabel;
            newNode.height = nodeHeight;
            newNode.inheritProb = inheritProb;

            branchNumber = getBranchNumber(newNodeNumber);
	    Log.debug.printf("Voici le newNodeNumber\n" + newNodeNumber + "\n");//CE
	    Log.debug.printf("Voici le branchNumber\n" + branchNumber + "\n");//CE

        } else {
	    Log.debug.printf("J ai un matchingNodeNr !!!\n");//CE
            newNode = nodes[matchingNodeNr];
            if (treeNode.getMetaDataNames().contains("gamma")){
		Log.debug.printf("Mon meta data contient gamma !!!\n");//CE
                newNode.inheritProb = 1.0 - (double)treeNode.getMetaData("gamma");
	    }

            branchNumber = getBranchNumber(matchingNodeNr) + 1;
        }

        for (Node child: treeNode.getChildren()) {
            final Integer childBranchNr = rebuildNetwork(child);
            newNode.childBranchNumbers.add(childBranchNr);
        }

        return branchNumber;
    }

    @Override
    public void initStateNodes() {
    		Log.debug.println("CE Je suis ds initStateNode de NetworkParser\n");
        if (networkInput.get() != null) {
        	Log.debug.println("CE Je passe ds le initStateNode de NetworkParser\n");
            networkInput.get().assignFrom(this);
        }
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
    		Log.debug.println("CE Je suis ds getInitialisedStateNodes de NetworkParser\n");
        if (networkInput.get() != null) {
        	Log.debug.println("CE Je passe ds le getInitialisedStateNodes de NetworkParser\n");
            stateNodes.add(networkInput.get());
        }
    }
}
