package snappNetProject.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.util.Randomizer;
import snappNetProject.core.Network;
import snappNetProject.core.NetworkNode;
import snappNetProject.core.SanityChecks;

/**
 * @author Chi Zhang
 */

@Description("Flip  the direction of a reticulation branch.")
public class FlipReticulation extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);

    // empty constructor to facilitate construction by XML + initAndValidate
    public FlipReticulation() {
    }

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());   
        
        final int nHybridNodes = speciesNetwork.getReticulationNodeCount();
        if (nHybridNodes == 0)  // there is no reticulation branch to flip
            return Double.NEGATIVE_INFINITY;

        // pick a reticulation branch randomly
        final Integer hybridBranchNr = Randomizer.nextInt(2 * nHybridNodes) + speciesNetwork.getReticulationOffset();
        
     
        final int hybridNodeNr = speciesNetwork.getNodeNumber(hybridBranchNr);
        // branch with hybridBranchNr is connecting hybridNode and parentNode
        NetworkNode hybridNode = speciesNetwork.getNode(hybridNodeNr);
        NetworkNode parentNode = hybridNode.getParentByBranch(hybridBranchNr);
        if (parentNode.isReticulation())  // cannot flip a branch connecting two reticulation nodes
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

        // do not deal with bubble
        if (parentNode == hNParentNode && hybridNode == pNChildNode)
            return Double.NEGATIVE_INFINITY;

        final double hNParentHeight = hNParentNode.getHeight();
        final double pNChildHeight = pNChildNode.getHeight();
        final double hNChildHeight = hNChildNode.getHeight();
        final double lower = Math.max(pNChildHeight, hNChildHeight);
        if (hNParentHeight <= lower)  // cannot flip
            return Double.NEGATIVE_INFINITY;
        // new speciation node height
        final double spHeight = lower + (hNParentHeight - lower) * Randomizer.nextDouble();

        final double pNParentHeight = pNParentNode.getHeight();
        final double upperF = Math.min(spHeight, pNParentHeight);
        // new hybridization node height
        final double hyHeight = pNChildHeight + (upperF - pNChildHeight) * Randomizer.nextDouble();

        // calculate proposal ratio
        final double upperB = Math.min(parentNode.getHeight(), hNParentHeight);
        final double logProposalRatio = Math.log(hNParentHeight - lower) + Math.log(upperF - pNChildHeight)
                                      - Math.log(pNParentHeight - lower) - Math.log(upperB - hNChildHeight);

        //start moving
        speciesNetwork.startEditing(this);

        parentNode.setHeight(spHeight);
        hybridNode.setHeight(hyHeight);
        hNParentNode.childBranchNumbers.remove(hNParentBranchNr);
        hNParentNode.childBranchNumbers.add(pNParentBranchNr);
        parentNode.childBranchNumbers.remove(pNChildBranchNr);
        parentNode.childBranchNumbers.add(hNChildBranchNr);
        pNParentNode.childBranchNumbers.remove(pNParentBranchNr);
        pNParentNode.childBranchNumbers.add(hNParentBranchNr);
        hybridNode.childBranchNumbers.remove(hNChildBranchNr);
        hybridNode.childBranchNumbers.add(pNChildBranchNr);
        hNParentNode.updateRelationships();
        parentNode.updateRelationships();
        hNChildNode.updateRelationships();
        pNParentNode.updateRelationships();
        hybridNode.updateRelationships();
        pNChildNode.updateRelationships();

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());
        return logProposalRatio;
    }
}
