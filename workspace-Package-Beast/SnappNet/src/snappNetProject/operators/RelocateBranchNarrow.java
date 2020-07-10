package snappNetProject.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.util.Randomizer;
import snappNetProject.core.Network;
import snappNetProject.core.NetworkNode;
import snappNetProject.core.SanityChecks;

/**
 * Pick an internal  network node randomly.
 * If bifurcation node: move the end of either the two parent branches.
 * If reticulation node: move the top of either the two child branches.
 * Propose a destination branch randomly from all possible branches, and attach the picked branch to the new position.
 *
 * @author Chi Zhang
 */

@Deprecated
@Description("Relocate the source of an edge starting with speciation node, " +
             "or the destination of an edge ending with hybridization node.")
public class RelocateBranchNarrow extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // pick an internal node randomly
        final NetworkNode[] internalNodes = speciesNetwork.getInternalNodes();
        int randomIndex = Randomizer.nextInt(internalNodes.length);
        NetworkNode pickedNode = internalNodes[randomIndex];

        // start moving
        speciesNetwork.startEditing(this);

        double logProposalRatio = 0.0;

        if (pickedNode.isReticulation()) {
            // move the end of either the two parent branches
            final Integer pickedParentBrNr, otherParentBrNr;
            if (Randomizer.nextBoolean()) {
                pickedParentBrNr = pickedNode.gammaBranchNumber;
                otherParentBrNr = pickedNode.gammaBranchNumber + 1;
            } else {
                pickedParentBrNr = pickedNode.gammaBranchNumber + 1;
                otherParentBrNr = pickedNode.gammaBranchNumber;
            }

            NetworkNode otherParent = pickedNode.getParentByBranch(otherParentBrNr);
       //     NetworkNode pickedParent = pickedNode.getParentByBranch(pickedParentBrNr);
            //final double upperLimit = pickedParent.getHeight();  // upper bound

            final Integer pickedChildBrNr = pickedNode.childBranchNumbers.get(0);
            NetworkNode pickedChild = pickedNode.getChildByBranch(pickedChildBrNr);

            // look for all the candidate branches to attach to
            List<Integer> candidateBranchNumbers = new ArrayList<>();
            /*
            if (isWide) {
                for (NetworkNode node : speciesNetwork.getAllNodesExceptOrigin()) {
                    // do not attach to the original position
                    if (node != pickedNode && node != pickedChild && node.getHeight() < upperLimit) {
                        candidateBranchNumbers.add(node.gammaBranchNumber);
                        if (node.isReticulation())
                            candidateBranchNumbers.add(node.gammaBranchNumber + 1);
                    }
                }
                if (upperLimit > otherParent.getHeight())
                    logProposalRatio -= Math.log(otherParent.getHeight() - pickedChild.getHeight());
                else
                    logProposalRatio -= Math.log(upperLimit - pickedChild.getHeight());
            } */
            for (NetworkNode node : speciesNetwork.getAllNodesExceptOrigin()) {
                // do not attach to the original position
                if (node != pickedNode && node != pickedChild && node.getHeight() < pickedNode.getHeight()) {
                    NetworkNode gParent = node.getParentByBranch(node.gammaBranchNumber);
                    NetworkNode oParent = node.getParentByBranch(node.gammaBranchNumber + 1);
                    if (gParent.getHeight() > pickedNode.getHeight())
                        candidateBranchNumbers.add(node.gammaBranchNumber);
                    if (node.isReticulation() && oParent.getHeight() > pickedNode.getHeight())
                        candidateBranchNumbers.add(node.gammaBranchNumber + 1);
                }
            }

            // pick a candidate branch randomly
            if (candidateBranchNumbers.size() == 0)
                return Double.NEGATIVE_INFINITY;
            randomIndex = Randomizer.nextInt(candidateBranchNumbers.size());
            final Integer attachBranchNr = candidateBranchNumbers.get(randomIndex);

            final Integer attachChildNr = speciesNetwork.getNodeNumber(attachBranchNr);
            NetworkNode attachChild = speciesNetwork.getNode(attachChildNr);
            NetworkNode attachParent = attachChild.getParentByBranch(attachBranchNr);

            /*
            if (isWide) {
                final double upper;
                if (upperLimit > attachParent.getHeight())
                    upper = attachParent.getHeight();
                else
                    upper = upperLimit;
                final double lower = attachChild.getHeight();
                // propose an attachment height
                final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();
                // set new height
                pickedNode.setHeight(newHeight);
                logProposalRatio += Math.log(upper - lower);
            } */

            // deal with the node relationships
            otherParent.childBranchNumbers.remove(otherParentBrNr);
            otherParent.childBranchNumbers.add(pickedChildBrNr);
            attachParent.childBranchNumbers.remove(attachBranchNr);
            attachParent.childBranchNumbers.add(otherParentBrNr);
            pickedNode.childBranchNumbers.remove(pickedChildBrNr);
            pickedNode.childBranchNumbers.add(attachBranchNr);
            otherParent.updateRelationships();
            pickedChild.updateRelationships();
            attachParent.updateRelationships();
            attachChild.updateRelationships();
            pickedNode.updateRelationships();

        } else {
            // move the top of either the two child branches
            final Integer pickedChildBrNr, otherChildBrNr;
            if (Randomizer.nextBoolean()) {
                pickedChildBrNr = pickedNode.childBranchNumbers.get(0);
                otherChildBrNr = pickedNode.childBranchNumbers.get(1);
            } else {
                pickedChildBrNr = pickedNode.childBranchNumbers.get(1);
                otherChildBrNr = pickedNode.childBranchNumbers.get(0);
            }

            NetworkNode otherChild = pickedNode.getChildByBranch(otherChildBrNr);
          //  NetworkNode pickedChild = pickedNode.getChildByBranch(pickedChildBrNr);
           // final double lowerLimit = pickedChild.getHeight();  // lower bound

            final Integer pickedParentBrNr = pickedNode.gammaBranchNumber;
            NetworkNode pickedParent = pickedNode.getParentByBranch(pickedParentBrNr);

            // look for all the candidate branches to attach to
            List<Integer> candidateBranchNumbers = new ArrayList<>();
            /*
            if (isWide) {
                for (NetworkNode node : speciesNetwork.getInternalNodesWithOrigin()) {
                    // do not attach to the original position
                    if (node != pickedNode && node.getHeight() > lowerLimit) {
                        for (Integer childBrNr : node.childBranchNumbers) {
                            if (node.getChildByBranch(childBrNr) != pickedNode)
                                candidateBranchNumbers.add(childBrNr);
                        }
                    }
                }
                if (lowerLimit < otherChild.getHeight())
                    logProposalRatio -= Math.log(pickedParent.getHeight() - otherChild.getHeight());
                else
                    logProposalRatio -= Math.log(pickedParent.getHeight() - lowerLimit);
            } */
            for (NetworkNode node : speciesNetwork.getInternalNodes()) {
                // do not attach to the original position
                if (node != pickedNode && node.getHeight() > pickedNode.getHeight()) {
                    for (Integer childBrNr : node.childBranchNumbers) {
                        NetworkNode gChild = node.getChildByBranch(childBrNr);
                        if (gChild != pickedNode && gChild.getHeight() < pickedNode.getHeight())
                            candidateBranchNumbers.add(childBrNr);
                    }
                }
            }

            // pick a candidate branch randomly
            if (candidateBranchNumbers.size() == 0)
                return Double.NEGATIVE_INFINITY;
            randomIndex = Randomizer.nextInt(candidateBranchNumbers.size());
            final Integer attachBranchNr = candidateBranchNumbers.get(randomIndex);

            final Integer attachChildNr = speciesNetwork.getNodeNumber(attachBranchNr);
            NetworkNode attachChild = speciesNetwork.getNode(attachChildNr);
            NetworkNode attachParent = attachChild.getParentByBranch(attachBranchNr);

            /*
            if (isWide) {
                // propose an attachment height
                final double lower;
                if (lowerLimit < attachChild.getHeight())
                    lower = attachChild.getHeight();
                else
                    lower = lowerLimit;
                final double upper = attachParent.getHeight();
                final double newHeight = lower + (upper - lower) * Randomizer.nextDouble();
                // set new height
                pickedNode.setHeight(newHeight);
                logProposalRatio += Math.log(upper - lower);
            } */

            // deal with the node relationships
            pickedParent.childBranchNumbers.remove(pickedParentBrNr);
            pickedParent.childBranchNumbers.add(otherChildBrNr);
            attachParent.childBranchNumbers.remove(attachBranchNr);
            attachParent.childBranchNumbers.add(pickedParentBrNr);
            pickedNode.childBranchNumbers.remove(otherChildBrNr);
            pickedNode.childBranchNumbers.add(attachBranchNr);
            pickedParent.updateRelationships();
            otherChild.updateRelationships();
            attachParent.updateRelationships();
            attachChild.updateRelationships();
            pickedNode.updateRelationships();
        }

        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        return logProposalRatio;
    }
}
