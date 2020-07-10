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
 * Change origin height using a multiplier with reflection.
 *
 * @author Chi Zhang
 */

@Description("Change origin height using a multiplier  with reflection.")
public class OriginMultiplier extends Operator {
    public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
    public final Input<RealParameter> originInput =
            new Input<>("origin", "The time when the process started.", Validate.REQUIRED);
    public final Input<Double> tuningInput =
            new Input<>("tuning", "A fine-tuning parameter (default is 1).", 1.0);
    public final Input<Boolean> optimiseInput = new Input<>("optimise",
            "Automatic tuning in order to achieve a good acceptance rate (default false)", false);

    private double tuning;

    @Override
    public void initAndValidate() {
        tuning = tuningInput.get();
    }

    @Override
    public double proposal() {
        final Network speciesNetwork = speciesNetworkInput.get();

        // determine the lower and upper bounds
        // double upper = Double.MAX_VALUE;
        double lower = speciesNetwork.getRoot().getHeight();

        // propose a new height, reflect it back if it's outside the boundary
        NetworkNode origin = speciesNetwork.getOrigin();
        final double oldHeight = origin.getHeight();
        double newHeight = oldHeight * Math.exp (tuning * (Randomizer.nextDouble() - 0.5));
        while (newHeight < lower ) {
            newHeight = lower * lower / newHeight;
        }

        // update the origin height
        final RealParameter originTime = originInput.get();
        if (outsideBounds(newHeight, originTime))
            return Double.NEGATIVE_INFINITY;
        originTime.setValue(newHeight);
        speciesNetwork.startEditing(this);
        origin.setHeight(newHeight);
        SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin());

        // log proposal ratio
        return Math.log(newHeight / oldHeight);
    }

    private boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
    }

    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            setCoercableParameterValue(tuning * Math.exp(delta));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return tuning;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        tuning = Math.max(Math.min(value, 20.0), 1e-4);
    }
}
