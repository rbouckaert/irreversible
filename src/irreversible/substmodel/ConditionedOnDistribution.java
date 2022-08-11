package irreversible.substmodel;

import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.core.Input.Validate;

@Description("Returns P(X|Y) based on P(X,Y) and P(Y) by calculating P(X|Y)=P(X,Y)/P(Y)")
public class ConditionedOnDistribution extends Distribution {
	final public Input<Distribution> distrXYInput = new Input<>("distrXY", "main distribution P(X,Y) with extra condioned on distribution included", Validate.REQUIRED);
	final public Input<Distribution> distrYInput = new Input<>("distrY", "distribution P(Y) to be conditioned on", Validate.REQUIRED);

	
	@Override
	public double calculateLogP() {		
		final double logPXY = distrXYInput.get().calculateLogP();
		final double logPY  = distrYInput.get().calculateLogP();
		final double logP = logPXY - logPY;
		if (Double.isNaN(logP)) {
//			System.err.println(logPXY + " - " + logPY);
			return Double.NEGATIVE_INFINITY;
		}
		return logP;
	}
	
	@Override
	public double getCurrentLogP() {
		final double logPXY = distrXYInput.get().getCurrentLogP();
		final double logPY  = distrYInput.get().getCurrentLogP();
		return logPXY - logPY;
	}
	
	
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}
