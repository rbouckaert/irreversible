package irreversible.substmodel;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;

@Description("Returns P(X|Y) based on P(X,Y) and P(Y) by calculating P(X|Y)=P(X,Y)/P(Y)")
public class ConditionedOnDistribution extends Distribution {
	final public Input<Distribution> distrXYInput = new Input<>("distrXY", "main distribution P(X,Y) with extra condioned on distribution included", Validate.REQUIRED);
	final public Input<Distribution> distrYInput = new Input<>("distrY", "distribution P(Y) to be conditioned on", Validate.REQUIRED);

	
	@Override
	public double calculateLogP() {		
		return distrXYInput.get().calculateLogP() - distrYInput.get().calculateLogP();
	}
	
	@Override
	public double getCurrentLogP() {
		return distrXYInput.get().getCurrentLogP() - distrYInput.get().getCurrentLogP();
	}
	
	
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}
