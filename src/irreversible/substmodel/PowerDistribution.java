package irreversible.substmodel;

import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;

@Description("Distribution taken to a power, thus up or down weighing its impact")
public class PowerDistribution extends Distribution {
	final public Input<Distribution> distributionInput = new Input<>("distribution", "individual probability distributions, e.g. the likelihood and prior making up a posterior", Validate.REQUIRED);
	final public Input<Function> powerInput = new Input<>("power", "fraction used to weight the distribution", Validate.REQUIRED);

	@Override
	public double calculateLogP() {
		logP = distributionInput.get().calculateLogP();
		double power = powerInput.get().getArrayValue();
		logP *= power;
		return logP;
	}
	
	@Override
	public double getCurrentLogP() {
		logP = distributionInput.get().getCurrentLogP();
		double power = powerInput.get().getArrayValue();
		logP *= power;
		return logP;
	}

	
	@Override
	public List<String> getArguments() {
		return null;
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}

}
