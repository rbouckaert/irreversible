package irreversible.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.MRCAPrior;

@Description("Prior on clade being present")
public class CladePrior extends MRCAPrior {
	final public Input<Double> cladeProbInput = new Input<>("cladeProb", "probability that the clade is present in a tree", 0.5);

	
	private double logCladeProb;
	private double logNotCladeProb;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		logCladeProb = Math.log(cladeProbInput.get());
		logNotCladeProb = Math.log(1.0 - cladeProbInput.get());
	}
	
	
	@Override
	public double calculateLogP() {
		super.calculateLogP();
		logP = isMonophyletic ? logCladeProb : logNotCladeProb;
		return logP;
	}
}
