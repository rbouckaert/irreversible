package irreversible.branchratemodel;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel.Base;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

@Description("Clock model with decreasing rate due to running out of unmutated sites")
public class DecreasingRateClockModel extends Base {
	public enum interpolation {linear, midexp, integratedexp}
    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree");
    final public Input<Function> endRateInput = new Input<>("endRate", "if specified, use this value as target rate, unless data is specified, then endRate is ignored");
    
    final public Input<interpolation> interpolationInput = new Input<>("interpolation", "interpolation mode: "
    		+ "linear = use linear interpolation (where midpoint equals integral) "
    		+ "midexp = use exponential interpolation at middle of a branch "
    		+ "integratedexp = use exponetial by integrating over a branch", interpolation.midexp, interpolation.values());

    private TreeInterface tree;
    private Function meanRate;
    private double endRate, logU;
    private interpolation mode;
    private Function endRateFunction;
    
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		meanRate = meanRateInput.get();
		endRateFunction = endRateInput.get();
		
		if (dataInput.get() != null) {
			// determine fraction of not mutated sites
			Frequencies f  = new Frequencies();
			f.initByName("data", dataInput.get());
			double [] freqs = f.getFreqs();

			// assuming encoding A=AA C=AG G=GA T=GG
			endRate = freqs[0] + 0.5*freqs[1]+0.5*freqs[2];
			endRate = freqs[0];
		} else {
			endRate = endRateInput.get().getArrayValue();
		}
		logU = Math.log(endRate);
		
		Log.warning("DecreasingRateClockModel::endRate = " + endRate);
		
		mode = interpolationInput.get();
	}
    
    
	@Override
	public double getRateForBranch(Node node) {
		if (node.isRoot()) {
			return meanRate.getArrayValue();
		}
		if (dataInput.get() != null) {
			endRate = endRateFunction.getArrayValue();
		}		
		double rootHeight = tree.getRoot().getHeight();
		double parentHeight = node.getParent().getHeight();
		double nodeHeight = node.getHeight();
		double mutationFraction = 0;
		switch (mode) {
		case linear:
			// linear interpolation between root where mutation fraction = 1
			// and leaf, where mutationFraction = unmutatedFraction
			mutationFraction = endRate + ((parentHeight + nodeHeight)/2.0)/rootHeight * (1-endRate); 
		break;
		case midexp:
			// use exponential interpolation instead
			mutationFraction = Math.exp(logU * (1.0-(((parentHeight + nodeHeight)/2.0)/rootHeight)));
			mutationFraction = (f(parentHeight, rootHeight) +
			 f((7*parentHeight + nodeHeight)/8.0, rootHeight) +
			 f((3*parentHeight + nodeHeight)/4.0, rootHeight) +
			 f((5*parentHeight + 3*nodeHeight)/8.0, rootHeight) +
			 f((parentHeight + nodeHeight)/2.0, rootHeight) +
			 f((3*parentHeight + 5*nodeHeight)/8.0, rootHeight) +
			 f((parentHeight + 3*nodeHeight)/4.0, rootHeight) +
			 f((parentHeight + 7*nodeHeight)/8.0, rootHeight) +
			 f(nodeHeight, rootHeight))/9.0;
		break;
		case integratedexp:
			// integrate over exponential
			// integral_a^b exp(c x) dx = (e^(b c) - e^(a c))/c
			// with 
			// a = node height
			// b = parent height
			// integral_a^b exp(logU * (1.0-x/rootHeight)) 
			// = integral_a^b exp(logU - x*logU/rootHeight) 
			// = exp(logU)* integral_a^b exp(-x*logU/rootHeight)
			// so c = -logU/rootHeight
			
			 double c = -logU/rootHeight;
			 mutationFraction = endRate * (Math.exp(parentHeight *c) - Math.exp(nodeHeight *c))/c;
			 mutationFraction /= (parentHeight - nodeHeight);
			 break;
		}
		return mutationFraction * meanRate.getArrayValue();
	}


	private double f(double nodeHeight, double rootHeight) {
		double r = Math.exp(logU * (1.0-nodeHeight/rootHeight)); 
		return r;
	}


}
