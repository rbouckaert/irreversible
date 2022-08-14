package test.irreversible.branchratemodel;

import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import irreversible.branchratemodel.DecreasingRateClockModel;
import test.beast.BEASTTestCase;

public class DecreasingRateClockModelTest {
	
	double endRate = 0.0666;
	
	@Test
	public void testRates() throws Exception {
		Alignment data = BEASTTestCase.getAlignment();
		Tree tree = BEASTTestCase.getTree(data);
		
		
		DecreasingRateClockModel clockModel = new DecreasingRateClockModel();
		clockModel.initByName("tree", tree, 
				"interpolation", "linear", 
				"endRate", endRate, 
				"clock.rate", "1.0");
		printRates(clockModel, tree);
		
		clockModel.initByName("tree", tree, 
				"interpolation", "midexp", 
				"endRate", endRate, 
				"clock.rate", "1.0");
		printRates(clockModel, tree);

		clockModel.initByName("tree", tree, 
				"interpolation", "integratedexp", 
				"endRate", endRate, 
				"clock.rate", "1.0");
		printRates(clockModel, tree);

		
	}

	private void printRates(DecreasingRateClockModel clockModel, Tree tree) {
		System.out.println("mode: " + clockModel.interpolationInput.get());
		System.out.println("rate                absolute-error");
		for (Node node : tree.getNodesAsArray()) {
			double target = f(node);
			double calculated = clockModel.getRateForBranch(node);
			System.out.println(calculated + " " + (calculated - target));
		}
	}

	private double f(Node node) {
		if (node.isRoot()) {
			return 1;
		}
		
		// numerical integration over interval
		double rootHeight = node.getTree().getRoot().getHeight();
		double parentHeight = node.getParent().getHeight();
		double nodeHeight = node.getHeight();
		int N = 1000000;
		double sum = 0;
		for (int i = 0; i <= N; i++) {
			sum += f((i*parentHeight + (N-i)*nodeHeight)/(double)N, rootHeight);
		}
		sum /= (N+1);
		return sum;
	}

	private double f(double nodeHeight, double rootHeight) {
		double logU = Math.log(endRate);
		double r = Math.exp(logU * (1.0-nodeHeight/rootHeight)); 
		return r;
	}

}
