package irreversible.util;


import java.io.File;
import java.io.PrintStream;
import java.util.List;

import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Runnable;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.TaxonSet;
import beastlabs.evolution.tree.RNNIMetric;
import beastlabs.evolution.tree.RobinsonsFouldMetric;
import beast.base.evolution.tree.ClusterTree;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeMetric;
import beast.base.parser.NexusParser;



@Description("Fits UPGMA tree to an alignemnt and ")
public class UPGMAFitter extends Runnable {
	final public Input<String> dirInput = new Input<>("dir", "directory  containing the alignments", "/tmp/hc/");
	final public Input<String> referenceNewickInput = new Input<>("reference", "file containing reference trees in NEXUS format", "/tmp/hc/truth.trees");
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("[[none]]"));
	final public Input<TreeFile> treeOutputInput = new Input<>("treeout", "output file for UPGMA treees, or stdout if not specified",
			new TreeFile("[[none]]"));
	final public Input<Double> rootHeightInput = new Input<>("rootHeight", "if specified, scale tree to get desired root height");
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
			out = new PrintStream(outputInput.get());
		}

		PrintStream outtree = System.out;
		if (treeOutputInput.get() != null && !treeOutputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + treeOutputInput.get().getPath());
			outtree = new PrintStream(treeOutputInput.get());
		}

		NexusParser parser = new NexusParser();
		parser.parseFile(new File(referenceNewickInput.get()));
		List<Tree> trees = parser.trees;
		TaxonSet taxonset = trees.get(0).getTaxonset();
		
				
		
		String dir = dirInput.get();
		out.println("sample\tRNNIDistance\tRFDistance");
		double sumRF = 0, sumRNNI = 0;
		for (int i = 0; i < 100; i++) {
			String fileName = dir + "Tree" + (i+2) + "_Barcodes.nex";
			NexusParser p2 = new NexusParser();
			p2.parseFile(new File(fileName));
			Alignment data = p2.m_alignment;
			ClusterTree upgmaTree = new ClusterTree();
			upgmaTree.initByName("taxa", data, "clusterType", "irreversible");
			
			TreeMetric metric = new RNNIMetric(taxonset);
			double RNNIDistance = metric.distance(trees.get(i+1), upgmaTree);
			TreeMetric metric2 = new RobinsonsFouldMetric(taxonset);
			double RFDistance = metric2.distance(trees.get(i+1), upgmaTree);
			out.println(i+"\t" + RNNIDistance + "\t" + RFDistance);
			
			if (rootHeightInput.get() != null) {
				double scale = rootHeightInput.get() / upgmaTree.getRoot().getHeight();
				upgmaTree.scale(scale);
			}
			outtree.println(upgmaTree.getRoot().toNewick() + ";");
			sumRF += RFDistance;
			sumRNNI += RNNIDistance;
		}
		
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}
		if (treeOutputInput.get() != null && !treeOutputInput.get().getName().equals("[[none]]")) {
			outtree.close();
		}
		Log.warning(sumRNNI/100 + " " + sumRF/100);
		Log.warning("Done");
	}

	public static void main(String[] args) throws Exception {
		new Application(new UPGMAFitter(), "UPGMA Fitter", args);
	}

}
