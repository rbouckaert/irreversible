package test.irreversible.substmodel;

import java.io.File;
import java.io.IOException;
import java.util.*;

import org.junit.jupiter.api.Test;

import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.parser.NexusParser;
import irreversible.substmodel.HyperCascadeAlignmentForTriplets;
import irreversible.substmodel.HyperCascadeDataTypeQuintuplet;

public class HyperCascadeAlignmentTest {
	
	@Test
	public void testHyperCascadeAlignment() throws IOException {
		NexusParser parser = new NexusParser();
		parser.parseFile(new File("examples/hypercascade_hiPSC_barcs_subset.nex"));
		Alignment nucleotideData = parser.m_alignment;
		StringBuilder b = new StringBuilder();
		for (Sequence seq : nucleotideData.sequenceInput.get()) {
			try {
				parseSeq(seq);
			} catch (Throwable e) {
				b.append("parsing failed for taxon " + seq.getTaxon() + " : " + e.getMessage() + "\n");
			}
		}
		Log.warning(b.toString());
		
//		Sequence seq = new Sequence("t1", "AGGGAAAGGG AGGAAAAGG AAAAAAAG AAAAAAA");
//		seq = new Sequence("t1", "GGGAGGGGGGGAAAAAAGA GGAGAAAGGGGAAAAAAA GAAAAAGAGGAAAAAAA AAAAAAGAAAAAAAAA");
//		parseSeq(seq);
				
	}

	private void parseSeq(Sequence seq) {
		List<Sequence> sequences = new ArrayList<>();
		sequences.add(seq);
		
		HyperCascadeDataTypeQuintuplet dataType = new HyperCascadeDataTypeQuintuplet();
		dataType.initByName("layers", 4);
		
		Alignment data = new Alignment(sequences, "nucleotide");
		
		HyperCascadeAlignmentForTriplets hyperData = new HyperCascadeAlignmentForTriplets();
		hyperData.initByName("data", data, "userDataType", dataType);
		
		hyperData.getSiteCount();
		for (int i = 0; i < hyperData.getPatternCount(); i++) {
			int [] pattern = hyperData.getPattern(i);		
			System.out.println(Arrays.toString(pattern) + " " + dataType.getCharacter(pattern[0]));
		}
	}

}
