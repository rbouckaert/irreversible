package test.irreversible.substmodel;

import java.util.*;

import org.junit.jupiter.api.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import irreversible.substmodel.HyperCascadeAlignment;
import irreversible.substmodel.HyperCascadeDataType;

public class HyperCascadeAlignmentTest {
	
	@Test
	public void testHyperCascadeAlignment() {
		Sequence seq = new Sequence("t1", "AGGGAGGAAA");
		List<Sequence> sequences = new ArrayList<>();
		sequences.add(seq);
		
		HyperCascadeDataType dataType = new HyperCascadeDataType();
		dataType.initByName("layers", 4);
		
		Alignment data = new Alignment(sequences, "nucleotide");
		
		HyperCascadeAlignment hyperData = new HyperCascadeAlignment();
		hyperData.initByName("data", data, "userDataType", dataType);
		
		hyperData.getSiteCount();
		int [] pattern = hyperData.getPattern(0);
		System.out.println(Arrays.toString(pattern));
	}

}
