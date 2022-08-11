package irreversible.substmodel;

import java.util.ArrayList;

import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import irreversible.util.HyperCascadeCounter;


@Description("Data type for hypercascade data representing states of n sites in n layers")
public class HyperCascadeDataTypeTriplet extends DataType.Base {

	private double [][]ambiguousPartials;
	
    public HyperCascadeDataTypeTriplet() {} // default c'tor
    
    @Override
    public void initAndValidate() {
    	List<String> states = getStates(2);

    	stateCount = states.size();
        codeLength = states.get(0).length();

		codeMap = "";
		mapCodeToStateSet = new int[stateCount][];
		int k = 0;
		for (String str : states) {
		    // process the code
			String code = str;
			codeMap += code;
			if (codeLength > 0) {
			    if (code.length() != codeLength) {
			        throw new IllegalArgumentException("Invalide code '" + code + "'. Expected code of length " + codeLength);
			    }
			} else {
			    codeMap += ",";
			}
			mapCodeToStateSet[k] = new int[k];
			k++;
		}
    }

    private List<String> getStates(Integer layerCount) {
    	// determine number of states;
    	int layerSize = 0;
    	for (int i = 0; i < layerCount; i++) {
    		layerSize += layerCount-i;
    	}
    	
    	int [] mem = new int[layerSize];
    	int n = 1 << layerSize;
    	
    	
    	// list all states
    	List<String> states = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			Arrays.fill(mem, 0);
			HyperCascadeCounter.toBinary(i, mem, layerSize);
			StringBuilder b = new StringBuilder();
			for (int j = 0; j < layerSize; j++) {
				if (mem[j] == 0) {
					b.append('0');
 				} else {
 					b.append('1');
 				}
			}
			states.add(b.toString());
		}
		return states;
	}

	@Override
	public String getCharacter(int code) {
        if (codeLength > 0) {
        	return codeMap.substring(code * codeLength, (code + 1) * codeLength);
        } else {
        	return String.valueOf(codeMap.split(",")[code]);
		}
    }

    @Override
    public boolean hasConstantCodeLength() {
    	return true;
    }
	
	@Override
	public String getTypeDescription() {
		return "HyperCascadeData";
	}

	public double[] getAmbiguousPartials(int state) {
		if (ambiguousPartials == null) {
			initAmbiguousPartials();
		}
		return ambiguousPartials[state];
	}	
	
	
	private void initAmbiguousPartials() {
		int stateCount = getStateCount();
		ambiguousPartials = new double[stateCount][stateCount];
		
		String [] code = new String[stateCount];
		for (int i = 0; i < stateCount; i++) {
			code[i] = getCharacter(i);
		}
		
		boolean [] isFirstInLayer = new boolean[code[0].length()];
		int layers = 2;
		int k = 0;
		for (int i = 0; i < layers; i++) {
			isFirstInLayer[k] = true;
			k+=layers-i;
		}
		
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				if (isCompatible(code[i], code[j], isFirstInLayer)) {
					ambiguousPartials[i][j] = 1.0;
				}
			}
		}
	}

	private boolean isCompatible(String state1, String state2, boolean [] isFirstInLayer) {
		for (int i = 0; i < state1.length(); i++) {
			if (!(isFirstInLayer[i] || state1.charAt(i) == state2.charAt(i))) {
				return false;
			}
		}
		return true;
	}


}
