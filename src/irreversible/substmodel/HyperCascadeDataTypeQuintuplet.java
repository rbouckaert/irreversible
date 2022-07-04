package irreversible.substmodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.DataType;
import irreversible.util.HyperCascadeCounter;


@Description("Data type for hypercascade data representing states of n sites in n layers")
public class HyperCascadeDataTypeQuintuplet extends DataType.Base {
	public Input<Integer> layersInput = new Input<>("layers", "number of layers in system", 4);

	private double [][]ambiguousPartials;
	
    public HyperCascadeDataTypeQuintuplet() {} // default c'tor
    
    @Override
    public void initAndValidate() {
    	List<String> states = getStates();

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

    private List<String> getStates() {
    	// determine number of states;
    	int layerSize = 5;
    	
    	int [] mem = new int[layerSize];
    	int n = 32;
    	
    	
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
		
		
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				if (isCompatible(code[i], code[j])) {
					ambiguousPartials[i][j] = 1.0;
				}
			}
		}
	}

	private boolean isCompatible(String state1, String state2) {
		return state1.charAt(3) == state2.charAt(3) && state1.charAt(4) == state2.charAt(4);		
	}


}
