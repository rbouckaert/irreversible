package irreversible.substmodel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.DataType;
import irreversible.util.HyperCascadeCounter;


@Description("Data type for hypercascade data representing states of n sites in n layers")
public class HyperCascadeDataType extends DataType.Base {
	public Input<Integer> layersInput = new Input<>("layers", "number of layers in system", 4);

    public HyperCascadeDataType() {} // default c'tor
    
    @Override
    public void initAndValidate() {
    	List<String> states = getStates(layersInput.get());

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
}
