package irreversible.substmodel;


import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;

@Description("Class representing hyper cascade alignment data")
public class HyperCascadeAlignment extends Alignment {
    final public Input<Alignment> alignmentInput = new Input<>("data", "alignment to be filtered", Validate.REQUIRED);

	
	@Override
	public void initAndValidate() {
        Alignment data = alignmentInput.get();
        m_dataType = userDataTypeInput.get();
        if (!(m_dataType instanceof HyperCascadeDataType)) {
        	throw new IllegalArgumentException("Expected user data type of type " + HyperCascadeDataType.class.getName());
        }
        
        taxaNames = data.getTaxaNames();
        int taxonCount = taxaNames.size();

        int layerCount = ((HyperCascadeDataType)m_dataType).layersInput.get();
        
        int siteCount = data.getSiteCount();
        for (int i = 1; i < layerCount; i++) {
        	siteCount += i;
        }
        siteCount /= layerCount;
        
        
        stateCounts = new ArrayList<>();
        maxStateCount = m_dataType.getStateCount();
    	for (int i = 0; i < siteCount; i++) {
            stateCounts.set(i, maxStateCount);
    	}

		siteWeights = new int[siteCount];
		Arrays.fill(siteWeights, 1);

		patternWeight = new int[siteCount];
		Arrays.fill(patternWeight, 1);
		
		patternIndex = new int[siteCount];
		for (int i = 0; i < siteCount; i++) {
			patternIndex[i] = i;
		}

		sitePatterns = new int[siteCount][taxonCount];
		for (int j = 0; j < taxonCount; j++) {
			List<Integer> states = calcSiteValue(siteCount, j, layerCount, data);
			for (int i = 0; i < siteCount; i++) {
				sitePatterns[i][j] = states.get(i);
			}
		}
        // calcPatterns();
	}


	private List<Integer> calcSiteValue(int siteCount, int taxon, int layerCount, Alignment data) {
		String seq = data.getSequenceAsString(taxaNames.get(taxon)).toUpperCase();
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < siteCount; i++) {
			int offset = i;
			for (int j = 0; j < layerCount; j++) {
				for (int k = 0; k < layerCount - j; k++) {
					char c = seq.charAt(offset + k);
					if (c == 'A') {
						b.append('0');
					} else {
						b.append('1');
					}						
				}
				b.append(' ');
				offset += siteCount - j;
			}
		}
		
		return m_dataType.stringToEncoding(b.toString());
	}
	
	
	
	
	
}
