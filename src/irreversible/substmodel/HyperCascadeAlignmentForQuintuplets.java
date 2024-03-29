package irreversible.substmodel;


import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import irreversible.util.HyperCascadeCounter;

@Description("Class representing hyper cascade alignment data for triplets a_{i}a_{i+1}a_{i+2}b_{i}b_{i+1} "
		+ "with ambiguities for a_{i}a_{i+1}a_{i+2}.")
public class HyperCascadeAlignmentForQuintuplets extends Alignment {
    final public Input<Alignment> alignmentInput = new Input<>("data", "alignment to be filtered", Validate.REQUIRED);
    final public Input<Integer> totalLayersInput = new Input<>("totalLayers", "total number of layers in input alignment.", Validate.REQUIRED);

	private double [][]ambiguousPartials;

	@Override
	public void initAndValidate() {
        Alignment data = alignmentInput.get();
        m_dataType = userDataTypeInput.get();
        if (!(m_dataType instanceof HyperCascadeDataTypeQuintuplet)) {
        	throw new IllegalArgumentException("Expected user data type of type " + HyperCascadeDataTypeQuintuplet.class.getName());
        }
        
        taxaNames = data.getTaxaNames();
        int taxonCount = taxaNames.size();

        int totalLayerCount = ((HyperCascadeDataTypeQuintuplet)m_dataType).layersInput.get();
        if (totalLayersInput.get() > 0) {
        	totalLayerCount = totalLayersInput.get();
        }
        
        int firstLayerSiteCount = data.getSiteCount();
        for (int i = 1; i < totalLayerCount; i++) {
        	firstLayerSiteCount += i;
        }
        firstLayerSiteCount /= totalLayerCount;
        
        int layerCount = ((HyperCascadeDataTypeQuintuplet)m_dataType).layersInput.get();
        int siteCount = 0;
        for (int k = 0; k < totalLayerCount - 2; k++) {
        	siteCount += firstLayerSiteCount - k - 2;
        }
        
        
        stateCounts = new ArrayList<>();
        maxStateCount = m_dataType.getStateCount();
    	for (int i = 0; i < siteCount; i++) {
            stateCounts.add(maxStateCount);
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
		int layerOffset = 0;
        for (int k = 0; k < totalLayerCount - 2; k++) {
			for (int j = 0; j < taxonCount; j++) {
				List<Integer> states = calcSiteValue(layerOffset, firstLayerSiteCount - k, j, 3, data);
				for (int i = 0; i < firstLayerSiteCount - k - 2; i++) {
					sitePatterns[i + layerOffset - k * 2][j] = states.get(i);
				}
			}
			layerOffset += firstLayerSiteCount  - k;
        }
        // calcPatterns();
	}


	private List<Integer> calcSiteValue(int layerOffset, int siteCount, int taxon, int layerCount, Alignment data) {
		String seq = data.getSequenceAsString(taxaNames.get(taxon)).toUpperCase();
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < siteCount - layerCount + 1; i++) {
			int offset = i + layerOffset;
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < layerCount - j; k++) {
					char c = seq.charAt(offset + k);
					if (c == 'A') {
						b.append('0');
					} else {
						b.append('1');
					}						
				}
				if (j<layerCount-1) {
					// b.append('_');
					offset += siteCount - j;
				}
			}
		}
		
		return m_dataType.stringToEncoding(b.toString());
	}
	
	
		
	@Override
	public double[] getTipLikelihoods(int taxonIndex, int patternIndex_) {
    	int state = getPattern(taxonIndex, patternIndex_);
		return ((HyperCascadeDataTypeQuintuplet)m_dataType).getAmbiguousPartials(state);
	}
	
	public double[] getAmbiguousPartials(int state) {
		if (ambiguousPartials == null) {
			initAmbiguousPartials();
		}
		return ambiguousPartials[state];
	}	
	
	
	private void initAmbiguousPartials() {
		int stateCount = 32;
		ambiguousPartials = new double[stateCount][stateCount];
		
		String [] code = new String[stateCount];
		int [] buf = new int[5];
		for (int i = 0; i < stateCount; i++) {
			HyperCascadeCounter.toBinary(i, buf, 5);
			code[i] = (buf[0] == 0 ? "0" : "1") +
					  (buf[1] == 0 ? "0" : "1") +
					  (buf[2] == 0 ? "0" : "1") +
					  (buf[3] == 0 ? "0" : "1") +
					  (buf[4] == 0 ? "0" : "1")
					  ;
		}
		
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				if (code[i].charAt(3) == code[j].charAt(3) && code[i].charAt(4) == code[j].charAt(4)) {
					ambiguousPartials[i][j] = 1.0;
				}
			}
		}
	}
	

}
