package irreversible.substmodel;


import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;

@Description("Class representing hyper cascade alignment data for triplets a_{i}a_{i+1}b_{i}.")
public class HyperCascadeAlignmentForTriplets extends Alignment {
    final public Input<Alignment> alignmentInput = new Input<>("data", "alignment to be filtered", Validate.REQUIRED);
    final public Input<Integer> totalLayersInput = new Input<>("totalLayers", "total number of layers in input alignment.", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
        Alignment data = alignmentInput.get();
        m_dataType = userDataTypeInput.get();
        if (!(m_dataType instanceof HyperCascadeDataTypeTriplet)) {
        	throw new IllegalArgumentException("Expected user data type of type " + HyperCascadeDataTypeTriplet.class.getName());
        }
        
        taxaNames = data.getTaxaNames();
        int taxonCount = taxaNames.size();

        int totalLayerCount = totalLayersInput.get();
        
        int firstLayerSiteCount = data.getSiteCount();
        for (int i = 1; i < totalLayerCount; i++) {
        	firstLayerSiteCount += i;
        }
        firstLayerSiteCount /= totalLayerCount;
        

        int layerCount = 2;
        int siteCount = 0;
        for (int k = 0; k <= totalLayerCount - layerCount; k++) {
        	siteCount += firstLayerSiteCount - k - 1;
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
        for (int k = 0; k < totalLayerCount - 1; k++) {
			for (int j = 0; j < taxonCount; j++) {
				List<Integer> states = calcSiteValue(layerOffset, firstLayerSiteCount - k, j, layerCount, data);
				for (int i = 0; i < firstLayerSiteCount - k - 1; i++) {
					sitePatterns[i + layerOffset][j] = states.get(i);
				}
			}
			layerOffset += firstLayerSiteCount - k - 1;
        }
        // calcPatterns();
	}


	private List<Integer> calcSiteValue(int layerOffset, int siteCount, int taxon, int layerCount, Alignment data) {
		String seq = data.getSequenceAsString(taxaNames.get(taxon)).toUpperCase();
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < siteCount - layerCount + 1; i++) {
			int offset = i + layerOffset;
			for (int j = 0; j < layerCount; j++) {
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
    	int stateCount = m_dataType.getStateCount();
         double [] partials = new double[stateCount];
         partials[state] = 1;
         return partials;
	}
}
