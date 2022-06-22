package irreversible.substmodel;


import java.util.*;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;

@Description("Class representing hyper cascade alignment data")
public class HyperCascadeAlignment extends Alignment {
    final public Input<Alignment> alignmentInput = new Input<>("data", "alignment to be filtered", Validate.REQUIRED);
    final public Input<Integer> offsetInput = new Input<>("offset", "offset of the alignment to be filtered = number of "
    		+ "sites at start and end to be removed. By default, all sites will be included.", 0);
	public Input<Boolean> addAmbiguitiesInput = new Input<>("addAmbiguities", "assume first entry in each layer is "
			+ "ambiguous. This can be useful for conditioning on.", false);

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
        
        int firstLayerSiteCount = data.getSiteCount();
        for (int i = 1; i < layerCount; i++) {
        	firstLayerSiteCount += i;
        }
        firstLayerSiteCount /= layerCount;
        
        int offset = offsetInput.get();
        firstLayerSiteCount -= 2 * offset;
        
        stateCounts = new ArrayList<>();
        maxStateCount = m_dataType.getStateCount();
    	for (int i = 0; i < firstLayerSiteCount; i++) {
            stateCounts.add(maxStateCount);
    	}

		siteWeights = new int[firstLayerSiteCount];
		Arrays.fill(siteWeights, 1);

		patternWeight = new int[firstLayerSiteCount];
		Arrays.fill(patternWeight, 1);
		
		patternIndex = new int[firstLayerSiteCount];
		for (int i = 0; i < firstLayerSiteCount; i++) {
			patternIndex[i] = i;
		}

		sitePatterns = new int[firstLayerSiteCount- layerCount + 1][taxonCount];
		for (int j = 0; j < taxonCount; j++) {
			List<Integer> states = calcSiteValue(firstLayerSiteCount + 2 * offset, j, layerCount, data);
			for (int i = 0; i < firstLayerSiteCount - layerCount + 1; i++) {
				sitePatterns[i][j] = states.get(i + offset);
			}
		}
        // calcPatterns();
	}


	private List<Integer> calcSiteValue(int siteCount, int taxon, int layerCount, Alignment data) {
		String seq = data.getSequenceAsString(taxaNames.get(taxon)).toUpperCase();
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < siteCount - layerCount + 1; i++) {
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
		if (addAmbiguitiesInput.get()) {
			return ((HyperCascadeDataType)m_dataType).getAmbiguousPartials(state);
		} else {
	         double [] partials = new double[stateCount];
	         partials[state] = 1;
	         return partials;
		}
	}
}
