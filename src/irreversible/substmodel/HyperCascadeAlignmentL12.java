package irreversible.substmodel;



import java.util.*;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import irreversible.util.HyperCascadeCounter;

@Description("Class representing hyper cascade alignment data for layers 1 & 2. "
		+ "It takes offset=1 and ambiguities over a_{i-1}b_{i-1}.")
public class HyperCascadeAlignmentL12 extends Alignment {
    final public Input<Alignment> alignmentInput = new Input<>("data", "alignment to be filtered", Validate.REQUIRED);
    final public Input<Integer> totalLayersInput = new Input<>("totalLayers", "total number of layers in input alignment.", Validate.REQUIRED);

	private double [][]ambiguousPartials;

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
        
        int layerCount = 2;//((HyperCascadeDataTypeTriplet)m_dataType).layersInput.get();
        int siteCount = firstLayerSiteCount - 2;
        
        
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
		for (int j = 0; j < taxonCount; j++) {
			List<Integer> states = calcSiteValue(firstLayerSiteCount, j, layerCount, data);
			for (int i = 0; i < siteCount; i++) {
				sitePatterns[i][j] = states.get(i);
			}
		}
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
		return ((HyperCascadeDataTypeTriplet)m_dataType).getAmbiguousPartials(state);
	}
	
	public double[] getAmbiguousPartials(int state) {
		if (ambiguousPartials == null) {
			initAmbiguousPartials();
		}
		return ambiguousPartials[state];
	}	
	
	
	private void initAmbiguousPartials() {
		int stateCount = 8;
		ambiguousPartials = new double[stateCount][stateCount];
		
		String [] code = new String[stateCount];
		int [] buf = new int[3];
		for (int i = 0; i < stateCount; i++) {
			HyperCascadeCounter.toBinary(i, buf, 3);
			code[i] = (buf[0] == 0 ? "0" : "1") +
					  (buf[1] == 0 ? "0" : "1") +
					  (buf[2] == 0 ? "0" : "1");
		}
		
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++) {
				if (code[i].charAt(1) == code[j].charAt(1)) {
					ambiguousPartials[i][j] = 1.0;
				}
			}
		}
	}
	
}
