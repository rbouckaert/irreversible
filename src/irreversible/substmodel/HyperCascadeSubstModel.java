package irreversible.substmodel;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.ComplexSubstitutionModel;

public class HyperCascadeSubstModel extends ComplexSubstitutionModel {
    final public Input<RealParameter> AGRateInput = new Input<RealParameter>("AGRate", "the rate of mutating "
    		+ "dimension 1 from GG-A to GG-G, "
    		+ "dimension 2 from GA-A to GA-G or AG-A to AG-G, "
    		+ "dimension 3 from AA-A to AA-G.", Validate.REQUIRED);
    final public Input<HyperCascadeDataType> dataTypeInput = new Input<>("dataType", "reference to hyper cascade datatype", Validate.REQUIRED);

    RealParameter AGRate;
    protected double[][] unnormalizedQ;
    protected double[][] storedUnnormalizedQ;
    
    // rates[][0] = from state index
    // rates[][1] = to state index
    // rates[][2] = rate index = number of incompatible parents
    private int [][] rates;

    public HyperCascadeSubstModel() {
        ratesInput.setRule(Validate.OPTIONAL);
        frequenciesInput.setRule(Validate.FORBIDDEN);
    }

    
    @Override
    public void initAndValidate() {
    	AGRate = AGRateInput.get();
    	if (AGRate.getDimension() != 3) {
    		throw new IllegalArgumentException("AGRate must be of dimension 3, but is dimension " + AGRate.getDimension());
    	}
    	updateMatrix = true;
    	
    	HyperCascadeDataType dataType = dataTypeInput.get();

        nrOfStates = dataType.getStateCount();
        
        try {
			eigenSystem = createEigenSystem();
		} catch (SecurityException e) {
			throw new IllegalArgumentException(e.getMessage());
		}

        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates * (nrOfStates-1)];
        storedRelativeRates = new double[nrOfStates * (nrOfStates-1)];
    	
        unnormalizedQ = new double[nrOfStates][nrOfStates];
        storedUnnormalizedQ = new double[nrOfStates][nrOfStates];
        initRateIndices(dataType);
    }
	

    private void initRateIndices(HyperCascadeDataType dataType) {
    	List<String> states = new ArrayList<>(nrOfStates);
    	for (int i = 0; i < nrOfStates; i++) {
    		states.add(dataType.getCharacter(i));
    	}
    	
    	List<int[]> compatiblePairs = new ArrayList<>();
    	for (int i = 0; i < nrOfStates; i++) {
    		for (int j = 0; j < nrOfStates; j++) {
    			int rateID = isCompatible(states, i, j);
    			if (rateID >= 0) {
    				compatiblePairs.add(new int[] {i, j, rateID});
    			}
    		}
    	}
    	
    	rates = new int[compatiblePairs.size()][3];
    	for (int i = 0; i < rates.length; i++) {
    		rates[i][0] = compatiblePairs.get(i)[0];
    		rates[i][1] = compatiblePairs.get(i)[1];
    		rates[i][2] = compatiblePairs.get(i)[2];
    	}
	}

    // two states are compatible if 
    // 1. the states differ by only a single 1
    // 2. both parent states are 1 
	private int isCompatible(List<String> states, int i, int j) {
		if (i == j) {
			return -1;
		}
		String state1 = states.get(i);
		String state2 = states.get(j);
		boolean diffBy1 = false;
		
		int layerCount = dataTypeInput.get().layersInput.get();
		int k = 0;
		int position = -1, layer = -1;
		for (int r = 0; r < layerCount; r++) {
			for (int s = 0; s < layerCount-r; s++) {
				char c1 = state1.charAt(k);
				char c2 = state2.charAt(k);
				if (c1 == '1') {
					if (c2 == '0') {
						return -1;
					}
				} else { // c1 == 0
					if (c2 == '1') {
						if (diffBy1) {
							return -1;
						} else {
							diffBy1 = true;
							position = k;
							layer = r;
						}
					}
				}
				k++;
			}
		}

		String state1_ = state1.replaceAll("([01][01][01])([01][01])([01])", "$1\\\\n$2\\\\n$3");
		String state2_ = state2.replaceAll("([01][01][01])([01][01])([01])", "$1\\\\n$2\\\\n$3");
		
		System.out.print("\"" + state1_ + "\" -> \"" + state2_ +"\"");
		if (layer == 0) {
			System.out.print("[color=\"blue\"];\n");
			return 0;
		}
		if (state1.charAt(position -layerCount+layer-1) == '1' && 
			state1.charAt(position -layerCount+layer) == '1') {
			System.out.print("[color=\"blue\"];\n");
			return 0;
		}

		if (state1.charAt(position -layerCount+layer-1) == '1' || 
			state1.charAt(position -layerCount+layer) == '1') {
			System.out.print("[color=\"black\"];\n");
			return 1;
		}
		
		System.out.print("[color=\"grey\"];\n");
		return 2;
	}



	@Override
    protected void setupRelativeRates() {
    }

    @Override
	public void setupRateMatrix() {
        setupUnnormalizedQMatrix();

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = unnormalizedQ[i][j];
            }
        }

        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }
    } // setupRateMatrix

    @Override
    public double[] getFrequencies() {
        double[] fFreqs = new double[4];
        fFreqs[0] = 1;
        return fFreqs;
    }

    protected void setupUnnormalizedQMatrix() {
    	Double [] r = AGRate.getValues();
    	for (int i = 0; i < rates.length; i++) {
	        unnormalizedQ[rates[i][0]][rates[i][1]] = r[rates[i][2]];
    	}
	        
    }
	
	@Override
	public boolean canHandleDataType(DataType dataType) {
	       return dataType instanceof HyperCascadeDataType;
	}
	
	
	
	public static void main(String[] args) {
		HyperCascadeDataType h = new HyperCascadeDataType();
		h.initByName("layers", 3);
		
		
		HyperCascadeSubstModel subst = new HyperCascadeSubstModel();
		subst.initByName("dataType", h, "AGRate", "1.0 1.0 1.0");
	}
}
