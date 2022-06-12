package irreversible.substmodel;

import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.ComplexSubstitutionModel;

public class HyperCascadeSubstModel extends ComplexSubstitutionModel {
    final public Input<RealParameter> AGRateInput = new Input<RealParameter>("AGRate", "the rate of mutating from AA to AG", Validate.REQUIRED);
    final public Input<HyperCascadeDataType> dataTypeInput = new Input<>("dataType", "reference to hyper cascade datatype", Validate.REQUIRED);

    RealParameter AGRate;
    protected double[][] unnormalizedQ;
    protected double[][] storedUnnormalizedQ;
    private int [][] rates;

    public HyperCascadeSubstModel() {
        ratesInput.setRule(Validate.OPTIONAL);
        frequenciesInput.setRule(Validate.FORBIDDEN);
    }

    
    @Override
    public void initAndValidate() {
    	AGRate = AGRateInput.get();
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
    			if (isCompatible(states, i, j)) {
    				compatiblePairs.add(new int[] {i, j});
    			}
    		}
    	}
    	
    	rates = new int[compatiblePairs.size()][2];
    	for (int i = 0; i < rates.length; i++) {
    		rates[i][0] = compatiblePairs.get(i)[0];
    		rates[i][1] = compatiblePairs.get(i)[1];
    	}
	}


	private boolean isCompatible(List<String> states, int i, int j) {
		if (i == j) {
			return false;
		}
		final String state1 = states.get(i);
		final String state2 = states.get(j);
		boolean diffBy1 = false;
		for (int k = 0; k < state1.length(); k++) {
			char c1 = state1.charAt(k);
			char c2 = state2.charAt(k);
			if (c1 == '1') {
				if (c2 == '0') {
				return false;
				}
			} else {
				if (c2 == '1') {
					if (diffBy1) {
						return false;
					} else {
						diffBy1 = true;
					}
				}
			}
		}
		System.out.println("\"" + state1 + "\" -> \"" + state2 +"\";");
		return true;
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
    	for (int i = 0; i < rates.length; i++) {
	        unnormalizedQ[rates[i][0]][rates[i][1]] = AGRate.getValue();
    	}
	        
    }
	
	@Override
	public boolean canHandleDataType(DataType dataType) {
	       return dataType instanceof HyperCascadeDataType;
	}
	
	
	
	public static void main(String[] args) {
		HyperCascadeDataType h = new HyperCascadeDataType();
		h.initByName("layers", 4);
		
		
		HyperCascadeSubstModel subst = new HyperCascadeSubstModel();
		subst.initByName("dataType", h, "AGRate", "1.0");
	}
}
