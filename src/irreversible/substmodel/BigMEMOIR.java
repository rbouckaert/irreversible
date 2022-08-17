package irreversible.substmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;

@Description("Irreversible substitution model to handle bigMEMOIR data."
		+ "It assumes a 4 state model, where the first state represents AA "
		+ "and the other three states AG, GA and GG respectively.")
public class BigMEMOIR extends ComplexSubstitutionModel {
	
    public Input<RealParameter> AGRateInput = new Input<RealParameter>("AGRate", "the rate of mutating from AA to AG", Validate.REQUIRED);
    public Input<RealParameter> GARateInput = new Input<RealParameter>("GARate", "the rate of mutating from AA to GA", Validate.REQUIRED);
    public Input<RealParameter> GGRateInput = new Input<RealParameter>("GGRate", "the rate of mutating from AA to GG", Validate.REQUIRED);

    RealParameter AGRate, GARate, GGRate;
    protected double[][] unnormalizedQ;
    protected double[][] storedUnnormalizedQ;

    public BigMEMOIR() {
        ratesInput.setRule(Validate.OPTIONAL);
        frequenciesInput.setRule(Validate.FORBIDDEN);
    }

    
    @Override
    public void initAndValidate() {
    	AGRate = AGRateInput.get();
    	GARate = GARateInput.get();
    	GGRate = GGRateInput.get();

        updateMatrix = true;
        nrOfStates = 4;
        
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
    }
	

    @Override
    public void setupRelativeRates() {
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
	        unnormalizedQ[0][1] = AGRate.getValue();
	        unnormalizedQ[0][2] = GARate.getValue();
	        unnormalizedQ[0][3] = GGRate.getValue();
	
	        unnormalizedQ[1][0] = 0.0;
	
	        unnormalizedQ[1][2] = 0.0;
	        unnormalizedQ[1][3] = 0.0;
	        
	        
	        unnormalizedQ[2][0] = 0.0;
	        unnormalizedQ[2][1] = 0.0;
	
	        unnormalizedQ[2][3] = 0.0;
	
	        unnormalizedQ[3][0] = 0.0;
	        unnormalizedQ[3][1] = 0.0;
	        unnormalizedQ[3][2] = 0.0;
	        
    }
	
	@Override
	public boolean canHandleDataType(DataType dataType) {
	       return dataType.getStateCount() == 4;
	}
}
