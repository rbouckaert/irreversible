package irreversible.substmodel;


import beast.base.core.Description;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.ComplexColtEigenSystem;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;

@Description("Irreversible pure birth substitution model")
public class PureBirthModel extends ComplexSubstitutionModel {

	
    public PureBirthModel() {
    	ratesInput.setRule(Validate.OPTIONAL);
    }
    

    @Override
    public void initAndValidate() {
        frequencies = frequenciesInput.get();
        double[] freqs = getFrequencies();
        nrOfStates = freqs.length;
        rateMatrix = new double[(nrOfStates)][(nrOfStates)];
        
        eigenSystem = new ComplexColtEigenSystem(nrOfStates);
    }


	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
        synchronized (this) {
            if (updateMatrix) {
            	double [][] rateMatrix = getRateMatrix();

                eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
                updateMatrix = false;
            }
        }
        return eigenDecomposition;
    }
    
	public void setupRelativeRates() {
		if (relativeRates == null) {
			relativeRates = new double[2];
		}
    	int k = 0;
		relativeRates[k++] = 1;
		relativeRates[k++] = 0;
    }

	public double[][] getRateMatrix() {
		double [][] rateMatrix = new double[2][2];
		rateMatrix[0][0] = -1;
		rateMatrix[0][1] = 1; 
		rateMatrix[1][0] = 0;
		rateMatrix[1][1] = 0;
		return rateMatrix;
	}

	@Override
	public void setupRateMatrix() {
		rateMatrix[0][0] = -1;
		rateMatrix[0][1] = 1; 
		rateMatrix[1][0] = 0;
		rateMatrix[1][1] = 0;
	}
	
	@Override
	public void restore() {
		updateMatrix = true;
		super.restore();
	}
	
    @Override
    protected boolean requiresRecalculation() {
        // we only get here if delParameter or mutationRate is dirty
    	updateMatrix = true;
        return true;
    }
    
    
    @Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType.getStateCount() == 2;
	}


}
