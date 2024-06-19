package irreversible.substmodel;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.ComplexColtEigenSystem;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;

@Description("HyperCascadeSubstModel that allows A->M->G and A->D transitions where"
		+ "D = deleted, "
		+ "M = merged between A and G. "
		+ "States: A=0, M=1, G=2, D=3. "
		+ "Rate[0] = A->M, Rate[1] = M->G, Rate[2] = A->D")
public class HyperCascadeSubstModelWithDeletions extends ComplexSubstitutionModel {

	private Function rates;
	
	@Override
    public void initAndValidate() {
        frequencies = frequenciesInput.get();
        double[] freqs = getFrequencies();
        nrOfStates = freqs.length;
        if (nrOfStates != 4) {
        	throw new IllegalArgumentException("Expected 4 states in frequencies");
        }
        rateMatrix = new double[(nrOfStates)][(nrOfStates)];
        
        eigenSystem = new ComplexColtEigenSystem(nrOfStates);
        rates = ratesInput.get();
        if (rates.getDimension() != 3) {
        	throw new IllegalArgumentException("Expected 3 dimensions for rates");
        }
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
			relativeRates = new double[12];
		}
		relativeRates[0] = rates.getArrayValue(0);
		relativeRates[4] = rates.getArrayValue(1);
		relativeRates[2] = rates.getArrayValue(2);
    }

	public double[][] getRateMatrix() {
		double [][] rateMatrix = new double[4][4];
		rateMatrix[0][0] = -rates.getArrayValue(0) - rates.getArrayValue(2);
		rateMatrix[0][1] = rates.getArrayValue(0); 
		rateMatrix[0][3] = rates.getArrayValue(2);
		rateMatrix[1][1] = -rates.getArrayValue(1);
		rateMatrix[1][2] = rates.getArrayValue(1);
		return rateMatrix;
	}

	@Override
	public void setupRateMatrix() {
		rateMatrix[0][0] = -rates.getArrayValue(0) - rates.getArrayValue(2);
		rateMatrix[0][1] = rates.getArrayValue(0); 
		rateMatrix[0][3] = rates.getArrayValue(2);
		rateMatrix[1][1] = -rates.getArrayValue(1);
		rateMatrix[1][2] = rates.getArrayValue(1);
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
		return dataType.getStateCount() == 4;
	}


}
