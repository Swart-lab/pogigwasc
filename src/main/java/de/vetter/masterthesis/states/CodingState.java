package de.vetter.masterthesis.states;

import de.vetter.masterthesis.Utilities;

public class CodingState extends HMMState {
	private boolean strand;

	/**
	 * 
	 * @param name
	 * @param strand whether on forward strand
	 */
	public CodingState(String name, boolean strand) {
		super(name);
		this.strand = strand;
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		// TODO Codon usage or some approximation (maybe always use empirical distr. of base 1, to avoid overfitting?)
		if(newEmission.length() != 3)
			return Double.NEGATIVE_INFINITY;
		
		if(!strand)
			newEmission = Utilities.reverseComplement(newEmission);
		
		double result = 0;
		for(int i = 0; i < 3; i++)
			result += Utilities.getLogBaseProbabilityCDS(newEmission.charAt(i), i);
		
		return result;
	}

}
