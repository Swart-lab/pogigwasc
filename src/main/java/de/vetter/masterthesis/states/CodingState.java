package de.vetter.masterthesis.states;

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
		
		
		return 0;
	}

}
