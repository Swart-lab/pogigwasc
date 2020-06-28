package de.vetter.masterthesis.states;

import de.vetter.masterthesis.Utilities;

public class NoncodingState extends HMMState {

	// TODO: preferably read this from parameter-file
	private double[] baseFrequencies = new double[] {0.46, 0.10, 0.355, 0.085};
	
	public NoncodingState(String name) { super(name); }

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() > 1) {
			return Double.NEGATIVE_INFINITY;
		}
		return Math.log(baseFrequencies[Utilities.baseToIndex(newEmission.charAt(0))]);
	}

}
