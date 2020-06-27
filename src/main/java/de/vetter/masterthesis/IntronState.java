package de.vetter.masterthesis;

/**
 * Implements a state for introns, which only considers the emmission-length
 * when computing the emmission-probability (at this point)
 * 
 * @author David Emanuel Vetter
 */
public class IntronState extends HMMState {

	public IntronState(String name) {
		super(name);
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		// TODO Auto-generated method stub: Implement
		// -> Use poisson.
		return 0;
	}

}
