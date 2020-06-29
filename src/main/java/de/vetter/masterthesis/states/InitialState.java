package de.vetter.masterthesis.states;

public class InitialState extends HMMState {

	/** Name of the initial state; reserved for this state alone */
	public static final String INITIAL_STATE_NAME = "initial";

	
	public InitialState() {
		super(INITIAL_STATE_NAME);
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		return newEmission.length() == 0 ? 0 : Double.NEGATIVE_INFINITY;
	}

}
