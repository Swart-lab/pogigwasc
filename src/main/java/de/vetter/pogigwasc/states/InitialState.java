package de.vetter.pogigwasc.states;

/**
 * The special initial state (cf. Stanke's dissertation): Can only emit the
 * empty word. It has the reserved name {@link #INITIAL_STATE_NAME}
 * 
 * @author David Emanuel Vetter
 */
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
