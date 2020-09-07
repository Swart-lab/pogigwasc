package de.vetter.pogigwasc.states;

public class TerminalState extends HMMState {
	
	/** Name of the terminal state; reserved for this state alone */
	public static final String TERMINAL_STATE_NAME = "terminal";

	public TerminalState() {
		super(TERMINAL_STATE_NAME);
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		return newEmission.length() == 0 ? 0 : Double.NEGATIVE_INFINITY;
	}

}
