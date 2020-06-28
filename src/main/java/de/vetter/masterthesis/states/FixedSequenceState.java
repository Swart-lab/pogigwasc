package de.vetter.masterthesis.states;

/**
 * States like the start-codon only emit exactly one sequence with probability one.
 * 
 * @author David Emanuel Vetter
 */
public class FixedSequenceState extends HMMState {
	private String sequence;
	
	public FixedSequenceState(String name, String sequence) { 
		super(name); 
		this.sequence = sequence;
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		return newEmission == sequence ? 0 : Double.NEGATIVE_INFINITY;
	}

	
	public String getSequence() {
		return sequence;
	}
	
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
}
