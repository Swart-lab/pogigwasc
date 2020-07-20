package de.vetter.masterthesis.states;

import de.vetter.masterthesis.ModelParameters;

public abstract class HMMStateWithStrandAndParameters extends HMMState {

	private boolean strand;
	protected ModelParameters parameters;
	
	public HMMStateWithStrandAndParameters(String name, boolean strand, ModelParameters parameters) {
		super(name);
		this.strand = strand;
		this.parameters = parameters;
	}

	/**
	 * @return
	 */
	public boolean isForward() {
		return strand;
	}
	
	/**
	 * @return whether this state is concerned with the reverse strand: Method introduced for readability
	 */
	public boolean isReverse() {
		return !strand;
	}

}
