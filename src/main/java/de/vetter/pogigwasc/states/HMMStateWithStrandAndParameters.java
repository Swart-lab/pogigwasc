package de.vetter.pogigwasc.states;

import de.vetter.pogigwasc.ModelParameters;

public abstract class HMMStateWithStrandAndParameters extends HMMState {

	private boolean strand;
	protected ModelParameters parameters;
	
	/**
	 * 
	 * @param name       this state's name: Can and should be informative, e.g.
	 *                   corresponding to the state-names used in an external graph
	 *                   of the HMM
	 * @param strand     whether on the forward strand (i.e.: whether working on the
	 *                   sequence/emissions directly, or on reverse complement
	 *                   thereof)
	 * @param parameters a {@link ModelParameters}-instance to retrieve
	 *                   parameter-values from (e.g. learnable transition probabilities)
	 */
	public HMMStateWithStrandAndParameters(String name, boolean strand, ModelParameters parameters) {
		super(name);
		this.strand = strand;
		this.parameters = parameters;
	}

	/**
	 * @return whether this state is concerned with the forward strand
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
