package de.vetter.masterthesis;

import java.util.List;

/**
 * Implements the viterbi algorithm as described by Stanke (AUGUSTUS)
 * 
 * @author David Emanuel Vetter
 *
 */
public class Viterbi {
	private double[][] viterbiVariables;
	private GHMM model;
	private String sequence;
	
	/**
	 * Constructor: Checks given model for transition-validity
	 * @param model
	 */
	public Viterbi(GHMM model, String sequence) {
		// TODO: check model validity and potentially throw exceptions
		this.model = model;
		this.sequence = sequence;
	}
	
	private void computeViterbiVariables() {
		// TODO: implement
	}
	
	public List<Pair<HMMState, Integer>> computeParse() {
		computeViterbiVariables();
		return null;
	}
	
}
