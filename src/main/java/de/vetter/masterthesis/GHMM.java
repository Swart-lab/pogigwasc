package de.vetter.masterthesis;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements a generalized Hidden Markov Model, which always has an initial
 * state (index 0) and a terminal state (index 1)
 * 
 * @author David Emanuel Vetter
 */
public class GHMM {

	private double[][] logTransitions;
	private List<HMMState> states;
	
	/**
	 * Sets up a GHMM with only a start and an end-state
	 */
	public GHMM() {
		states = new ArrayList<HMMState>();
		// TODO: implement
	}
	
	/**
	 * Builds an initial, valid Transition-Matrix (which will be the deterministic
	 * linear march from start through all the states to the terminal state)
	 * 
	 * Call this 
	 */
	public void initialiseTransitionMatrix() {
		// TODO: implement
	}
	
	/**
	 * @return whether the transition matrix has been initialised. 
	 */
	public boolean isSetTransitionMatrix() {
		// TODO: implement
		return true;
	}
	
	/**
	 * Checks whether the transition-matrix contains valid entries
	 * 
	 * @return Whether the log-transition-probability matrix is actually that, i.e.
	 *         whether the entries are normalised
	 */
	public boolean checkTransitions() {
		// TODO: implement
		return true;
	}
	
	/**
	 * Sets the probability (!) of transitioning from->to.
	 * Call {@link initialiseTransitionMatrix} before calling this. 
	 * 
	 * @param from
	 * @param to
	 * @param probability a probability (NOT logarithm thereof)
	 */
	public void setTransitionProbability(int from, int to, double probability) {
		// TODO: implement
	}
	
	/**
	 * Normalises the given column of the probability matrix
	 * @param state
	 */
	public void normaliseExitProbabilities(int state) {
		// TODO: implement
	}
	
	/**
	 * Does NOT check whether the matrix is valid. Assert that before calling this method!
	 * 
	 * @param from: index of a state
	 * @param to: index of a state
	 * @return: log_2 of the probability of transitioning from->to.
	 */
	public double getLogTransitionProbability(int from, int to) {
		return logTransitions[from][to];
	}
	
	/**
	 * Computes the log of the emission-probability e (as specified by Stanke) of
	 * emitting sequence newEmission when in state current, having before been in
	 * state previous and having hitherto (from the start) emitted priorEmission
	 * 
	 * @param previous
	 * @param current
	 * @param priorEmission
	 * @param newEmission
	 * @return
	 */
	public double getLogEmissionProbability(int previous, int current, String priorEmission, String newEmission) {
		// TODO: implement, ask current state
		return 0;
	}
	
	
	
	/**
	 * Adds a new State
	 * @param newState
	 */
	public void addState(HMMState newState) {
		// TODO: implement
	}
	
	/**
	 * @param index
	 * @return the State specified by that index (0 is the initial state, 1 is the terminal state)
	 */
	public HMMState getState(int index) {
		return states.get(index);
	}
	
	/**
	 * Removes the state specified by index, but CANNOT remove the start nor the
	 * end-state
	 * 
	 * @param index of the state to be removed
	 */
	public void removeState(int index) {
		// TODO: implement
	}
	
}
