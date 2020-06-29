package de.vetter.masterthesis;

import java.util.ArrayList;
import java.util.List;

import de.vetter.masterthesis.states.HMMState;
import de.vetter.masterthesis.states.InitialState;
import de.vetter.masterthesis.states.TerminalState;

/**
 * Implements a generalized Hidden Markov Model, which always has an initial
 * state (index 0) and a terminal state (index 1).
 * 
 * Usage: Construct, then add all desired states via
 * {@link #addState(HMMState)}, then call {@link #initialiseTransitionMatrix()}
 * and implement the transition-probabilities via calls of
 * {@link #setTransitionProbability(int, int, double)}, and potentially
 * normalise the matrix via {{@link #normaliseExitProbabilities(int)}. Finally,
 * you might want to check the validity of the matrix via
 * {{@link #checkTransitions()}
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
		addState(new InitialState());
		addState(new TerminalState());
	}
	
	/**
	 * Builds an initial, valid Transition-Matrix (which will be the deterministic
	 * linear march from start through all the states to the terminal state)
	 * 
	 * Call this after adding all the states
	 */
	public void initialiseTransitionMatrix() {
		// TODO: add option to initialise all-zero-matrix
		int n = states.size();
		logTransitions = new double[n][n];
		for(int from = 0; from < n; from++) {
			for(int to = 0; to < n; to++) {
				if(from+1 == to && from > 1) {
					logTransitions[from][to] = 0;
				} else {
					logTransitions[from][to] = Double.NEGATIVE_INFINITY;
				}
			}
		}
		logTransitions[0][Math.min(2, n-1)] = 0; // if size=2 -> have to set 0->1 to log(1)
		logTransitions[1][1] = 0; // terminal state is never left
		logTransitions[n-1][1] = 0; // from last state, potentially terminal state, always to terminal state
	}
	
	/**
	 * Sets ALL entries to -infty (i.e. probability 0).
	 * Use carefully!
	 */
	public void clearTransitionMatrix() {
		int n = states.size();
		for(int from = 0; from < n; from++) {
			for(int to = 0; to < n; to++) {
				logTransitions[from][to] = Double.NEGATIVE_INFINITY;
			}
		}
	}
	
	/**
	 * @return whether the transition matrix has been initialised. 
	 */
	public boolean isSetTransitionMatrix() {
		return logTransitions != null;
	}
	
	/**
	 * Checks whether the transition-matrix is valid. Requires that it actually be set.
	 * 
	 * @return Whether the log-transition-probability matrix is actually that, i.e.
	 *         whether the entries are normalised
	 */
	public boolean checkTransitions() {
		for(int from = 0; from < states.size(); from++) {
			double summedExits = 0;
			for(int to = 0; to < states.size(); to++) {
				summedExits += Math.exp(logTransitions[from][to]);
			}
			if(Math.abs(summedExits - 1) > 1e-9) {
				System.out.println("Invalid at: " + from + " (sum=" + summedExits + ")");
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Sets the probability (!, not log probability) of transitioning from->to. Can
	 * instead give weights here and then normalise via
	 * {@link #normaliseExitProbabilities(int)}. Call
	 * {@link initialiseTransitionMatrix} before calling this.
	 * 
	 * 
	 * @param from
	 * @param to
	 * @param probability a probability (NOT logarithm thereof)
	 */
	public void setTransitionProbability(int from, int to, double probability) {
		logTransitions[from][to] = Math.log(probability);
	}
	
	/**
	 * Normalises the given column of the probability matrix, unless all entries are -Infty. in that case turns the state into simple forward.
	 * @param state
	 */
	public void normaliseExitProbabilities(int state) {
		double normaliser = 0;
		for(int to = 0; to < states.size(); to++) {
			normaliser += Math.exp(logTransitions[state][to]);
		}
		if(normaliser > 0) {
			normaliser = Math.log(normaliser);
			for(int to = 0; to < states.size(); to++) {
				logTransitions[state][to] -= normaliser; 
			}
		} else {
			// The row-sum is -infty, then turn state into simple forward (to next state i+1, or terminal, if i+1 does not exist)
			logTransitions[state][state+1 == states.size() ? 1 : state+1] = 0;
		}
	}
	
	/**
	 * Does NOT check whether the matrix is valid, nor whether it is set. Assert that before calling this method!
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
		return states.get(current).computeLogEmissionProbability(previous, priorEmission, newEmission);
	}
	
	
	
	/**
	 * Adds a new State; The transition matrix will be unset, since it will have to be reinitialised anyways
	 * @param newState: Must not be null!
	 * @throws NullPointerException if state is null
	 */
	public void addState(HMMState newState) {
		if(newState != null) {
			logTransitions = null;
			states.add(newState);
		} else {
			throw new NullPointerException("Expected a non-null state to add to the states-list, but was null");
		}
	}
	
	/**
	 * @param index
	 * @return the State specified by that index (0 is the initial state, 1 is the terminal state)
	 */
	public HMMState getState(int index) {
		return states.get(index);
	}
	
	/**
	 * @return the total number of states; this includes the initial and terminal
	 *         state. Number of nontrivial states is result-2
	 */
	public int getNumberOfStates() {
		return states.size();
	}
	
	/**
	 * Removes the state specified by index, but CANNOT remove the start nor the
	 * end-state. Like {@link #addState(HMMState)}, unsets the transition matrix
	 * 
	 * @param index of the state to be removed (must be 2 or greater)
	 * @return returns the removed state (e.g. for shifting a state to the end)
	 * @throws IndexOutOfBoundsException if you try to remove the start or end-state
	 */
	public HMMState removeState(int index) {
		if(index < 2) {
			throw new IndexOutOfBoundsException("Tried to remove state " + index
					+ ", but only indices > 1 are allowed (start and end-state cannot be removed)");
		}
		HMMState removed = states.remove(index);
		logTransitions = null;
		return removed;
	}
	
	
	public String toString() {
		StringBuffer result = new StringBuffer();
		result.append("GHMM\n");
		result.append("states: " + states + "\n");
		result.append("Transition-probabilities are" + (isSetTransitionMatrix() ? " " : "not ") + "set");
		if(isSetTransitionMatrix()) {
			result.append(", and " + (checkTransitions() ? "valid" : "invalid"));
			
			for(int from = 0; from < states.size(); from++) {
				result.append("\n");
				for(int to = 0; to < states.size(); to++) {
					result.append("  ");
					result.append(Math.exp(logTransitions[from][to]));
				}
			}
		}
		
		return result.toString();
	}
}
