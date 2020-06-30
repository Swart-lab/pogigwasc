package de.vetter.masterthesis.states;

import java.util.Iterator;

/**
 * Abstract class to generalise over different types of states (mainly
 * differentiated by their emmission distributions)
 * 
 * @author David Emanuel Vetter
 *
 */
public abstract class HMMState {
	
	private String name; // mainly for output/human-readability
	
	public HMMState(String name) {
		this.name = name;
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	/**
	 * Override this in inheriting states to limit the iteration to only sensible lPrimes.
	 * 
	 * @param l from viterbi-recursion: this is the upper limit (exclusive) of the values that are to be iterated.
	 * @return an iterator iterating from 0 to l-1 (inclusive), in the basic form
	 */
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = 0;
					
					@Override
					public boolean hasNext() {
						return currentLPrime < l;
					}

					@Override
					public Integer next() {
						return currentLPrime++;
					}
				};
			}
			
		};
	}
	
	/**
	 * Computes the log_2 of the probability of emitting the newEmission when in
	 * current State, having before been in previousState and having (from the
	 * start) previously emitted emissionHistory.
	 * 
	 * To be implemented by concrete instances (probably ad-hoc-inheritors)
	 * 
	 * @param previousState
	 * @param emissionHistory
	 * @param newEmission
	 * @return log_2 of e as specified by Stanke.
	 */
	public abstract double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission);

}
