package de.vetter.pogigwasc.states;

import java.util.Iterator;

/**
 * States like the old start-codon only emit exactly one sequence with probability one.
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
		return newEmission.equals(sequence) ? 0 : Double.NEGATIVE_INFINITY;
	}

	
	public String getSequence() {
		return sequence;
	}
	
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	/**
	 * Allows only the length of this states fixed sequence
	 * @see de.vetter.pogigwasc.states.HMMState#getSupremumPermissibleEmissionLength()
	 */
	@Override
	public int getSupremumPermissibleEmissionLength() {
		return sequence.length();
	}
	
	/**
	 * Can only make step of size |fixedSequence|
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLPrimes(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - sequence.length());
					
					@Override
					public boolean hasNext() {
						return currentLPrime < l;
					}

					@Override
					public Integer next() {
						int result = currentLPrime;
						currentLPrime = l;
						return result;
					}
				};
			}
			
		};
	}
}
