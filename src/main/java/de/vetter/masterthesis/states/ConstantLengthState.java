package de.vetter.masterthesis.states;

import java.util.Iterator;

/**
 * HMM state with one fixed emission length, and thus with a well-known
 * {@link #iteratePermissibleLengths(int)}. The exact emission probability is
 * still to be implemented by the heirs
 * 
 * @author David Emanuel Vetter
 */
public abstract class ConstantLengthState extends HMMState {
	private int length;
	
	/**
	 * @param name this state's name
	 * @param length the one permissible length
	 */
	public ConstantLengthState(String name, int length) {
		super(name);
		this.length = length;
	}
	
	/**
	 * @return the one permissible length
	 */
	public int getLength() {
		return length;
	}
	
	/**
	 * constant length -> only one previous l' to check
	 * @see de.vetter.masterthesis.states.HMMState#iteratePermissibleLengths(int)
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - length);
					
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
