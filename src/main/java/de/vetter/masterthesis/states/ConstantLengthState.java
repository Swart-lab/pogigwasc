package de.vetter.masterthesis.states;

import java.util.Iterator;

public abstract class ConstantLengthState extends HMMState {
	protected int length;
	
	public ConstantLengthState(String name, int length) {
		super(name);
		this.length = length;
	}
	
	public int getLength() {
		return length;
	}
	
	/**
	 * constant length -> only one previous l' to check
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
