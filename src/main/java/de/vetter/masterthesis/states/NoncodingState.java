package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.Utilities;

public class NoncodingState extends HMMState {

	// TODO: preferably read this from parameter-file, also make symmetric?
	private final double[] BASE_FREQUENCIES = new double[] {0.46, 0.10, 0.355, 0.085};
	
	public NoncodingState(String name) { super(name); }

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() > 1) {
			return Double.NEGATIVE_INFINITY;
		}
		return Math.log(BASE_FREQUENCIES[Utilities.baseToIndex(newEmission.charAt(0))]);
	}
	
	/**
	 * NCS allows only steps of size 1 -> thus, iterate only once
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - 1);
					
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

}
