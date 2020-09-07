package de.vetter.pogigwasc.states;

import java.util.Iterator;

import de.vetter.pogigwasc.ModelParameters;

public class NoncodingState extends HMMStateWithStrandAndParameters {

	public NoncodingState(String name, ModelParameters params) { super(name, true, params); }

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() > 1) {
			return Double.NEGATIVE_INFINITY;
		}
		
		return parameters.getLogBaseProbabilityNCS(newEmission.charAt(0));
	}
	
	/**
	 * NCS allows only a length of 1
	 * @see de.vetter.pogigwasc.states.HMMState#getSupremumPermissibleEmissionLength()
	 */
	@Override
	public int getSupremumPermissibleEmissionLength() {
		return 1;
	}
	
	/**
	 * NCS allows only steps of size 1 -> thus, iterate only once
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLPrimes(final int l) {
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
