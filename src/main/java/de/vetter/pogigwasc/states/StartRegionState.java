package de.vetter.pogigwasc.states;

import java.util.Iterator;

import de.vetter.pogigwasc.ModelParameters;
import de.vetter.pogigwasc.Utilities;

public class StartRegionState extends HMMStateWithStrandAndParameters {
	
	public StartRegionState(String name, boolean strand, ModelParameters parameters) {
		super(name, strand, parameters);
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() != parameters.getStartRegionSize()) {
			return Double.NEGATIVE_INFINITY;
		}
		
		if(isReverse()) {
			newEmission = Utilities.reverseComplement(newEmission);
		}
		
		if(!newEmission.endsWith("ATG")) {
			return Double.NEGATIVE_INFINITY;
		}
		
		double result = 0;
		for(int i = 0; i < parameters.getStartRegionSize() - 3; i++) {
			result += parameters.getLogBaseProbabilityStartRegion(newEmission.charAt(i));
		}
		
		return result;
	}

	/**
	 * @see de.vetter.pogigwasc.states.HMMState#getSupremumPermissibleEmissionLength()
	 */
	@Override
	public int getSupremumPermissibleEmissionLength() {
		return parameters.getStartRegionSize();
	}
	
	/**
	 * StartRegion allows only size 6 nt (3 before ATG, then ATG)
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLPrimes(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - parameters.getStartRegionSize());
					
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
