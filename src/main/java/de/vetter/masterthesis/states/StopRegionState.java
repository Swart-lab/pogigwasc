package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.ModelParameters;
import de.vetter.masterthesis.Utilities;

public class StopRegionState extends HMMStateWithStrandAndParameters {
	
	public StopRegionState(String name, boolean strand, ModelParameters parameters) { 
		super(name, strand, parameters);
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() != parameters.getStopRegionSize())
			return Double.NEGATIVE_INFINITY;
		
		if(isReverse()) {
			newEmission = Utilities.reverseComplement(newEmission);
		}
		
		if(!newEmission.endsWith("TGA"))
			return Double.NEGATIVE_INFINITY;
		
		double result = 0;
		for(int i = 0; i < parameters.getStopRegionSize() - 3; i += 3)
			result += parameters.getLogCodonProbabilityStopRegion(newEmission.substring(i, i+3));
		
		return result;
	}
	
	/**
	 * StopRegion allows only steps of size 24 -> thus, iterate only once
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - parameters.getStopRegionSize());
					
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
