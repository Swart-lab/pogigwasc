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
		
		/* old mode:
		double[] empirical = { 0.027093596059113302, 0.027093596059113302, 0.03201970443349754, 0.022167487684729065,
				0.019704433497536946, 0.007389162561576354, 0.014778325123152709, 0.0, 0.029556650246305417,
				0.014778325123152709, 0.0, 0.017241379310344827, 0.0024630541871921183, 0.012315270935960592, 0.0,
				0.017241379310344827, 0.034482758620689655, 0.012315270935960592, 0.024630541871921183,
				0.0049261083743842365, 0.014778325123152709, 0.009852216748768473, 0.017241379310344827,
				0.0024630541871921183, 0.0024630541871921183, 0.0049261083743842365, 0.027093596059113302,
				0.007389162561576354, 0.0024630541871921183, 0.0024630541871921183, 0.0, 0.0, 0.03201970443349754,
				0.017241379310344827, 0.019704433497536946, 0.027093596059113302, 0.022167487684729065,
				0.0024630541871921183, 0.014778325123152709, 0.0, 0.027093596059113302, 0.022167487684729065,
				0.05665024630541872, 0.03201970443349754, 0.007389162561576354, 0.0024630541871921183,
				0.03940886699507389, 0.0049261083743842365, 0.029556650246305417, 0.014778325123152709,
				0.007389162561576354, 0.0049261083743842365, 0.03201970443349754, 0.0049261083743842365,
				0.022167487684729065, 0.0024630541871921183, 0.024630541871921183, 0.019704433497536946,
				0.046798029556650245, 0.027093596059113302, 0.0049261083743842365, 0.0049261083743842365,
				0.012315270935960592, 0.0024630541871921183 };
		result = 0;
		for(int i = 0; i < parameters.getStopRegionSize() - 3; i += 3) {
			int codonIndex = Utilities.baseToIndex(newEmission.charAt(i)) * 16;
			codonIndex += Utilities.baseToIndex(newEmission.charAt(i+1)) * 4;
			codonIndex += Utilities.baseToIndex(newEmission.charAt(i+2));
			result += empirical[codonIndex];
		}
		// END old mode */
		
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
