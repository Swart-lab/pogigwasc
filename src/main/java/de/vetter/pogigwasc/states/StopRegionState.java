package de.vetter.pogigwasc.states;

import java.util.Iterator;

import de.vetter.pogigwasc.ModelParameters;
import de.vetter.pogigwasc.Utilities;

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
		double[] empirical = { 0.03837953091684435, 0.029850746268656716, 0.0255863539445629, 0.019189765458422176,
				0.014925373134328358, 0.008528784648187633, 0.010660980810234541, 0.0021321961620469083,
				0.029850746268656716, 0.0255863539445629, 0.0042643923240938165, 0.017057569296375266,
				0.010660980810234541, 0.008528784648187633, 0.0, 0.01279317697228145, 0.023454157782515993,
				0.014925373134328358, 0.021321961620469083, 0.0042643923240938165, 0.029850746268656716,
				0.01279317697228145, 0.023454157782515993, 0.006396588486140725, 0.008528784648187633,
				0.006396588486140725, 0.019189765458422176, 0.0042643923240938165, 0.008528784648187633,
				0.0021321961620469083, 0.0021321961620469083, 0.0021321961620469083, 0.0255863539445629,
				0.014925373134328358, 0.006396588486140725, 0.021321961620469083, 0.017057569296375266,
				0.006396588486140725, 0.008528784648187633, 0.0042643923240938165, 0.0255863539445629,
				0.017057569296375266, 0.05543710021321962, 0.042643923240938165, 0.008528784648187633,
				0.006396588486140725, 0.042643923240938165, 0.008528784648187633, 0.02771855010660981,
				0.01279317697228145, 0.006396588486140725, 0.006396588486140725, 0.031982942430703626,
				0.0021321961620469083, 0.021321961620469083, 0.006396588486140725, 0.019189765458422176,
				0.019189765458422176, 0.02771855010660981, 0.017057569296375266, 0.008528784648187633,
				0.006396588486140725, 0.021321961620469083, 0.006396588486140725 };
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
	 * @see de.vetter.pogigwasc.states.HMMState#getSupremumPermissibleEmissionLength()
	 */
	@Override
	public int getSupremumPermissibleEmissionLength() {
		return parameters.getStopRegionSize();
	}
	
	/**
	 * StopRegion allows only steps of size 24 -> thus, iterate only once
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLPrimes(final int l) {
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
