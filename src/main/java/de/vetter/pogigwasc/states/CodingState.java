package de.vetter.pogigwasc.states;

import java.util.Iterator;

import de.vetter.pogigwasc.ModelParameters;
import de.vetter.pogigwasc.Utilities;

/**
 * GHMM-State representation of a CDS: Always emits a single codon (permissible lengths={3}) 
 * @author David Emanuel Vetter
 */
public class CodingState extends HMMStateWithStrandAndParameters {

	/**
	 * @see HMMStateWithStrandAndParameters#HMMStateWithStrandAndParameters(String, boolean, ModelParameters)
	 */
	public CodingState(String name, boolean strand, ModelParameters parameters) {
		super(name, strand, parameters);
	}

	/**
	 * @see de.vetter.pogigwasc.states.HMMState#computeLogEmissionProbability(int, java.lang.String, java.lang.String)
	 */
	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() != 3)
			return Double.NEGATIVE_INFINITY;
		
		if(isReverse())
			newEmission = Utilities.reverseComplement(newEmission);
		
		double result = 0;
		for(int i = 0; i < 3; i++)
			result += parameters.getLogBaseProbabilityCDS(newEmission.charAt(i), i);
		
		return result;
	}
	
	/**
	 * CDS allows only a length of 3
	 * @see de.vetter.pogigwasc.states.HMMState#getSupremumPermissibleEmissionLength()
	 */
	@Override
	public int getSupremumPermissibleEmissionLength() {
		return 3;
	}
	
	/**
	 * CDS allows only steps of size 3 -> thus, iterate only once
	 * @see de.vetter.pogigwasc.states.HMMState#iteratePermissibleLPrimes(int)
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLPrimes(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - 3);
					
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
