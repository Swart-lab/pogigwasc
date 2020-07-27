package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.ModelParameters;
import de.vetter.masterthesis.Utilities;

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
	 * @see de.vetter.masterthesis.states.HMMState#computeLogEmissionProbability(int, java.lang.String, java.lang.String)
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
	 * CDS allows only steps of size 3 -> thus, iterate only once
	 * @see de.vetter.masterthesis.states.HMMState#iteratePermissibleLengths(int)
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
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
