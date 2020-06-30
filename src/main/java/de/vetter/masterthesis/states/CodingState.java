package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.Utilities;

public class CodingState extends HMMState {
	private boolean strand;

	/**
	 * 
	 * @param name
	 * @param strand whether on forward strand
	 */
	public CodingState(String name, boolean strand) {
		super(name);
		this.strand = strand;
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		// TODO Codon usage or some approximation (maybe always use empirical distr. of base 1, to avoid overfitting?)
		if(newEmission.length() != 3)
			return Double.NEGATIVE_INFINITY;
		
		if(!strand)
			newEmission = Utilities.reverseComplement(newEmission);
		
		double result = 0;
		for(int i = 0; i < 3; i++)
			result += Utilities.getLogBaseProbabilityCDS(newEmission.charAt(i), i);
		
		return result;
	}
	
	/**
	 * CDS allows only steps of size 3 -> thus, iterate only once
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
