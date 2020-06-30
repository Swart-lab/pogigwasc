package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.Utilities;

public class StopRegionState extends HMMState {
	private boolean strand;
	
	public StopRegionState(String name, boolean strand) { 
		super(name);
		this.strand = strand;
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() != 24)
			return Double.NEGATIVE_INFINITY;
		
		if(!strand) {
			newEmission = Utilities.reverseComplement(newEmission);
		}
		
		if(!newEmission.endsWith("TGA"))
			return Double.NEGATIVE_INFINITY;
		
		double result = 0;
		for(int i = 0; i < 21; i += 3)
			result += Utilities.getLogCodonProbabilityStopRegion(newEmission.substring(i, i+3));
		
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
					private int currentLPrime = Math.max(0, l - 24);
					
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
