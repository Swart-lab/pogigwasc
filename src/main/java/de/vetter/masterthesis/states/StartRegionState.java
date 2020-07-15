package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.Utilities;

public class StartRegionState extends HMMState {

	private boolean strand;
	public static final int START_REGION_LENGTH = 6;
	
	/** TCAG */
	private static final double[] UPSTREAM_BASE_USAGE = new double[] {0.15, 0.02, 0.78, 0.05};
	
	public StartRegionState(String name, boolean strand) {
		super(name);
		this.strand = strand;
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		// TODO Auto-generated method stub
		if(newEmission.length() != START_REGION_LENGTH) {
			return Double.NEGATIVE_INFINITY;
		}
		
		if(!strand) {
			newEmission = Utilities.reverseComplement(newEmission);
		}
		
		if(!newEmission.endsWith("ATG")) {
			return Double.NEGATIVE_INFINITY;
		}
		
		double result = 0;
		for(int i = 0; i < START_REGION_LENGTH - 3; i++) {
			result += Math.log(UPSTREAM_BASE_USAGE[Utilities.baseToIndex(newEmission.charAt(i))]);
		}
		return result;
	}

	/**
	 * StartRegion allows only size 6 nt (3 before ATG, then ATG)
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - START_REGION_LENGTH);
					
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
