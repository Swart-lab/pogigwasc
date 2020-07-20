package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.ModelParameters;
import de.vetter.masterthesis.Utilities;

/**
 * Implements a state for introns, which only considers the emmission-length
 * when computing the emmission-probability (at this point)
 * 
 * @author David Emanuel Vetter
 */
public class IntronState extends HMMStateWithStrandAndParameters {
		
	public IntronState(String name, boolean strand, ModelParameters parameters) {
		super(name, strand, parameters);
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() < 4 || newEmission.length() > parameters.getMaxIntronSize())
			return Double.NEGATIVE_INFINITY;
		
		if(isReverse())
			newEmission = Utilities.reverseComplement(newEmission);
		
		if(! (newEmission.startsWith("GT") && newEmission.endsWith("AG")))
			return Double.NEGATIVE_INFINITY;
		
		int length = newEmission.length();
		
		double lengthProb = parameters.getLogProbabilityIntronLength(length);
		
		double baseUsage = 0;
		/* GT AG have probability 1 */
		for(int i = 0; i < length; i++) {
			if(i < parameters.getSDSSize()) {
				baseUsage += parameters.getLogBaseProbabilitySDS(newEmission.charAt(i), i); 
			} else if (i >= length - parameters.getSASSize()) {
				baseUsage += parameters.getLogBaseProbabilitySAS(newEmission.charAt(i), i - (length - parameters.getSASSize()));
			} else {
				baseUsage += parameters.getLogBaseProbabilityIntron(newEmission.charAt(i));
			}
		}
		
		return lengthProb + baseUsage;
	}

	
	/**
	 * introns could maybe take arbitrary lengths >= 4, for performance, limit it to 50
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - 33);
					
					@Override
					public boolean hasNext() {
						return currentLPrime < l - 3;
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
