package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.ModelParameters;
import de.vetter.masterthesis.Utilities;

/**
 * Implements a state for introns, which only allows a certain range of intron
 * sizes (cf {@link ModelParameters#getMinIntronSize()} and
 * {@link ModelParameters#getMaxIntronSize()})<br>
 * 
 * Further, sequence information at the ends of the intron (the 5' splice donor
 * site SDS and the 3' splice acceptor site SAS) factor into the emission probability
 * 
 * @author David Emanuel Vetter
 */
public class IntronState extends HMMStateWithStrandAndParameters {
		
	public IntronState(String name, boolean strand, ModelParameters parameters) {
		super(name, strand, parameters);
	}

	/**
	 * Intron emission probability: Consider both length and exact sequence of the
	 * new emission: Require edges 5': GT and 3': AG, and take sequence preferences
	 * at the SDS and SAS into account. Internal base frequency (inside intron,
	 * outside SDS and SAS) will be similar to NCS base frequency. While the
	 * base-distribution depends on the position in the emission, each position is
	 * modeled as an independent random variable.
	 * 
	 * @see de.vetter.masterthesis.states.HMMState#computeLogEmissionProbability(int,
	 *      java.lang.String, java.lang.String)
	 */
	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if (newEmission.length() < parameters.getMinIntronSize()
				|| newEmission.length() > parameters.getMaxIntronSize())
			return Double.NEGATIVE_INFINITY;
		
		if(isReverse())
			newEmission = Utilities.reverseComplement(newEmission);
		
		// TODO: redundant
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
		/* Compute without respect to intron-structure! This is the old style
		baseUsage = 0;
		for(int i = 2; i < length - 2; i++) {
			baseUsage += parameters.getLogBaseProbabilityNCS(newEmission.charAt(i));
		}
		*/
		
		return lengthProb + baseUsage;
	}

	
	/**
	 * introns could maybe take arbitrary lengths >= 4, for performance, limit it to
	 * the allowed range of {l-max, ..., l-min} (note the inclusive upper end)
	 * 
	 * @see de.vetter.masterthesis.states.HMMState#iteratePermissibleLengths(int)
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLengths(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - parameters.getMaxIntronSize());
					
					@Override
					public boolean hasNext() {
						return currentLPrime < l - parameters.getMinIntronSize() + 1;
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
