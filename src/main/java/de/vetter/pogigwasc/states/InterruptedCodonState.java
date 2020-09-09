package de.vetter.pogigwasc.states;

import java.util.Iterator;

import de.vetter.pogigwasc.ModelParameters;
import de.vetter.pogigwasc.Utilities;

/**
 * State for modeling codons interrupted by an intron: Such a codo-part has a
 * length (1 or 2 nt) and is either to the left (pre-intron) or to the right
 * (post-intron) of the intron. Since the emission is a codon, it uses the same
 * parameters as CDS does.
 * 
 * @author David Emanuel Vetter
 *
 */
public class InterruptedCodonState extends HMMStateWithStrandAndParameters {

	private boolean left;
	private int length;
	
	public InterruptedCodonState(String name, boolean strand, boolean left, int length, ModelParameters parameters) {
		super(name, strand, parameters);
		this.left = left;
		this.length = length;
	}
	
	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if (newEmission.length() != length)
			return Double.NEGATIVE_INFINITY;

		double logProbability = 0;
		
		if(isReverse()) {
			newEmission = Utilities.reverseComplement(newEmission);
			for (int i = 0; i < length; i++) {
				logProbability += parameters.getLogBaseProbabilityCDS(newEmission.charAt(i), left ? i + 3 - getLength() : i);
			}
		} else {
			for (int i = 0; i < length; i++) {
				logProbability += parameters.getLogBaseProbabilityCDS(newEmission.charAt(i), left ? i : (i + 3 - getLength()));
			}
		}
		return logProbability;
	}

	/**
	 * @return the one permissible length
	 */
	public int getLength() {
		return length;
	}
	
	/**
	 * It is the constant length
	 * @see de.vetter.pogigwasc.states.HMMState#getSupremumPermissibleEmissionLength()
	 */
	@Override
	public int getSupremumPermissibleEmissionLength() {
		return length;
	}
	
	/**
	 * constant length -> only one previous l' to check
	 * @see de.vetter.pogigwasc.states.HMMState#iteratePermissibleLPrimes(int)
	 */
	@Override
	public Iterable<Integer> iteratePermissibleLPrimes(final int l) {
		return new Iterable<Integer>() {

			@Override
			public Iterator<Integer> iterator() {
				return new Iterator<Integer>() {
					private int currentLPrime = Math.max(0, l - length);
					
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
