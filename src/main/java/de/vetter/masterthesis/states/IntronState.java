package de.vetter.masterthesis.states;

import java.util.Iterator;

import de.vetter.masterthesis.Utilities;

/**
 * Implements a state for introns, which only considers the emmission-length
 * when computing the emmission-probability (at this point)
 * 
 * @author David Emanuel Vetter
 */
public class IntronState extends HMMState {

	private boolean strand;
	private static final int LAMBDA = 19;
	private double[] BASE_FREQUENCIES = new double[] {0.46, 0.10, 0.355, 0.085};
	
	private final static int min = 12; 
	private final static int mid = 19;
	private final static int max = 30;
	private final static double r = 0.5;
	private final static double s = 0.3;
	private final static double p = 0.45;
	private double logA;
	
	// TODO: this might very well be overfitting!
	private double[] LENGTH_PROBABILITIES = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 0.2188, 0.2812, 0.1042, 0.2292, 0.1146, 0.0208, 0.0104, 0.0104, 0.0, 0.0, 0.0104, 0.0,
			0.0, 0.0 };
	
	private void computeLogBaseTermForFiniteBidirectionalGeometricDistribution() {
		double leftSum = 1 - Math.pow(s, mid - min + 1);
		leftSum = leftSum / (1-s);
		double rightSum = 1 - Math.pow(r, max - mid + 1);
		rightSum = rightSum / (1-r);
		
		logA = leftSum + p*rightSum;
		
		logA = -Math.log(logA);
	}
	
	private double logProbability(int length) {
		double result = Double.NEGATIVE_INFINITY;
		if(length == mid)
			result = logA + Math.log(1 + p);
		
		if(length >= min && length < mid)
			result = logA + (length - min) * Math.log(s);
		
		if(length > mid && length <= max)
			result = logA + Math.log(p) + (length - mid) * Math.log(r);
		
		return result;
	}
	
	public IntronState(String name, boolean strand) {
		super(name);
		this.strand = strand;
		computeLogBaseTermForFiniteBidirectionalGeometricDistribution();
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() < 4 || newEmission.length() >= LENGTH_PROBABILITIES.length)
			return Double.NEGATIVE_INFINITY;
		
		if(!strand)
			newEmission = Utilities.reverseComplement(newEmission);
		
		if(! (newEmission.startsWith("GT") && newEmission.endsWith("AG")))
			return Double.NEGATIVE_INFINITY;
		
		// lambda = 19
		// log Poisson = k log(lambda) - lambda - log(k!)
		/*
		double lengthProb = newEmission.length() * Math.log(LAMBDA) - LAMBDA;
		poisson -= Utilities.logFactorial(newEmission.length()); // TODO: might abbreviate here
		*/
		// double lengthProb = Math.log(LENGTH_PROBABILITIES[newEmission.length()]);
		double lengthProb = logProbability(newEmission.length());
		
		double baseUsage = 0;
		// GT AG have probability 1
		for(char b : newEmission.substring(2, newEmission.length() - 2).toCharArray())
			baseUsage += BASE_FREQUENCIES[Utilities.baseToIndex(b)];
		
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
