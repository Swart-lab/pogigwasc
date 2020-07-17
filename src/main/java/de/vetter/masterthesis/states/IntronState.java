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
	
	private double[] BASE_FREQUENCIES = new double[] {0.46, 0.10, 0.355, 0.085};
	
	private int SDS_SIZE = 5;
	/** 5' end of intron (splice donor site): start with GT~A~A; (pos, base) */
	private double[][] BASE_FREQUENCIES_SDS = new double[][] {
		{0d, 0d, 0d, 1d},
		{1d, 0d, 0d, 0d},
		{0.02, 0.02, 0.94, 0.02},
		{0.02, 0.02, 0.94, 0.02},
		{0.35, 0.05, 0.15, 0.45}
	};
	
	private int SAS_SIZE = 4;
	/** 3' end of intron: end with ~T/A~T; (pos, base) */
	private double[][] BASE_FREQUENCIES_SAS = new double[][] {
		{0.45, 0.05, 0.35, 0.15},
		{0.65, 0.25, 0.05, 0.05},
		{0d, 0d, 1d, 0d},
		{0d, 0d, 0d, 1d}
	};
	
	private final static int MIN = 12; 
	private final static int MID = 18;
	private final static int MAX = 30;
	private final static double R = 0.6;
	private final static double S = 0.4;
	private final static double P = 0.55;
	private double logA;
	

	
	private void computeLogBaseTermForFiniteBidirectionalGeometricDistribution() {
		double leftSum = 1 - Math.pow(S, MID - MIN); //  + 1);
		leftSum = leftSum / (1-S);
		double rightSum = 1 - Math.pow(R, MAX - MID + 1);
		rightSum = rightSum / (1-R);
				
		logA = leftSum + P*rightSum;
		
		logA = -Math.log(logA);
		
		/* 
		System.out.println("Intron-state: length-distribution:");
		System.out.println("l,p,poisson,empirical");
		double sum = 0;
		for(int k = 10; k <= MAX + 1; k++) {
			double poisson = k * Math.log(LAMBDA) - LAMBDA;
			poisson -= Utilities.logFactorial(k);
			System.out.println(k + "," + Math.exp(logProbability(k)) + "," + Math.exp(poisson) + ","
					+ (k < LENGTH_PROBABILITIES.length ? LENGTH_PROBABILITIES[k] : 0));
			sum += Math.exp(logProbability(k));
		}
		System.out.println("sum=" + sum);
		*/
	}
	
	private double logProbability(int length) {
		double result = Double.NEGATIVE_INFINITY;
		//if(length == MID)
		//	result = logA + Math.log(1 + P);
		// and 2nd if: MID < length ...
		
		if(MIN <= length && length < MID)
			result = logA + (MID - length - 1) * Math.log(S);
		
		if(MID <= length && length <= MAX)
			result = logA + Math.log(P) + (length - MID) * Math.log(R);
		
		return result;
	}
	
	public IntronState(String name, boolean strand) {
		super(name);
		this.strand = strand;
		computeLogBaseTermForFiniteBidirectionalGeometricDistribution();
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() < 4 || newEmission.length() > MAX)
			return Double.NEGATIVE_INFINITY;
		
		if(!strand)
			newEmission = Utilities.reverseComplement(newEmission);
		
		if(! (newEmission.startsWith("GT") && newEmission.endsWith("AG")))
			return Double.NEGATIVE_INFINITY;
		
		int length = newEmission.length();
		
		// lambda = 19
		// log Poisson = k log(lambda) - lambda - log(k!)
		/*
		double lengthProb = newEmission.length() * Math.log(LAMBDA) - LAMBDA;
		lengthProb -= Utilities.logFactorial(newEmission.length()); 
		*/
		// double lengthProb = Math.log(LENGTH_PROBABILITIES[newEmission.length()]);
		double lengthProb = logProbability(length);
		
		double baseUsage = 0;
		/* GT AG have probability 1 */
		/*
		for(char b : newEmission.substring(2, newEmission.length() - 2).toCharArray())
			baseUsage += Math.log(BASE_FREQUENCIES[Utilities.baseToIndex(b)]);
		*/
		
		for(int i = 0; i < length; i++) {
			if(i < SDS_SIZE) {
				baseUsage += Math.log(BASE_FREQUENCIES_SDS[i][Utilities.baseToIndex(newEmission.charAt(i))]);
			} else if (i >= length - SAS_SIZE) {
				baseUsage += Math.log(BASE_FREQUENCIES_SAS[i - (length - SAS_SIZE)][Utilities.baseToIndex(newEmission.charAt(i))]);
			} else {
				baseUsage += Math.log(BASE_FREQUENCIES[Utilities.baseToIndex(newEmission.charAt(i))]);
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
