package de.vetter.masterthesis.states;

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
	
	public IntronState(String name, boolean strand) {
		super(name);
		this.strand = strand;
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		if(newEmission.length() < 4)
			return Double.NEGATIVE_INFINITY;
		
		if(!strand)
			newEmission = Utilities.reverseComplement(newEmission);
		
		if(! (newEmission.startsWith("GT") && newEmission.endsWith("AG")))
			return Double.NEGATIVE_INFINITY;
		
		// lambda = 19
		// log Poisson = k log(lambda) - lambda - log(k!) 
		double poisson = newEmission.length() * Math.log(LAMBDA) - LAMBDA;
		poisson -= Utilities.logFactorial(newEmission.length()); // TODO: might abbreviate here
		
		double baseUsage = 0;
		for(char b : newEmission.toCharArray())
			baseUsage += BASE_FREQUENCIES[Utilities.baseToIndex(b)];
		
		return poisson + baseUsage;
	}

}
