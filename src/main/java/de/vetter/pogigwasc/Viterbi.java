package de.vetter.pogigwasc;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements the Viterbi algorithm as described by Stanke (AUGUSTUS): Construct
 * a Viterbi-instance for a GHMM and a given sequence, to compute the
 * Viterbi-variables of that model on that sequence, and yield the MAP-parses
 * 
 * @author David Emanuel Vetter
 *
 */
public class Viterbi {
	private double[][] viterbiVariables;
	private GHMM model;
	private String sequence;
	
	private boolean abbreviating = false;
	
	/**
	 * Constructor: Checks given model for transition-validity
	 * @param model: Has to have its transition matrix set and valid
	 * @throws IllegalArgumentException if the model is (trivially) invalid.
	 */
	public Viterbi(GHMM model, String sequence) {
		// check model validity and potentially throw exceptions
		if(!(model.isSetTransitionMatrix() && model.checkTransitions())) {
			throw new IllegalArgumentException("The given GHMM does not have a valid transition matrix (not set, or not normalised).");
		}
		this.model = model;
		this.sequence = sequence;
	}
	
	/**
	 * This determines whether the emission history is constructed and handed to the
	 * states when computing the emission probabilities -- this is often not needed,
	 * in particular, i do not make use of the emission history here (so far)
	 * 
	 * @param abbreviate whether to omit the emission history (make sure the model
	 *                   does not need it!)
	 */
	public void setAbbreviating(boolean abbreviate) {
		this.abbreviating = abbreviate;
	}
	
	/**
	 * @return whether this viterbi-instance is abbreviating the
	 *         emission-probability-computation by omitting the emission-history
	 */
	public boolean isAbbreviating() {
		return abbreviating;
	}
	
	private void computeViterbiVariables() {
		int stateCount = model.getNumberOfStates();
		// NOTE! Compute in logarithm, i.e. probability 1 is entry 0 etc; Addition instead of multiplication
		viterbiVariables = new double[stateCount][sequence.length() + 1];
		
		// Initialisation:
		viterbiVariables[0][0] = 0; 
		for(int state = 1; state < stateCount; state++) {
			viterbiVariables[state][0] = Double.NEGATIVE_INFINITY;
		}
		
		// collect the longest possible emission that could start the sequence:
		int maxInitialEmissionLength = 0;
		for(int q = 2; q < stateCount; q++) {
			if(model.getLogTransitionProbability(0, q) > Double.NEGATIVE_INFINITY)
				maxInitialEmissionLength = Math.max(maxInitialEmissionLength, model.getState(q).getSupremumPermissibleEmissionLength());
		}
		
		System.out.println(" Viterbi: Maximum initial emission length = " + maxInitialEmissionLength + "\n");
		
		// 'Recursion'
		for(int l = 1; l < sequence.length() + 1; l++) {
			for(int q = 0; q < stateCount; q++) {
				if(q == 1)
					continue;
				double max = Double.NEGATIVE_INFINITY;
				for (int qPrime = 0; qPrime < stateCount; qPrime++) {
					if(model.getLogTransitionProbability(qPrime, q) > Double.NEGATIVE_INFINITY && qPrime != 1) {
						if(qPrime == 0 && l < maxInitialEmissionLength) {
							max = Math.max(max,
									viterbiVariables[qPrime][0] + model.getLogTransitionProbability(qPrime, q)
											+ model.getLogEmissionProbability(0, q, "",
													sequence.substring(0, l)));
						} else {
							// q' \in Q, i.e. neither initial nor terminal state (cannot come from the terminal state)
							for(int lPrime : model.getState(q).iteratePermissibleLPrimes(l)) {
								max = Math.max(max,
										viterbiVariables[qPrime][lPrime] + model.getLogTransitionProbability(qPrime, q)
												+ model.getLogEmissionProbability(qPrime, q,
														abbreviating ? null : sequence.substring(0, lPrime),
														sequence.substring(lPrime, l)));
							}
						}
					}
				}
				
				viterbiVariables[q][l] = max;
			}
			if (l % 5000 == 0) {
				System.out.println(" progress: l=" + l + "/" + (sequence.length() + 1));
			}
				
		}
		
		System.out.println("Computed variables!");
	}
	
	/**
	 * Computes all most likely parses (and reports progress in System.out)
	 * @return A List of Parses (where each parse is a list of pairs of a state and the emission-length from that state)
	 */
	public List<Parse> computeParses() {
		computeViterbiVariables();
		ArrayList<ViterbiSeed> workLoad = new ArrayList<ViterbiSeed>();
		ArrayList<ViterbiSeed> finished = new ArrayList<ViterbiSeed>(); // here put all seeds that reach l=0
		
		// Initialisation
		double initialMax = Double.NEGATIVE_INFINITY;
		for(int q1 = 2; q1 < model.getNumberOfStates(); q1++) {
			double candidate = viterbiVariables[q1][sequence.length()] + model.getLogTransitionProbability(q1, 1);
			if(initialMax <= candidate) {
				if(initialMax < candidate) {
					// TODO: is hard exact maximum sensible, or should there be a little double-tolerance?
					initialMax = candidate;
					workLoad.clear();
				}
				workLoad.add(new ViterbiSeed(model, sequence, q1, sequence.length(), null, abbreviating)); // the ends of the parses point to null.
			}
		}
		
		// 'Recursion': working through the viterbi-variable-matrix
		while(!workLoad.isEmpty()) {
			ViterbiSeed current = workLoad.remove(0);
			List<ViterbiSeed> stepped = current.step(viterbiVariables);
			if(stepped.size() > 1)
				System.out.println("Encountered ambiguous parse: |Workload|=" + workLoad.size());
			for(ViterbiSeed s : stepped) {
				if(s.isFinished()) {
					finished.add(s);
				} else {
					workLoad.add(s);
				}
			}
		}
		
		// This corresponds to the finalisation described by Stanke on page 15.
		ArrayList<Parse> parses = new ArrayList<Parse>();
		for (ViterbiSeed startSeed : finished) {
			Parse parse = new Parse();
			ViterbiSeed current = startSeed;
			while (current.getPrevious() != null) {
				parse.add(model.getState(current.getPrevious().getQ()), current.getPrevious().getL() - current.getL());
				current = current.getPrevious();
			}

			parses.add(parse);
		}
		
		return parses;
	}
	
}
