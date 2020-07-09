package de.vetter.masterthesis;

import java.util.ArrayList;
import java.util.List;

import de.vetter.masterthesis.states.HMMState;

/**
 * Implements the viterbi algorithm as described by Stanke (AUGUSTUS)
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
	
	public void setAbbreviating(boolean abbreviate) {
		this.abbreviating = abbreviate;
	}
	
	/**
	 * 
	 * @return whether this viterbi-instance is abbreviating the emission-probability-computation by omitting the emission-history
	 */
	public boolean isAbbreviating() {
		return abbreviating;
	}
	
	private void computeViterbiVariables() {
		int stateCount = model.getNumberOfStates();
		// NOTE! Compute in logarithm, i.e. probability 1 is entry 0 etc; Addition instead of multiplication
		viterbiVariables = new double[stateCount][sequence.length() + 1];
		
		/*long tBefore, tAfter, tBeforeFindingMax, tAfterFindingMax;
		tBefore = System.currentTimeMillis();*/
		// Initialisation:
		viterbiVariables[0][0] = 0; 
		for(int state = 1; state < stateCount; state++) {
			viterbiVariables[state][0] = Double.NEGATIVE_INFINITY;
		}
		/*tAfter = System.currentTimeMillis();
		System.out.println("Initialisation took " + (tAfter - tBefore) + "ms");*/
		
		// 'Recursion'
		for(int l = 1; l < sequence.length() + 1; l++) {
			/*tBefore = System.currentTimeMillis();
			if (l % 100 == 0) {
				System.out.println("0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22");
			}*/
			for(int q = 0; q < stateCount; q++) {
				/*tBeforeFindingMax = System.currentTimeMillis();*/
				if(q == 1)
					continue;
				double max = Double.NEGATIVE_INFINITY;
				for (int qPrime = 0; qPrime < stateCount; qPrime++) {
					if(model.getLogTransitionProbability(qPrime, q) > Double.NEGATIVE_INFINITY && qPrime != 1) {
						if(qPrime == 0 && l < 10) {
							// this causes slowdown: for large L, no point in doing this;
							max = Math.max(max,
									viterbiVariables[qPrime][0] + model.getLogTransitionProbability(qPrime, q)
											+ model.getLogEmissionProbability(0, q, "",
													sequence.substring(0, l)));
						} else {
							// q' \in Q, i.e. neither initial nor terminal state (cannot come from the terminal state)
							// for(int lPrime = Math.max(0, l - 100); lPrime < l; lPrime++) {
							for(int lPrime : model.getState(q).iteratePermissibleLengths(l)) {
								max = Math.max(max,
										viterbiVariables[qPrime][lPrime] + model.getLogTransitionProbability(qPrime, q)
												+ model.getLogEmissionProbability(qPrime, q,
														abbreviating ? null : sequence.substring(0, lPrime),
														sequence.substring(lPrime, l)));
							}
						}
					}
				}
				/*tAfterFindingMax = System.currentTimeMillis();
				if (l % 100 == 0) 
					System.out.print((tAfterFindingMax - tBeforeFindingMax) + "  ");
				*/
				
				viterbiVariables[q][l] = max;
			}
			if (l % 1000 == 0) {
				/* tAfter = System.currentTimeMillis(); */
				System.out.println(" done l=" + l + "/" + (sequence.length() + 1)); // + "\ttook "  + (tAfter - tBefore) + "ms");
			}
				
		}
		
		System.out.println("Computed variables!");
	}
	
	public List<List<Pair<HMMState, Integer>>> computeParses() {
		computeViterbiVariables();
		ArrayList<ViterbiSeed> workLoad = new ArrayList<ViterbiSeed>();
		ArrayList<ViterbiSeed> finished = new ArrayList<ViterbiSeed>(); // here put all seeds that reach l=0
		
		double initialMax = Double.NEGATIVE_INFINITY;
		for(int q1 = 2; q1 < model.getNumberOfStates(); q1++) {
			double candidate = viterbiVariables[q1][sequence.length()] + model.getLogTransitionProbability(q1, 1);
			if(initialMax <= candidate) {
				if(initialMax < candidate) {
					initialMax = candidate;
					workLoad.clear();
				}
				workLoad.add(new ViterbiSeed(model, sequence, q1, sequence.length(), null, abbreviating)); // the ends of the parses point to null.
			}
		}
		
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
		
		ArrayList<List<Pair<HMMState, Integer>>> parses = new ArrayList<List<Pair<HMMState, Integer>>>();
		for(ViterbiSeed startSeed : finished) {
			ArrayList<Pair<HMMState, Integer>> parse = new ArrayList<Pair<HMMState, Integer>>();
			ViterbiSeed current = startSeed;
			while(current.getPrevious() != null) {
				parse.add(new Pair<HMMState, Integer>(model.getState(current.getPrevious().getQ()),
						current.getPrevious().getL() - current.getL()));
				current = current.getPrevious();
			}
			
			parses.add(parse);
		}
		
		return parses;
	}
	
}
