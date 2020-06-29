package de.vetter.masterthesis;

import java.util.ArrayList;
import java.util.List;

/**
 * Class for traceback/parse-retrieval
 * 
 * @author David Emanuel Vetter
 */
public class ViterbiSeed {
	
	private GHMM model;
	private String sequence;
	private int q, l;
	private ViterbiSeed previous;

	public ViterbiSeed(GHMM model, String sequence, int q, int l, ViterbiSeed previous) {
		this.model = model;
		this.sequence = sequence;
		this.q = q;
		this.l = l;
		this.previous = previous;
	}
	
	public ViterbiSeed getPrevious() {
		return previous;
	}
	
	public boolean isFinished() {
		return l == 0;
	}
	
	public int getQ() {
		return q;
	}
	
	public int getL() {
		return l;
	}
	
	public List<ViterbiSeed> step(double[][] viterbiVariables) {
		ArrayList<ViterbiSeed> result = new ArrayList<ViterbiSeed>();
		
		double max = Double.NEGATIVE_INFINITY;
		ArrayList<Pair<Integer, Integer>> argmaxes = new ArrayList<Pair<Integer, Integer>>();
		for (int qPrime = 0; qPrime < viterbiVariables.length; qPrime++) {
			if (qPrime == 0) {
				double candidate = viterbiVariables[qPrime][0] + model.getLogTransitionProbability(qPrime, q)
						+ model.getLogEmissionProbability(qPrime, q, "", sequence.substring(0, l));
				if (max <= candidate) {
					if (max < candidate) {
						argmaxes.clear();
						max = candidate;
					}
					argmaxes.add(new Pair<Integer, Integer>(0, 0));
				}
			} else {
				for(int lPrime = 1; lPrime < l; lPrime++) {
					double candidate = viterbiVariables[qPrime][lPrime] + model.getLogTransitionProbability(qPrime, q)
							+ model.getLogEmissionProbability(qPrime, q, sequence.substring(0, lPrime),
									sequence.substring(lPrime, l));
					if(max <= candidate) {
						if (max < candidate) {
							argmaxes.clear();
							max = candidate;
						}
						argmaxes.add(new Pair<Integer, Integer>(qPrime, lPrime));
					}
				}
			}
		}

		for(Pair<Integer, Integer> p : argmaxes) {
			result.add(new ViterbiSeed(model, sequence, p.getFirst(), p.getSecond(), this));
		}
		
		return result;
	}
	
}
