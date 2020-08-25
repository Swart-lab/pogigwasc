package de.vetter.masterthesis;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

/**
 * Class for traceback/parse-retrieval: Key-method is {@link #step(double[][])},
 * which will 'spawn' a number (at least 1) of Viterbi-Seeds that have regressed
 * backwards into the sequence: To compute the viterbi-parses, step the initial
 * seeds through the matrix, until they reach the start of the sequence<br>
 * Notice that this produces a very specialized list (or, since there may be
 * multiple ViterbiSeeds pointing at the same previous seed, a tree), with the
 * {@link #getPrevious()} comparable to the {@link Iterator#next()} method.
 * 
 * @author David Emanuel Vetter
 */
public class ViterbiSeed {
	
	private GHMM model;
	private String sequence;
	private boolean abbreviating = false;
	private int q, l;
	private ViterbiSeed previous;

	/**
	 * @param model      the gHMM which is assumed to have produced the sequence
	 * @param sequence   the total sequence emitted by that gHMM (could also pass
	 *                   this just to step)
	 * @param q          current state
	 * @param l          current position in the sequence
	 * @param previous   the parent-seed that spawned this one; may be {@code null}
	 *                   (i.e. this seed is a final seed)
	 * @param abbreviate whether to abbreviate the computation of the
	 *                   emission-probabilities by omitting the emission history,
	 *                   cf. {@link Viterbi#setAbbreviating(boolean)}; defaults to
	 *                   {@code false}
	 */
	public ViterbiSeed(GHMM model, String sequence, int q, int l, ViterbiSeed previous, boolean abbreviate) {
		this(model, sequence, q, l, previous);
		this.abbreviating = abbreviate;
	}
	
	/**
	 * Use default-value for abbreviating, otherwise same as {@link #ViterbiSeed(GHMM, String, int, int, ViterbiSeed, boolean)}
	 * @param model
	 * @param sequence
	 * @param q
	 * @param l
	 * @param previous
	 */
	public ViterbiSeed(GHMM model, String sequence, int q, int l, ViterbiSeed previous) {
		this.model = model;
		this.sequence = sequence;
		this.q = q;
		this.l = l;
		this.previous = previous;
	}
	
	/**
	 * @return the {@link ViterbiSeed} that spawned this one; think of this like the {@link Iterator#next()} method
	 */
	public ViterbiSeed getPrevious() {
		return previous;
	}
	
	/**
	 * @return Whether this Seed has reached the start of the sequence -- trying to
	 *         call {@link #step(double[][])} will then no longer yield results of
	 *         any value
	 */
	public boolean isFinished() {
		return l == 0;
	}
	
	/**
	 * @return q, the state corresponding to this seed/step
	 */
	public int getQ() {
		return q;
	}
	
	/**
	 * @return l, the position in the sequence corresponding to this seed/step
	 */
	public int getL() {
		return l;
	}
	
	/**
	 * Does a single step of traceback: Will collect all MAP ViterbiSeeds
	 * (corresponding to the pairs (q, l) of states and sequence-positions, whence
	 * the current entry in the viterbi-matrix most probably came). <br>
	 * For the intended usage, refer to {@link Viterbi#computeParses()}
	 * 
	 * @param viterbiVariables
	 * @return a list, containing equally probable seeds
	 */
	public List<ViterbiSeed> step(double[][] viterbiVariables) {
		ArrayList<ViterbiSeed> result = new ArrayList<ViterbiSeed>();
		
		double max = Double.NEGATIVE_INFINITY;
		// here need same optimisations as in Viterbi!
		ArrayList<Pair<Integer, Integer>> argmaxes = new ArrayList<Pair<Integer, Integer>>();
		for (int qPrime = 0; qPrime < viterbiVariables.length; qPrime++) {
			if(qPrime == 1) {
				// cannot leave terminal state, so no need to check there.
				continue;
			}
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
				for(int lPrime : model.getState(q).iteratePermissibleLPrimes(l)) {
					double candidate = viterbiVariables[qPrime][lPrime] + model.getLogTransitionProbability(qPrime, q)
							+ model.getLogEmissionProbability(qPrime, q, 
									abbreviating ? null : sequence.substring(0, lPrime),
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
			result.add(new ViterbiSeed(model, sequence, p.getFirst(), p.getSecond(), this, abbreviating));
		}
		
		return result;
	}
	
}
