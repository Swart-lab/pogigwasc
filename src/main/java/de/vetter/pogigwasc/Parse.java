package de.vetter.pogigwasc;

import java.util.ArrayList;
import java.util.List;


import de.vetter.pogigwasc.states.HMMState;

/**
 * Implements the notion of a <i>parse</i> as used by Viterbi's algorithm: A
 * parse is a list/vector of pairs of a state and a length (a natural
 * number)<br>
 * 
 * It does not immediately translate into GFF-format, so this class provides the
 * translation-scheme
 * 
 * @author David Emanuel Vetter
 *
 */
public class Parse {
	
	private List<Pair<HMMState, Integer>> parse;
	
	public Parse() {
		parse = new ArrayList<Pair<HMMState, Integer>>();
	}
	
	public void add(Pair<HMMState, Integer> step) {
		parse.add(step);
	}
	
	public void add(HMMState state, int length) {
		parse.add(new Pair<HMMState, Integer>(state, length));
	}
	
	public Pair<HMMState, Integer> get(int index) {
		return parse.get(index);
	}
	
	/**
	 * @return the <i>length</i> of the parse, i.e. the sum of the lengths as listed in the pairs
	 */
	public int getLength() {
		int length = 0;
		for(Pair<HMMState, Integer> step : parse) {
			length += step.getSecond();
		}
		return length;
	}
	
	/**
	 * @return The number of steps -- this corresponds to the {@link List#size()}-method.
	 */
	public int getNumberOfSteps() {
		return parse.size();
	}
}
