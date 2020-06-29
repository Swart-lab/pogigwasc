package de.vetter.masterthesis;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

import de.vetter.masterthesis.states.HMMState;

public class TestViterbi {

	@Test
	public void testComputeParses() {
		GHMM ghmm = new GHMM();
		ghmm.addState(new HMMState("A") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				return newEmission == "1" ? Math.log(0.8) : Math.log(0.2);
			}
		});
		
		ghmm.addState(new HMMState("B") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				if(newEmission.contains("1"))
					return Double.NEGATIVE_INFINITY;
				if(newEmission.length() == 4)
					return Math.log(0.2);
				else if(newEmission.length() == 6)
					return Math.log(0.35);
				else if(newEmission.length() == 8)
					return Math.log(0.3);
				else if(newEmission.length() == 10)
					return Math.log(0.15);
				
				return Double.NEGATIVE_INFINITY;
			}
		});
		
		ghmm.initialiseTransitionMatrix(); 
		
		ghmm.setTransitionProbability(2, 2, 0.7);
		ghmm.setTransitionProbability(2, 3, 0.25);
		ghmm.setTransitionProbability(2, 1, 0.05);
		ghmm.setTransitionProbability(3, 2, 1);
		ghmm.setTransitionProbability(3, 1, 0);
		
		System.out.println(ghmm);
		
		Viterbi viterbi = new Viterbi(ghmm, "11011010110000000111111011");
		for(List<Pair<HMMState, Integer>> parse : viterbi.computeParses()) {
			System.out.print("\nParse:");
			for(Pair<HMMState, Integer> s : parse) {
				System.out.print(s.getFirst().getName() + ":" + s.getSecond() + " - ");
			}
		}
		// Because there is a sequence of seven 0 in the sequence: can either have the first or last be emitted by A
		assertEquals(2, viterbi.computeParses().size());
	}

}
