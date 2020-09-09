package de.vetter.pogigwasc;

import static org.junit.Assert.*;

import org.junit.Test;

import de.vetter.pogigwasc.states.FixedSequenceState;
import de.vetter.pogigwasc.states.HMMState;

public class TestParse {

	@Test
	public void testAddPairOfHMMStateInteger() {
		Parse parse = new Parse();
		HMMState dummy = new FixedSequenceState("name", "ATG");
		
		assertEquals(0, parse.getNumberOfSteps());
		parse.add(new Pair<HMMState, Integer>(dummy, 3));
		assertEquals(1, parse.getNumberOfSteps());
		
		Pair<HMMState, Integer> pair = parse.get(0);
		assertEquals("name", pair.getFirst().getName());
		assertEquals(0, pair.getFirst().computeLogEmissionProbability(0, null, "ATG"), 1e-6);
		assertEquals(3, pair.getSecond().intValue());
	}

	@Test
	public void testAddHMMStateIntAndGet() {
		Parse parse = new Parse();
		parse.add(new FixedSequenceState("first", "ATG"), 3);
		parse.add(new FixedSequenceState("second", "AAAA"), 4);
		parse.add(new FixedSequenceState("third", "G"), 1);
		
		assertEquals("first", parse.get(0).getFirst().getName());
		assertEquals("second", parse.get(1).getFirst().getName());
		assertEquals("third", parse.get(2).getFirst().getName());
		
		assertEquals(3, parse.get(0).getSecond().intValue());
		assertEquals(4, parse.get(1).getSecond().intValue());
		assertEquals(1, parse.get(2).getSecond().intValue());
	}

	@Test
	public void testGetLengthEmpty() {
		Parse parse = new Parse();
		assertEquals(0, parse.getLength());
	}
	
	@Test
	public void testGetLength() {
		Parse parse = new Parse();
		parse.add(new FixedSequenceState("first", "ATG"), 3);
		assertEquals(3, parse.getLength());
		parse.add(new FixedSequenceState("second", "CTCTC"), 5);
		assertEquals(8, parse.getLength());
		parse.add(new FixedSequenceState("third", "GG"), 2);
		assertEquals(10, parse.getLength());
		parse.add(new FixedSequenceState("fourth", "TGA"), 3);
		assertEquals(13, parse.getLength());
	}

	@Test
	public void testGetNumberOfStepsEmpty() {
		Parse parse = new Parse();
		assertEquals(0, parse.getNumberOfSteps());
	}
	
	@Test
	public void testGetNumberOfStepsMultiple() {
		Parse parse = new Parse();
		parse.add(new FixedSequenceState("first", "ATG"), 3);
		assertEquals(1, parse.getNumberOfSteps());
		parse.add(new FixedSequenceState("second", "CTCTC"), 5);
		assertEquals(2, parse.getNumberOfSteps());
		parse.add(new FixedSequenceState("third", "GG"), 2);
		assertEquals(3, parse.getNumberOfSteps());
		parse.add(new FixedSequenceState("fourth", "TGA"), 3);
		assertEquals(4, parse.getNumberOfSteps());
	}
}
