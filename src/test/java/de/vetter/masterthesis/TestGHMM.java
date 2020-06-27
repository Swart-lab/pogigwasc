package de.vetter.masterthesis;

import static org.junit.Assert.*;

import org.junit.Test;

public class TestGHMM {

	@Test
	public void testGHMM() {
		// -> check initial and terminal state!
		GHMM ghmm = new GHMM();
		assertEquals(2, ghmm.getNumberOfStates());
		
		assertEquals(InitialState.INITIAL_STATE_NAME, ghmm.getState(0).getName());
		assertTrue(ghmm.getState(0) instanceof InitialState);
		assertEquals(TerminalState.TERMINAL_STATE_NAME, ghmm.getState(1).getName());
		assertTrue(ghmm.getState(1) instanceof TerminalState);
	}

	@Test
	public void testInitialiseTransitionMatrix() {
		GHMM ghmm = new GHMM();
		ghmm.initialiseTransitionMatrix();
		
		assertTrue(ghmm.checkTransitions());
		
		ghmm.addState(new IntronState("intron a"));
		ghmm.addState(new IntronState("intron b"));
		ghmm.initialiseTransitionMatrix();
		
		assertTrue(ghmm.checkTransitions());
		assertEquals(0, ghmm.getLogTransitionProbability(0, 2), 1e-9);
		assertEquals(0, ghmm.getLogTransitionProbability(2, 3), 1e-9);
		// from last 'proper' state to terminal
		assertEquals(0, ghmm.getLogTransitionProbability(3, 1), 1e-9);
	}

	@Test
	public void testIsSetTransitionMatrix() {
		GHMM ghmm = new GHMM();
		assertFalse(ghmm.isSetTransitionMatrix());
		
		ghmm.initialiseTransitionMatrix();
		assertTrue(ghmm.isSetTransitionMatrix());
		
		ghmm.addState(new IntronState("intron b"));
		assertFalse(ghmm.isSetTransitionMatrix());
		
		ghmm.initialiseTransitionMatrix();
		assertTrue(ghmm.isSetTransitionMatrix());
	}

	@Test
	public void testCheckTransitions() {
		GHMM ghmm = new GHMM();
		ghmm.initialiseTransitionMatrix();
		assertTrue(ghmm.checkTransitions());
		
		ghmm.addState(new IntronState("a"));
		ghmm.initialiseTransitionMatrix();
		
		assertTrue(ghmm.checkTransitions());
		
		ghmm.setTransitionProbability(2, 2, 0.3);
		assertFalse(ghmm.checkTransitions());
		
		ghmm.setTransitionProbability(2, 1, 0.7);
		assertTrue(ghmm.checkTransitions());
	}

	@Test
	public void testSetTransitionProbability() {
		GHMM ghmm = new GHMM();
		ghmm.addState(new IntronState("intron a"));
		
		ghmm.initialiseTransitionMatrix();
		assertEquals(0, Math.exp(ghmm.getLogTransitionProbability(1, 2)), 1e-9);
		
		ghmm.setTransitionProbability(2, 2, 0.6);
		ghmm.setTransitionProbability(2, 1, 0.4);
		// Notice the exponential
		assertEquals(0.6, Math.exp(ghmm.getLogTransitionProbability(2, 2)), 1e-9);
		assertEquals(0.4, Math.exp(ghmm.getLogTransitionProbability(2, 1)), 1e-9);
		// this is not to be changed:
		assertEquals(0, Math.exp(ghmm.getLogTransitionProbability(1, 2)), 1e-9);
		
		ghmm.setTransitionProbability(2, 2, 0); // notice: only doing this yields an invalid matrix.
		assertEquals(0, Math.exp(ghmm.getLogTransitionProbability(2, 2)), 1e-9);
	}

	@Test
	public void testNormaliseExitProbabilities() {
		GHMM ghmm = new GHMM();
		ghmm.addState(new IntronState("intron a"));
		ghmm.initialiseTransitionMatrix();
		
		ghmm.setTransitionProbability(2, 2, 7);
		ghmm.setTransitionProbability(2, 1, 3);
		
		assertFalse(ghmm.checkTransitions());
		ghmm.normaliseExitProbabilities(2);
		assertTrue(ghmm.checkTransitions());

		assertEquals(0.7, Math.exp(ghmm.getLogTransitionProbability(2, 2)), 1e-9);
		assertEquals(0.3, Math.exp(ghmm.getLogTransitionProbability(2, 1)), 1e-9);
	}
	
	@Test
	public void testNormaliseExitProbabilitiesAllZero() {
		GHMM ghmm = new GHMM();
		ghmm.addState(new IntronState("intron a"));
		ghmm.initialiseTransitionMatrix();
		
		ghmm.setTransitionProbability(2, 1, 0); // PROBABILITY!

		assertFalse("Transitions should no longer be valid", ghmm.checkTransitions());
		ghmm.normaliseExitProbabilities(2);
		assertTrue("Transitions should be valid again", ghmm.checkTransitions());

		assertEquals("Shoudl have introduced a forwarding-state", 1, Math.exp(ghmm.getLogTransitionProbability(2, 1)), 1e-9);
	}

	@Test
	public void testGetLogEmissionProbability() {
		GHMM ghmm = new GHMM();
		IntronState intron = new IntronState("intron a");
		ghmm.addState(intron);
		ghmm.initialiseTransitionMatrix();
		
		assertEquals(intron.computeLogEmissionProbability(0, "GGGGG", "ATACCTTAAAAC"),
				ghmm.getLogEmissionProbability(0, 2, "GGGGG", "ATACCTTAAAAC"), 1e-9);
	}

	@Test
	public void testAddState() {
		GHMM ghmm = new GHMM();
		
		assertEquals(2, ghmm.getNumberOfStates());
		
		ghmm.addState(new IntronState("intron a"));
		assertEquals(3, ghmm.getNumberOfStates());
		assertTrue(ghmm.getState(2) instanceof IntronState);
		
		assertTrue(ghmm.getState(0) instanceof InitialState);
		assertTrue(ghmm.getState(1) instanceof TerminalState);
	}

	@Test
	public void testGetState() {
		GHMM ghmm = new GHMM();
		assertTrue(ghmm.getState(0) instanceof InitialState);
		assertTrue(ghmm.getState(1) instanceof TerminalState);
		ghmm.addState(new IntronState("intron a"));
		assertTrue(ghmm.getState(2) instanceof IntronState);
		assertEquals("intron a", ghmm.getState(2).getName());
	}

	@Test(expected = IndexOutOfBoundsException.class)
	public void testRemoveStateCannotRemoveInitial() {
		GHMM ghmm = new GHMM();
		ghmm.removeState(0);
	}
	
	@Test(expected = IndexOutOfBoundsException.class)
	public void testRemoveStateCannotRemoveTerminal() {
		GHMM ghmm = new GHMM();
		ghmm.removeState(1);
	}
	
	@Test
	public void testRemoveState() {
		GHMM ghmm = new GHMM();
		ghmm.addState(new IntronState("intron a"));
		ghmm.addState(new IntronState("intron b"));
		assertEquals("intron a", ghmm.getState(2).getName());
		
		ghmm.initialiseTransitionMatrix();
		assertTrue(ghmm.isSetTransitionMatrix());
		
		assertEquals("intron a", ghmm.removeState(2).getName());
		assertFalse(ghmm.isSetTransitionMatrix());
		
		assertEquals("intron b", ghmm.getState(2).getName());
	}

}
