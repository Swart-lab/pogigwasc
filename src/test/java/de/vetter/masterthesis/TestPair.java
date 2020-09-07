package de.vetter.masterthesis;

import static org.junit.Assert.*;

import org.junit.Test;

import de.vetter.pogigwasc.Pair;

public class TestPair {

	@Test
	public void testGetSetFirst() {
		Pair<Integer, String> pair = new Pair<Integer, String>(12, "example");
		assertEquals(12, (int) pair.getFirst());
		
		pair.setFirst(-34);
		assertEquals(-34, (int) pair.getFirst());
	}

	@Test
	public void testGetSetSecond() {
		Pair<Integer, String> pair = new Pair<Integer, String>(12, "example");
		assertEquals("example", pair.getSecond());
		
		pair.setFirst(-34);
		assertEquals("example", pair.getSecond());
		
		pair.setSecond("new");
		assertEquals("new", pair.getSecond());
	}

}
