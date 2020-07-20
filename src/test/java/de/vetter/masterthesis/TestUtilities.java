package de.vetter.masterthesis;

import static org.junit.Assert.*;

import org.junit.Test;

public class TestUtilities {

	@Test
	public void testBaseToIndex() {
		String bases = "TCAG";
		for(int i = 0; i < 4; i++) {
			assertEquals(i, Utilities.baseToIndex(bases.charAt(i)));
			assertEquals(i, Utilities.baseToIndex(bases.toLowerCase().charAt(i)));	
		}
	}

	@Test
	public void testReverseComplement() {
		String example = "GAATTC";
		assertEquals(example, Utilities.reverseComplement(example));
		example = "TCCATGGGGA";
		assertEquals("TCCCCATGGA", Utilities.reverseComplement(example));
	}
	
	@Test
	public void testReverseComplementUnspecificSequence() {
		String example = "GTNNNANNAG";
		assertEquals("CTNNTNNNAC", Utilities.reverseComplement(example));
	}

	@Test
	public void testLogFactorial() {
		double runningSum = 0;
		for(int i = 1; i < 15; i++) {
			runningSum += Math.log(i);
			assertEquals(runningSum, Utilities.logFactorial(i), 1e-10);
		}
	}
	
	@Test
	public void testLogFactorialLarge() {
		assertEquals(148.47776695177305, Utilities.logFactorial(50), 1e-10);
	}

}
