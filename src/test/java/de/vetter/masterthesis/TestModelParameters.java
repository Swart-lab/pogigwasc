package de.vetter.masterthesis;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.junit.Before;
import org.junit.Test;

public class TestModelParameters {

	ModelParameters mp;
	private static final double TOLERANCE = 1e-10;
	
	// could reasonably switch this to BeforeClass
	@Before
	public void setup() throws IOException {
		FileReader reader = new FileReader(
				new File("resources//de//vetter//masterthesis//parameter//parameters-examplefile.properties"));
		mp = new ModelParameters(reader);
	}
	
	@Test
	public void testParseMatrix() {
		try {
			double[][] matrix = ModelParameters.parseMatrix("{   {3, 1, - 5}  {2.37, 2, 0}}", 3);
			assertEquals(3d, matrix[0][0], TOLERANCE);
			assertEquals(1d, matrix[0][1], TOLERANCE);
			assertEquals(-5d, matrix[0][2], TOLERANCE);
			assertEquals(2.37, matrix[1][0], TOLERANCE);
			assertEquals(2d, matrix[1][1], TOLERANCE);
			assertEquals(0d, matrix[1][2], TOLERANCE);
			assertEquals(2, matrix.length);
			assertEquals(3, matrix[0].length);
			assertEquals(3, matrix[1].length);
		} catch (IllegalArgumentException e) {
			fail("Should not cause an exception");
			e.printStackTrace();
		}
	}
	
	@Test(expected=IllegalArgumentException.class)
	public void testParseMatrixInvalid() {
		double[][] matrix = ModelParameters.parseMatrix("{   {3, 1, - 5}  {2.37, 2}}", 3);
		matrix[0] = matrix[1];
	}
	
	@Test(expected=IllegalArgumentException.class)
	public void testParseMatrixInvalidLengthButSquareMatrix() {
		double[][] matrix = ModelParameters.parseMatrix("{   {3, 1}  {2.37, 2}}", 3); // Length does not match matrix-length
		matrix[0] = matrix[1];
	}
	
	@Test
	public void testParseRowVector() {
		double[] array = ModelParameters.parseRowVector("{-1.034, 2. 7, 8, -0, .1}", 5);
		assertEquals(5, array.length);
		assertEquals(-1.034, array[0], TOLERANCE);
		assertEquals(2.7, array[1], TOLERANCE); // NOTE the tolerance for interrupting whitespaces-- the commata already delimit the entries
		assertEquals(8d, array[2], TOLERANCE);
		assertEquals(0d, array[3], TOLERANCE);
		assertEquals(0.1, array[4], TOLERANCE);
	}
	
	@Test(expected=IllegalArgumentException.class)
	public void testParseRowVectorInvalidLength() {
		double[] array = ModelParameters.parseRowVector("{-1.034, 2. 7, 8, -0, .1}", 3);
		assertEquals(5, array.length);
	}
	
	@Test(expected=NumberFormatException.class)
	public void testParseRowVectorInvalidContent() {
		double[] array = ModelParameters.parseRowVector("{-1.034, , a}", 3);
		assertEquals(5, array.length);
	}

	@Test
	public void testGetProbabilityOfStayingInNCS() {
		assertEquals(0.997588, mp.getProbabilityOfStayingInNCS(), TOLERANCE);
	}

	@Test
	public void testGetProbabilityOfStayingInCDS() {
		assertEquals(0.99995, mp.getProbabilityOfStayingInCDS(), TOLERANCE);
	}

	@Test
	public void testGetProbabilityGeneEnds() {
		assertEquals(0.00005 * 0.62, mp.getProbabilityGeneEnds(), TOLERANCE);
	}

	@Test
	public void testGetProbabilityCDSToIntron() {
		assertEquals(0.00005 * 0.38 / 3d, mp.getProbabilityCDSToIntron(), TOLERANCE);
	}

	@Test
	public void testGetStartRegionSize() {
		assertEquals(6, mp.getStartRegionSize());
	}

	@Test
	public void testGetStopRegionSize() {
		assertEquals(24, mp.getStopRegionSize());
	}

	@Test
	public void testGetLogProbabilityIntronLength() {
		assertEquals(Double.NEGATIVE_INFINITY, mp.getLogProbabilityIntronLength(11), TOLERANCE);
		assertEquals(Double.NEGATIVE_INFINITY, mp.getLogProbabilityIntronLength(31), TOLERANCE);
		assertTrue(Double.NEGATIVE_INFINITY < mp.getLogProbabilityIntronLength(12));
		assertTrue(Double.NEGATIVE_INFINITY < mp.getLogProbabilityIntronLength(30));
		assertEquals(mp.getLogProbabilityIntronLength(17), 0.55*mp.getLogProbabilityIntronLength(18), TOLERANCE);
		// TODO: might actually check the entire distribution
	}

	@Test
	public void testGetLogBaseProbabilityCDS() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetLogCodonProbabilityStopRegion() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetLogBaseProbabilityNCS() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetLogBaseProbabilityIntron() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetLogBaseProbabilityStartRegion() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetLogBaseProbabilitySDS() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetLogBaseProbabilitySAS() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetMaxIntronSize() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetMinIntronSize() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetSDSSize() {
		fail("Not yet implemented");
	}

	@Test
	public void testGetSASSize() {
		fail("Not yet implemented");
	}

}
