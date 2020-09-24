package de.vetter.pogigwasc;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
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
				new File("resources//de//vetter//pogigwasc//parameter//parameters-examplefile.properties"));
		mp = new ModelParameters(reader);
	}
	
	@Test
	public void testInvalidParameterFile() throws FileNotFoundException {
		FileReader reader = new FileReader(
				new File("resources//de//vetter//pogigwasc//parameter//parameters-invalid.properties"));
		
		try {
			mp = new ModelParameters(reader);
		} catch (IOException e) {
			assertTrue(e.getMessage().contains("base_frequencies_NCS"));
		}
	}
	
	@Test
	public void testInvalidCodonCorrection() throws FileNotFoundException {
		FileReader reader = new FileReader(
				new File("resources//de//vetter//pogigwasc//parameter//parameters-invalid-codon-correction.properties"));
		
		try {
			mp = new ModelParameters(reader);
		} catch (IOException e) {
			System.out.println(e.getMessage());
			assertTrue(e.getMessage().contains("explicit_codon_probabilities_stop_region"));
			assertTrue(e.getMessage().contains("MESS"));
		}
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
	}

	@Test
	public void testGetLogBaseProbabilityCDS() {
		assertEquals(0.24, Math.exp(mp.getLogBaseProbabilityCDS('T', 0)), TOLERANCE);
		assertEquals(0.13, Math.exp(mp.getLogBaseProbabilityCDS('C', 0)), TOLERANCE);
		assertEquals(0.31, Math.exp(mp.getLogBaseProbabilityCDS('A', 0)), TOLERANCE);
		assertEquals(0.32, Math.exp(mp.getLogBaseProbabilityCDS('G', 0)), TOLERANCE);
		
		assertEquals(0.30, Math.exp(mp.getLogBaseProbabilityCDS('T', 1)), TOLERANCE);
		assertEquals(0.18, Math.exp(mp.getLogBaseProbabilityCDS('C', 1)), TOLERANCE);
		assertEquals(0.37, Math.exp(mp.getLogBaseProbabilityCDS('A', 1)), TOLERANCE);
		assertEquals(0.15, Math.exp(mp.getLogBaseProbabilityCDS('G', 1)), TOLERANCE);
		
		assertEquals(0.36, Math.exp(mp.getLogBaseProbabilityCDS('T', 2)), TOLERANCE);
		assertEquals(0.12, Math.exp(mp.getLogBaseProbabilityCDS('C', 2)), TOLERANCE);
		assertEquals(0.38, Math.exp(mp.getLogBaseProbabilityCDS('A', 2)), TOLERANCE);
		assertEquals(0.14, Math.exp(mp.getLogBaseProbabilityCDS('G', 2)), TOLERANCE);
	}
	
	@Test(expected = IndexOutOfBoundsException.class)
	public void testGetLogBaseProbabilityCDSOutOfCodon() {
		assertEquals(0.30, Math.exp(mp.getLogBaseProbabilityCDS('T', 3)), TOLERANCE);
	}
	
	@Test(expected = IndexOutOfBoundsException.class)
	public void testGetLogBaseProbabilityCDSNotABase() {
		assertEquals(0.30, Math.exp(mp.getLogBaseProbabilityCDS('-', 1)), TOLERANCE);
	}

	@Test
	public void testGetLogCodonProbabilityStopRegion() {
		// check the four specified codons:
		assertEquals(0, Math.exp(mp.getLogCodonProbabilityStopRegion("TGA")), TOLERANCE);
		assertEquals(0.0018, Math.exp(mp.getLogCodonProbabilityStopRegion("TAA")), TOLERANCE);
		assertEquals(0.0472, Math.exp(mp.getLogCodonProbabilityStopRegion("AAG")), TOLERANCE);
		assertEquals(0.063968, Math.exp(mp.getLogCodonProbabilityStopRegion("AAA")), TOLERANCE);
		
		// Check SOME other codons:
		assertEquals(0.14*0.36*0.14, Math.exp(mp.getLogCodonProbabilityStopRegion("CAG")), TOLERANCE);
		assertEquals(0.355*0.16*0.12, Math.exp(mp.getLogCodonProbabilityStopRegion("ACC")), TOLERANCE);
		assertEquals(0.24*0.33*0.36, Math.exp(mp.getLogCodonProbabilityStopRegion("TTT")), TOLERANCE);
		assertEquals(0.265*0.16*0.14, Math.exp(mp.getLogCodonProbabilityStopRegion("GCG")), TOLERANCE);
	}
	
	@Test(expected = IllegalArgumentException.class)
	public void testGetLogCodonProbabilityStopRegionCodonTooLong() {
		assertEquals(0, Math.exp(mp.getLogCodonProbabilityStopRegion("GAATTC")), TOLERANCE);
	}
	
	@Test(expected = IllegalArgumentException.class)
	public void testGetLogCodonProbabilityStopRegionCodonTooShort() {
		assertEquals(0, Math.exp(mp.getLogCodonProbabilityStopRegion("GT")), TOLERANCE);
	}

	@Test
	public void testGetLogBaseProbabilityNCS() {
		double[] baseFrequencies = {0.46, 0.10, 0.355, 0.085};
		String bases = "TCAG";
		for(int i = 0; i < 4; i++) {
			assertEquals(baseFrequencies[i], Math.exp(mp.getLogBaseProbabilityNCS(bases.charAt(i))), TOLERANCE);	
		}
	}

	@Test
	public void testGetLogBaseProbabilityIntron() {
		double[] baseFrequencies = {0.46, 0.10, 0.355, 0.085};
		String bases = "TCAG";
		for(int i = 0; i < 4; i++) {
			assertEquals(baseFrequencies[i], Math.exp(mp.getLogBaseProbabilityIntron(bases.charAt(i))), TOLERANCE);	
		}
	}

	@Test
	public void testGetLogBaseProbabilityStartRegion() {
		double[] baseFrequencies = {0.15, 0.02, 0.78, 0.05};
		String bases = "TCAG";
		for(int i = 0; i < 4; i++) {
			assertEquals(baseFrequencies[i], Math.exp(mp.getLogBaseProbabilityStartRegion(bases.charAt(i))), TOLERANCE);
		}
	}

	@Test
	public void testGetLogBaseProbabilitySDS() {
		double[][] baseFrequencies = { { 0, 0, 0, 1 }, { 1, 0, 0, 0 }, { 0.02, 0.02, 0.94, 0.02 },
				{ 0.02, 0.02, 0.94, 0.02 }, { 0.35, 0.05, 0.15, 0.45 } };
		String bases = "TCAG";
		for (int pos = 0; pos < 5; pos++) {
			for (int i = 0; i < 4; i++) {
				assertEquals(baseFrequencies[pos][i], Math.exp(mp.getLogBaseProbabilitySDS(bases.charAt(i), pos)),
						TOLERANCE);
			}
		}
	}

	@Test
	public void testGetLogBaseProbabilitySAS() {
		double[][] baseFrequencies = { { 0.45, 0.05, 0.35, 0.15 }, { 0.65, 0.25, 0.05, 0.05 }, { 0, 0, 1, 0 },
				{ 0, 0, 0, 1 } };
		String bases = "TCAG";
		for (int pos = 0; pos < 4; pos++) {
			for (int i = 0; i < 4; i++) {
				assertEquals(baseFrequencies[pos][i], Math.exp(mp.getLogBaseProbabilitySAS(bases.charAt(i), pos)),
						TOLERANCE);
			}
		}
	}

	@Test
	public void testGetMaxIntronSize() {
		assertEquals(30, mp.getMaxIntronSize());
	}

	@Test
	public void testGetMinIntronSize() {
		assertEquals(12, mp.getMinIntronSize());
	}

	@Test
	public void testGetSDSSize() {
		assertEquals(5, mp.getSDSSize());
	}

	@Test
	public void testGetSASSize() {
		assertEquals(4, mp.getSASSize());
	}

}
