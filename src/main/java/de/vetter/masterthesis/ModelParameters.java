package de.vetter.masterthesis;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Properties;

/**
 * 
 * Class for reading and accessing model parameters from a parameter-file. This
 * file can be produced/tweaked manually or by some training-program -- which
 * need not be written by me or in java
 * 
 * @author David Emanuel Vetter
 *
 */
public class ModelParameters {
	// A list of the parameters of interest here:
	private double transitionCDS2CDS;
	private double transitionCDS2StopGivenLeavingCDS;
	
	private double transitionNCS2NCS;
	
	private int startregionSize, stopregionSize;
	
	private int intronMin, intronMax, intronB;
	private double intronS, intronR, intronP;
	private double intronLogA;
	
	private double[] baseProbabilitiesNCS, baseProbabilitiesIntron, baseProbabilitiesStartregion;
	
	/**	(pos, base) */
	private double[][] baseMarginalsCDS, baseMarginalsStopregion, baseProbabilitiesSDS, baseProbabilitiesSAS;
	
	/** Contains the PROBABILITIES (not log) of some specific codons (could, but shouldn't, be all codons) */
	private HashMap<String, Double> stopRegionExplicitCodons;
	
	/**
	 * 
	 * @param parameterFileReader a reader of a .properties-file containing all required parameters
	 * @throws IOException
	 */
	public ModelParameters(FileReader parameterFileReader) throws IOException {
		
		Properties properties = new Properties();
		properties.load(parameterFileReader);
		
		// These two are used to give somewhat informative problem-reports to the user (see the catch-response below)
		String currentParameter = "";
		String problemDetails = "";
		
		try {

			currentParameter = "transition_probability_of_staying_in_CDS";
			transitionCDS2CDS = Double.parseDouble(properties.getProperty(currentParameter));
			
			currentParameter = "transition_probability_of_CDS_to_stop_given_that_CDS_is_being_left";
			transitionCDS2StopGivenLeavingCDS = Double.parseDouble(
					properties.getProperty(currentParameter));
			
			currentParameter = "transition_probability_of_staying_in_NCS";
			transitionNCS2NCS = Double.parseDouble(properties.getProperty(currentParameter));

			currentParameter = "start_region_size";
			startregionSize = Integer.parseInt(properties.getProperty(currentParameter));
			
			currentParameter = "stop_region_size";
			stopregionSize = Integer.parseInt(properties.getProperty(currentParameter));

			currentParameter = "intron_minimum_length";
			intronMin = Integer.parseInt(properties.getProperty(currentParameter));
			currentParameter = "intron_maximum_length";
			intronMax = Integer.parseInt(properties.getProperty(currentParameter));
			currentParameter = "intron_distribution_breakpoint_b";
			intronB = Integer.parseInt(properties.getProperty(currentParameter));
			currentParameter = "intron_distribution_leftward_dropoff_s";
			intronS = Double.parseDouble(properties.getProperty(currentParameter));
			currentParameter = "intron_distribution_rightward_dropoff_r";
			intronR = Double.parseDouble(properties.getProperty(currentParameter));
			currentParameter = "intron_distribution_peak_proportion_p";
			intronP = Double.parseDouble(properties.getProperty(currentParameter));
			computeLogBaseTermForFiniteBidirectionalGeometricDistribution();

			
			currentParameter = "base_frequencies_NCS";
			problemDetails = ", malformatted or incomplete (not of length 4)";
			baseProbabilitiesNCS = parseRowVector(properties.getProperty(currentParameter), 4);

			currentParameter = "intron_base_frequencies";
			baseProbabilitiesIntron = parseRowVector(properties.getProperty(currentParameter), 4);

			currentParameter = "base_frequencies_start_region_upstream";
			baseProbabilitiesStartregion = parseRowVector(properties.getProperty(currentParameter), 4);

			currentParameter = "base_frequency_marginals_CDS";
			problemDetails = ", malformatted or incomplete (not all rows are of length 4)";
			baseMarginalsCDS = parseMatrix(properties.getProperty(currentParameter), 4);

			currentParameter = "base_frequency_marginals_stop_region";
			baseMarginalsStopregion = parseMatrix(properties.getProperty(currentParameter), 4);

			currentParameter = "intron_base_frequencies_splice_donor_site";
			baseProbabilitiesSDS = parseMatrix(properties.getProperty(currentParameter), 4);

			currentParameter = "intron_base_frequencies_splice_acceptor_site";
			baseProbabilitiesSAS = parseMatrix(properties.getProperty(currentParameter), 4);

			currentParameter = "explicit_codon_probabilities_stop_region";
			problemDetails = "";
			stopRegionExplicitCodons = new HashMap<String, Double>();
			String stopCorrections = properties.getProperty(currentParameter);
			stopCorrections.replaceAll("\\s+", "");
			stopCorrections = stopCorrections.substring(1, stopCorrections.length() - 1); // remove {}-frame
			for (String c : stopCorrections.split(",")) {
				String[] pair = c.split(":");
				if (pair.length != 2) {
					problemDetails = " or malformatted: '" + c
							+ "' cannot be analysed as a CODON:Probability-pair";
					throw new IOException(); // this will be caught by the catch just below.
				}
				stopRegionExplicitCodons.put(pair[0].replaceAll("\\s+", ""), Double.parseDouble(pair[1]));
			}
		} catch (Exception e) {
			throw new IOException(currentParameter + " is missing" + problemDetails);
		}
		
		checkProbabilities();
	}
	
	/**
	 * Decodes a matrix encoded as a string of the following format:
	 * <ul>
	 * <li>a pair of curly brackets surrounds the entire matrix</li>
	 * <li>the rows of the matrix are only separated by optional whitespaces, <b>not</b> by commata</li>
	 * <li>each row begins with {, ends with }, and has a comma-separated list of numbers in it</li>
	 * <li>all rows have the same length</li>
	 * <ul>
	 * 
	 * @param encoded the string encoding of a matrix
	 * @param rowLength the length of all the rows
	 * @return a double-array containing all the entries of the matrix at the correct positions
	 * @throws IllegalArgumentException if not all the rows have the correct length
	 */
	public static double[][] parseMatrix(String encoded, int rowLength) throws IllegalArgumentException {
		double[][] result;
		String cleaned = encoded.replaceAll("\\s+", "");
		cleaned = cleaned.substring(2, cleaned.length()-2); // remove outer { }, and first row's {, and last row's }
		String[] rows = cleaned.split("\\}\\{");
		result = new double[rows.length][rowLength];
		for(int r = 0; r < rows.length; r++) {
			result[r] = parseRowVector("{" + rows[r] + "}", rowLength);
		}
		
		return result;
	}
	
	/**
	 * Helper method to parse a String into a double-array. The string has to be
	 * correctly formatted, by containing exactly one comma-separated list of the
	 * given number of decimal-numbers, framed altogether by curly brackets
	 * 
	 * @param encoded        the encoding of the vector
	 * @param demandedLength the length of the vector
	 * @return the decoded vector
	 * @throws IllegalArgumentException if the vector does not have the specified
	 *                                  length
	 */
	 public static double[] parseRowVector(String encoded, int demandedLength) throws IllegalArgumentException {
		double[] result = null;
		String cleaned = encoded.replaceAll("\\s+", "");
		
		cleaned = cleaned.substring(1, cleaned.length() - 1); // remove { } frame
		String[] entries = cleaned.split(",");
		if(entries.length == demandedLength) {
			result = new double[demandedLength];
			for(int i = 0; i < demandedLength; i++) {
				result[i] = Double.parseDouble(entries[i]);
			}
		} else {
			throw new IllegalArgumentException("not of demanded length");
		}
		return result;
	}
	
	/**
	 * Checks whether extracted probabilities are valid
	 * @param properties
	 */
	private void checkProbabilities() {
		// TODO: implement: Check that probabilities are valid
	}
	
	
	// All the getters for the states/dependents to call and access the parameter-values
	
	/**
	 * @return transition probability NCS -> NCS
	 */
	public double getProbabilityOfStayingInNCS() {
		return transitionNCS2NCS;
	}
	
	/**
	 * @return transition probability CDS -> CDS
	 */
	public double getProbabilityOfStayingInCDS() {
		return transitionCDS2CDS;
	}
	
	/**
	 * @return Transition probability CDS -> +stop/-start
	 */
	public double getProbabilityGeneEnds() {
		return (1-transitionCDS2CDS) * transitionCDS2StopGivenLeavingCDS;
	}
	
	/**
	 * @return the probability of switching from CDS to one particular intron state
	 *         (NOT: of switching to any intron state; difference is a factor of 3)
	 */
	public double getProbabilityCDSToIntron() {
		return (1 - transitionCDS2CDS - getProbabilityGeneEnds()) / 3d;
	}
	
	/**
	 * @return size of the start-region in nt; the start-codon here is always the
	 *         last 3 bases of this start-region
	 */
	public int getStartRegionSize() {
		return startregionSize;
	}
	
	/**
	 * @return size of the stop-region in nt; the stop-codon here is always the last
	 *         3 bases of this stop-region
	 */
	public int getStopRegionSize() {
		return stopregionSize;
	}
	
	/**
	 * Requires that all intron-distribution-parameters be already available. Then
	 * computes the normaliser a (in logarithm)
	 */
	private void computeLogBaseTermForFiniteBidirectionalGeometricDistribution() {
		double leftSum = 1 - Math.pow(intronS, intronB - intronMin);
		leftSum = leftSum / (1-intronS);
		double rightSum = 1 - Math.pow(intronR, intronMax - intronB + 1);
		rightSum = rightSum / (1-intronR);
				
		intronLogA = leftSum + intronP*rightSum;
		
		intronLogA = -Math.log(intronLogA);
		
		// For plotting in latex:
		/*
		 * private double[] LENGTH_PROBABILITIES =
		 * new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 * 0.0, 0.0, 0.0, 0.0, 0.2188, 0.2812, 0.1042, 0.2292, 0.1146, 0.0208, 0.0104,
		 * 0.0104, 0.0, 0.0, 0.0104, 0.0, 0.0, 0.0 }; int lambda = 19;
		 * System.out.println("Intron-state: length-distribution:");
		 * System.out.println("l,p,poisson,empirical"); double sum = 0; for(int k = 10;
		 * k <= MAX + 1; k++) { double poisson = k * Math.log(LAMBDA) - LAMBDA; poisson
		 * -= Utilities.logFactorial(k); System.out.println(k + "," +
		 * Math.exp(logProbability(k)) + "," + Math.exp(poisson) + "," + (k <
		 * LENGTH_PROBABILITIES.length ? LENGTH_PROBABILITIES[k] : 0)); sum +=
		 * Math.exp(logProbability(k)); } System.out.println("sum=" + sum);
		 */
	}
	
	/**
	 * @param length the intron's length (including GT and AG, and all the rest of SAS and SDS)
	 * @return log probability of seeing the given intron length
	 */
	public double getLogProbabilityIntronLength(int length) {
		double result = Double.NEGATIVE_INFINITY;
		
		double poisson = length * Math.log(intronB) - intronB; 
		poisson -= Utilities.logFactorial(length);
		
		if(intronMin <= length && length < intronB)
			result = intronLogA + (intronB- length - 1) * Math.log(intronS);
		
		if(intronB <= length && length <= intronMax)
			result = intronLogA + Math.log(intronP) + (length - intronB) * Math.log(intronR);
		
		//return result;
		return poisson;
	}
	
	/**
	 * Probability of seeing given base at given position (0,1,2) in a codon in a
	 * CDS. <br>
	 * <b>Modelling Assumption/Approximation:</b> the bases in a codon are
	 * independent random variables<br>
	 * They are in fact not, but this is used to avoid overfitting to the empirical
	 * codon-distribution.
	 * 
	 * @param base one of: 'T', 'C', 'A', 'G'
	 * @param positionInCodon 0 (first base), 1, 2 (last base)
	 * @return <b>log</b> probability of seeing this base at that position
	 */
	public double getLogBaseProbabilityCDS(char base, int positionInCodon) {
		return Math.log(baseMarginalsCDS[positionInCodon][Utilities.baseToIndex(base)]);
	}
	
	/**
	 * @param codon a string of length 3 
	 * @return log probability of seeing this codon in the stop region
	 * @throws IllegalArgumentException if the codon has a length other than 3
	 */
	public double getLogCodonProbabilityStopRegion(String codon) {
		if(codon.length() != 3) 
			throw new IllegalArgumentException("Given codon '" + codon + "' is not of length 3");
		
		double result = 0d;
		// Special cases: Determined from StopRegions.ipynb and prior: UGA, UAA are
		// depleted, use AAG and AAA for correcting this in the distribution
		
		if(stopRegionExplicitCodons.containsKey(codon)) {
			return Math.log(stopRegionExplicitCodons.get(codon));
		}
		
		// for all the rest: use base-wise (and rounded) approximation to avoid
		// overfitting
		for (int i = 0; i < 3; i++) {
			result += Math.log(baseMarginalsStopregion[i][Utilities.baseToIndex(codon.charAt(i))]);
		}
		
		return result;
	}
	
	/**
	 * @param base any of 'T', 'C', 'A', 'G'
	 * @return <b>log</b> probability of seeing that base in the NCS
	 */
	public double getLogBaseProbabilityNCS(char base) {
		return Math.log(baseProbabilitiesNCS[Utilities.baseToIndex(base)]);
	}
	
	/**
	 * @param base any of 'T', 'C', 'A', 'G'
	 * @return <b>log</b> probability of seeing that base in the middle part of an
	 *         intron (outside the SAS and SDS)
	 */
	public double getLogBaseProbabilityIntron(char base) {
		return Math.log(baseProbabilitiesIntron[Utilities.baseToIndex(base)]);
	}
	
	/**
	 * @param base any of 'T', 'C', 'A', 'G'
	 * @return <b>log</b> probability of seeing that base in the start-region
	 *         upstream of the start-AUG. Because of the particular
	 *         Kozak-Consensus-sequence (AAA<u>ATG</u>) found in ciliates (Loxodes in particular),
	 *         one categorical distribuiton suffices for the three bases upstream of
	 *         the start-AUG
	 */
	public double getLogBaseProbabilityStartRegion(char base) {
		return Math.log(baseProbabilitiesStartregion[Utilities.baseToIndex(base)]);
	}
	
	/**
	 * @param base any of 'T', 'C', 'A', 'G'
	 * @param i    index within SDS
	 * @return <b>log</b> probability of seeing base at index i in the SDS (splice
	 *         donor-site)
	 */
	public double getLogBaseProbabilitySDS(char base, int i) {
		return Math.log(baseProbabilitiesSDS[i][Utilities.baseToIndex(base)]);
	}
	
	/**
	 * @param base any of 'T', 'C', 'A', 'G'
	 * @param i    index within SAS -- not within the intron; transform index w.r.t.
	 *             intron-start into SAS-index by: index - (intron_length -
	 *             sas_size)
	 * @return <b>log</b> probability of seeing base at index i in the SAS (splice
	 *         acceptor-size)
	 */
	public double getLogBaseProbabilitySAS(char base, int i) {
		return Math.log(baseProbabilitiesSAS[i][Utilities.baseToIndex(base)]);
	}
	
	/**
	 * @return maximum allowed intron size (in nt)
	 */
	public int getMaxIntronSize() {
		return intronMax;
	}
	
	/**
	 * @return minimum allowed intron size (in nt)
	 */
	public int getMinIntronSize() {
		return intronMin;
	}
	
	/**
	 * @return size of the splice donor site (GT, potentially followed by some other
	 *         conserved sequence)
	 */
	public int getSDSSize() {
		return baseProbabilitiesSDS.length;
	}

	/**
	 * @return size of the splice acceptor site (AG, potentially preceded by some
	 *         other conserved sequence)
	 */
	public int getSASSize() {
		return baseProbabilitiesSAS.length;
	}
} 
