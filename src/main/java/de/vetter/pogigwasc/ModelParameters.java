package de.vetter.pogigwasc;

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
	
	private int intronMin, intronMax;
	private double intronMean, intronTruncatedPoissonNormaliser;
	
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
			// old start-region (i.e. just deterministic ATG):
			// startregionSize = 3;
			
			currentParameter = "stop_region_size";
			stopregionSize = Integer.parseInt(properties.getProperty(currentParameter));

			currentParameter = "intron_minimum_length";
			intronMin = Integer.parseInt(properties.getProperty(currentParameter));
			currentParameter = "intron_maximum_length";
			intronMax = Integer.parseInt(properties.getProperty(currentParameter));
			currentParameter = "intron_mean_length";
			intronMean = Double.parseDouble(properties.getProperty(currentParameter));
			computePoissonNormaliser();

			
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
	 * Checks whether extracted probabilities are valid and consistent
	 * @param properties
	 */
	private void checkProbabilities() {
		double tolerance = 1e-6;
		
		if(getProbabilityOfStayingInNCS() < 0 || getProbabilityOfStayingInNCS() >= 1) {
			throw new IllegalArgumentException("transition_probability_of_staying_in_NCS has to be in [0,1), but was "
					+ getProbabilityOfStayingInNCS());
		}
		
		if(Math.abs(getProbabilityOfStayingInCDS() + getProbabilityGeneEnds() + 3*getProbabilityCDSToIntron() - 1.0) > tolerance) {
			throw new IllegalArgumentException("transitions leaving CDS do not sum up to 1");
		}
		
		if(!(intronMin <= intronMean && intronMean <= intronMax)) {
			throw new IllegalArgumentException("intron min, mean and max length do not obey min <= mean <= max");
		}
		
		Utilities.assertNormalisation(baseProbabilitiesNCS, tolerance, 
				"The base-frequencies for NCS are not normalised, they sum up to ");
		
		Utilities.assertNormalisation(baseProbabilitiesIntron, tolerance,
				"The base-frequencies for the middle-part of intron are not normalised, they sum up to ");
		
		Utilities.assertNormalisation(baseProbabilitiesStartregion, tolerance,
				"The base-frequencies for the start region are not normalised, they sum up to ");

		for (int row = 0; row < baseMarginalsCDS.length; row ++) {
			Utilities.assertNormalisation(baseMarginalsCDS[row], tolerance, "The base-frequencies in row " + row
					+ " of base_frequency_marginals_CDS do not sum up to 1, but to ");
		}
		
		for (int row = 0; row < baseMarginalsStopregion.length; row++) {
			Utilities.assertNormalisation(baseMarginalsStopregion[row], tolerance, "The base-frequencies in row " + row
					+ " of base_frequency_marginals_stop_region do not sum up to 1, but to ");
		}
		
		for (int row = 0; row < baseProbabilitiesSDS.length; row ++) {
			Utilities.assertNormalisation(baseProbabilitiesSDS[row], tolerance, "The base-frequencies in row " + row
					+ " of intron_base_frequencies_splice_donor_site do not sum up to 1, but to ");
		}
		
		for (int row = 0; row < baseProbabilitiesSAS.length; row ++) {
			Utilities.assertNormalisation(baseProbabilitiesSAS[row], tolerance, "The base-frequencies in row " + row
					+ " of intron_base_frequencies_splice_acceptor_site do not sum up to 1, but to ");
		}
		
		double correctedSum = 0d;
		double basewiseSum = 0d;
		for(String codon : stopRegionExplicitCodons.keySet()) {
			correctedSum += stopRegionExplicitCodons.get(codon);
			double codonProbability = 1d;
			for (int i = 0; i < 3; i++) {
				codonProbability *= baseMarginalsStopregion[i][Utilities.baseToIndex(codon.charAt(i))];
			}
			basewiseSum += codonProbability;
		}
		if(Math.abs(correctedSum - basewiseSum) > tolerance) {
			throw new IllegalArgumentException(
					"The probabilities listed in explicit_codon_probabilities_stop_region do not sum up to the same "
							+ "value as under the base-wise model: " + correctedSum + " vs " + basewiseSum);
		}
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
	 * Computes the "normaliser" (an additive constant) for the truncated poisson:
	 * Since a hard minimum and maximum are enforced, while the poisson-distribution
	 * generally allows smaller and larger values, just take the porbability of
	 * picking a value outside the allowed range, and add it uniformly to all
	 * allowed values (there are max-min+1 allowed values)
	 */
	private void computePoissonNormaliser() {
		double sumOverValid= 0;

		double lambda = intronMean - intronMin;
		for(int l = 0; l <= intronMax - intronMin; l++) {
			sumOverValid += Math.exp(l * Math.log(lambda) - Utilities.logFactorial(l));
		}

		intronTruncatedPoissonNormaliser = Math.log(sumOverValid) - lambda;
		intronTruncatedPoissonNormaliser = 1-Math.exp(intronTruncatedPoissonNormaliser);
		intronTruncatedPoissonNormaliser = intronTruncatedPoissonNormaliser/(intronMax - intronMin + 1);
		// This value will be ADDED to the probability of seeing any length within min-max.
		/*
		System.out.println("Intron-state: length-distribution:");
		System.out.println("l,poisson,empirical");

		double[] LENGTH_PROBABILITIES = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				0.0, 0.0, 0.0, 0.0, 0.2188, 0.2812, 0.1042, 0.2292, 0.1146, 0.0208, 0.0104,
				0.0104, 0.0, 0.0, 0.0104, 0.0, 0.0, 0.0 };
		
		double sum = 0;
		for(int k = 0; k < intronMax + 2; k++) {
			sum += Math.exp(getLogProbabilityIntronLength(k));
			System.out.println(k + "," + Math.exp(getLogProbabilityIntronLength(k)) + ","
					+ (k < LENGTH_PROBABILITIES.length ? LENGTH_PROBABILITIES[k] : 0));
		}
		System.out.println("\n sum = " + sum + " (should be 1)");
		*/
	}
	
	/**
	 * @param length the intron's length (including GT and AG, and all the rest of SAS and SDS)
	 * @return log probability of seeing the given intron length
	 */
	public double getLogProbabilityIntronLength(int length) {
		// Poisson is truncated!
		if(length < intronMin || intronMax < length)
			return Double.NEGATIVE_INFINITY;
		
		// to get poisson with smaller variance: shift min to 0 (truncating anyways)
		double lambda = intronMean - intronMin;
		int k = length - intronMin;
		
		// This would be the normal poisson-probability
		double poisson = k * Math.log(lambda) - lambda; 
		poisson -= Utilities.logFactorial(k);
		// However: Need/Want to add a little bit of normalisation to make up for normalisation:
		poisson = Math.exp(poisson) + intronTruncatedPoissonNormaliser;
		// this is no longer in log, so return to log:
		poisson = Math.log(poisson);
		
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
