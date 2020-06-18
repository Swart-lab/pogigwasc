package de.vetter.simple;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Pattern;

public class ReplaceUGA {

	private static Pattern standardCodon = Pattern.compile("[TGAC]{3}?");

	/** UTR/IGR-model: empirical low-CG */
	private static double[] baseFrequenciesUTR = { 0.46, 0.1, 0.36, 0.08 };

	/** CDS-base-frequency-model */
	private static double[] baseFrequenciesCDS = { 0.3, 0.14, 0.36, 0.34 };

	private static double[] codonFrequenciesRestricted = { 0.03759398496240601, 0.02593984962406015,
			0.02725563909774436, 0.01992481203007519, 0.018796992481203006, 0.0031954887218045114, 0.014661654135338346,
			0.004511278195488722, 0.02518796992481203, 0.012406015037593985, 0.002255639097744361, 0.014661654135338346,
			0.010902255639097745, 0.007894736842105263, 0.0, 0.0231203007518797, 0.025, 0.008270676691729323,
			0.01992481203007519, 0.005827067669172932, 0.015977443609022556, 0.003947368421052632, 0.012030075187969926,
			0.0011278195488721805, 0.010338345864661654, 0.006766917293233083, 0.017481203007518795,
			0.008082706766917294, 0.004511278195488722, 0.0018796992481203006, 0.00037593984962406017,
			0.0005639097744360903, 0.03176691729323308, 0.017481203007518795, 0.02462406015037594, 0.02030075187969925,
			0.014285714285714285, 0.005639097744360902, 0.011842105263157895, 0.002819548872180451,
			0.031390977443609025, 0.016541353383458645, 0.06259398496240602, 0.04981203007518797, 0.014097744360902255,
			0.008270676691729323, 0.04003759398496241, 0.005827067669172932, 0.027819548872180452, 0.005827067669172932,
			0.016917293233082706, 0.011090225563909775, 0.02612781954887218, 0.007142857142857143, 0.014285714285714285,
			0.0031954887218045114, 0.02744360902255639, 0.011842105263157895, 0.04398496240601504, 0.014849624060150377,
			0.010526315789473684, 0.005639097744360902, 0.02387218045112782, 0.0016917293233082707 };
	private static double[] codonFrequenciesGeneral = { 0.03107769423558897, 0.016791979949874688, 0.028446115288220553,
			0.015977443609022556, 0.014223057644110276, 0.0030075187969924814, 0.013345864661654135,
			0.002380952380952381, 0.02769423558897243, 0.013095238095238096, 0.011654135338345865, 0.015977443609022556,
			0.009774436090225564, 0.008270676691729323, 0.0011278195488721805, 0.013533834586466165,
			0.017669172932330827, 0.004761904761904762, 0.014411027568922305, 0.004260651629072682,
			0.015413533834586466, 0.0018170426065162908, 0.015162907268170427, 0.0010025062656641604,
			0.014536340852130326, 0.0061403508771929825, 0.017731829573934838, 0.0057644110275689225,
			0.002443609022556391, 0.0008145363408521303, 0.0020676691729323306, 0.00018796992481203009,
			0.035213032581453634, 0.009273182957393484, 0.026566416040100252, 0.02280701754385965, 0.02218045112781955,
			0.004323308270676692, 0.017606516290726817, 0.0033208020050125315, 0.03790726817042606,
			0.017606516290726817, 0.05369674185463659, 0.027192982456140352, 0.020050125313283207, 0.008333333333333333,
			0.032518796992481204, 0.004323308270676692, 0.02268170426065163, 0.004448621553884711, 0.02149122807017544,
			0.007142857142857143, 0.02362155388471178, 0.005200501253132832, 0.02274436090225564, 0.002380952380952381,
			0.04317042606516291, 0.014849624060150377, 0.05946115288220551, 0.017481203007518795, 0.014160401002506266,
			0.0053884711779448626, 0.038784461152882206, 0.0035087719298245615 };

	public static int baseToIndex(char base) {
		switch (base) {
		case 'T':
			return 0;
		case 'C':
			return 1;
		case 'A':
			return 2;
		case 'G':
			return 3;
		default:
			System.out.println("Encountered nonstandard: " + base);
			return -1;
		}
	}

	public static int codonToIndex(String codon) {
		int index = baseToIndex(codon.charAt(0)) * 16;
		index += baseToIndex(codon.charAt(1)) * 4;
		index += baseToIndex(codon.charAt(2));
		return index;
	}

	/**
	 * @param codon
	 * @return log Probability of observing the given base outside a CDS (in particular: in 3'UTR)
	 */
	public static double getLogProbabilityBaseUTR(char base) {
		if (!"TGAC".contains("" + base)) {
			System.out.println("\tNonstandard base " + base);
			return 0; // probability of encountering any base at all is 1
		}
		return Math.log(baseFrequenciesUTR[baseToIndex(base)]);
	}

	/**
	 * @param codon
	 * @return log Probability of observing the given base in a CDS
	 */
	public static double getLogProbabilityBaseCDS(char base) {
		if (!"TGAC".contains("" + base)) {
			System.out.println("\tNonstandard base " + base);
			return 0; // probability of encountering any base at all is 1
		}
		return Math.log(baseFrequenciesCDS[baseToIndex(base)]);
	}

	/**
	 * @param codon
	 * @return log Probability of observing the given codon in the CDS right before
	 *         the stop-UGA, where no more UGA may occur
	 */
	public static double getLogProbabilityRestricted(String codon) {
		if (!standardCodon.matcher(codon).matches()) {
			// Should/Could marginalise here!
			// This currently is a fallback model
			System.out.println("\tNonstandard codon " + codon);
			return getLogProbabilityBaseCDS(codon.charAt(0)) + getLogProbabilityBaseCDS(codon.charAt(1))
					+ getLogProbabilityBaseCDS(codon.charAt(2));
		}
		return Math.log(codonFrequenciesRestricted[codonToIndex(codon)]);
	}

	/**
	 * @param codon
	 * @return log Probability of observing the given codon in a CDS
	 */
	public static double getLogProbabilityGeneral(String codon) {
		if (!standardCodon.matcher(codon).matches()) {
			// TODO ^
			System.out.println("\tNonstandard codon " + codon);
			return getLogProbabilityBaseCDS(codon.charAt(0)) + getLogProbabilityBaseCDS(codon.charAt(1))
					+ getLogProbabilityBaseCDS(codon.charAt(2));
		}
		return Math.log(codonFrequenciesGeneral[codonToIndex(codon)]);
	}

	/**
	 * @param sequence
	 * @return reverse complement of the given sequence (standard biology)
	 */
	public static String reverseComplement(String sequence) {
		StringBuffer result = new StringBuffer();
		for (char c : sequence.toCharArray()) {
			switch (c) {
			case 'T':
				result.append('A');
				break;
			case 'C':
				result.append('G');
				break;
			case 'A':
				result.append('T');
				break;
			case 'G':
				result.append('C');
				break;
			default:
				result.append(c);
			}
		}
		return result.reverse().toString();
	}

	/**
	 * @param elements
	 * @return Sum of the elements in the array
	 */
	public static double sum(double[] elements) {
		double result = 0;
		for (double e : elements) {
			result += e;
		}
		return result;
	}

	/**
	 * Noncoding-model probability: what is the log of the probability of observing
	 * the given sequences before and after a TGA, assuming that the TGA is neither
	 * stop nor sense, but outside a coding region (e.g. intergenic or 3'UTR=
	 * 
	 * @param before
	 * @param after
	 * @return log probability
	 */
	public static double noncoding(String before, String after) {
		double logProbability = 0;
		for (char c : before.toCharArray()) {
			logProbability += getLogProbabilityBaseUTR(c);
		}
		for (char c : after.toCharArray()) {
			logProbability += getLogProbabilityBaseUTR(c);
		}
		return logProbability;
	}

	/**
	 * Sense-model probability: what is the log of the probability of observing the
	 * given sequences before and after a TGA, assuming that the TGA is a sense-TGA
	 * (coding for W)
	 * 
	 * @param before
	 * @param after
	 * @return log probability
	 */
	public static double sense(String before, String after) {
		double logProbability = 0;
		/*
		 * for(char c : before.toCharArray()) { logProbability +=
		 * getLogProbabilityBaseCDS(c); } for(char c : after.toCharArray()) {
		 * logProbability += getLogProbabilityBaseCDS(c); }
		 */

		for (int i = 0; i < before.length() - 3; i += 3) {
			logProbability += getLogProbabilityGeneral(before.substring(i, i + 3));
		}
		for (int i = 0; i < after.length() - 3; i += 3) {
			logProbability += getLogProbabilityGeneral(after.substring(i, i + 3));
		}
		return logProbability;
	}

	/**
	 * Stop-model probability: what is the log of the probability of observing the
	 * given sequences before and after a TGA, assuming that the TGA is a stop
	 * 
	 * @param before
	 * @param after
	 * @return log probability
	 */
	public static double stop(String before, String after) {
		double logProbability = 0;
		/*
		 * for(char c : before.toCharArray()) { logProbability +=
		 * getLogProbabilityBaseCDS(c); } for(char c : after.toCharArray()) {
		 * logProbability += getLogProbabilityBaseUTR(c); }
		 */
		for (int i = 0; i < before.length() - 3; i += 3) {
			logProbability += getLogProbabilityRestricted(before.substring(i, i + 3));
		}
		for (char c : after.toCharArray()) {
			logProbability += getLogProbabilityBaseUTR(c);
		}

		return logProbability;
	}

	/**
	 * @param occurrences this list will be filled with occurrence-indices: They can
	 *                    be both negative and positive, with negative indices
	 *                    simply denoting that the codon in question is TCA, not
	 *                    TGA. The magnitude of an occurrence o denotes the leftmost
	 *                    index where the codon (TGA or TCA) is found in the
	 *                    <b>forward</b> strand
	 * @param sequence    the sequence wherein to find stop occurrences (occurrences
	 *                    of TGA as stop-TGA)
	 * @param upstream    size of the upstream window
	 * @param downstream  size of the downstream window
	 * @param forward     whether forward strand: If reverse strand, the occurrences
	 *                    will be marked by a negative sign.
	 */
	public static void findStopOccurrences(List<Integer> occurrences, String sequence, int upstream, int downstream,
			boolean forward) {
		// TODO: this should better be: upstream is in codons, and currentUpstream is in
		// nt (no possibility of unthreeven windows)
		int currentUpstream = upstream;
		int currentDownstream = downstream;
		for (int i = sequence.indexOf("TGA"); i >= 0; i = sequence.indexOf("TGA", i + 3)) {
			// https://stackoverflow.com/questions/5034442/indexes-of-all-occurrences-of-character-in-a-string
			if (i < upstream) {
				currentUpstream = Math.floorDiv(i, 3) * 3;
			}
			if (i + 3 >= sequence.length() - downstream) {
				currentDownstream = Math.floorDiv(sequence.length() - (i + 3), 3) * 3;
			}

			String before = sequence.substring(i - currentUpstream, i);
			String after = sequence.substring(i + 3, i + 3 + currentDownstream);

			double logProbabilityStop = stop(before, after);
			double logProbabilitySense = sense(before, after);
			double logProbabilityNCS = noncoding(before, after);

			if (logProbabilityStop > logProbabilityNCS && logProbabilityStop > logProbabilitySense) {
				occurrences.add(forward ? i : (i - sequence.length() + 3));
			}

			currentUpstream = upstream;
			currentDownstream = downstream;
		}
	}

	/**
	 * Classifies all TGAs in the given sequence (both on forward and reverse
	 * strand, i.e. TGA and TCA) as either stop-TGA, sense-TGA or TGA in an
	 * untranslated region. Replaces all but stop-TGAs with TGG/CCA respectively
	 * (classification is highly imperfect)
	 * 
	 * @param sequence TODO: are there sideeffects?
	 * @param utrs:    a list of annotation-elements for comparison; these will be
	 *                 used to print information about whether correct stops were
	 *                 indeed found to std.out
	 * @return the modified sequence
	 */
	public static String replaceNonStopUGA(String sequence, List<ThreeUTR> utrs) {
		// step 1: replace stop-UGA/ACU with special code: OPQ/QPO (indicating frame by
		// direction)
		ArrayList<Integer> occurrences = new ArrayList<Integer>();

		// 15 seems to be underpowered
		int upstream = 24; // 18;
		int downstream = 18;

		// step 1.a)
		findStopOccurrences(occurrences, sequence, upstream, downstream, true);
		findStopOccurrences(occurrences, reverseComplement(sequence), upstream, downstream, false);

		boolean found = false;

		// Step 1.b: replace all stop-UGA-occurrences by OPQ
		char[] sequenceChars = sequence.toCharArray();
		for (int o : occurrences) {
			// TODO: check against annotation: are these enough and the right ones?
			// -> how many of the stop-TGA/TCA have been found?
			if (utrs != null) {

				for (ThreeUTR utr : utrs) {
					if (o > 0 && utr.isForward() && utr.getStopPosition() == o) {
						found = true;
						System.out.println("Found correct stop at " + o + " (forward)");
					} else if (o < 0 && !utr.isForward() && utr.getStopPosition() == -o) {
						found = true;
						System.out.println("Found correct stop at " + o + " (reverse)");
					}
				}

			}

			if (o > 0) {
				sequenceChars[o] = 'O';
				sequenceChars[o + 1] = '>';
				sequenceChars[o + 2] = 'Q';
			} else {
				sequenceChars[-o] = 'O';
				sequenceChars[1 - o] = '<';
				sequenceChars[2 - o] = 'Q';
			}
		}
		if (utrs != null && !found) {
			System.out.println("! Did not find a single correct stop\n");
		} else if (utrs != null) {
			System.out.println("+ + + + Found > 0 correct stops\n+ + + +\n");
		}

		sequence = String.valueOf(sequenceChars);

		// step 2: replace all remaining TGA->TGG and rev.comp. ACT->CCT
		sequence = sequence.replaceAll("TGA", "TGG");
		sequence = sequence.replaceAll("TCA", "CCA");

		// step 3: unreplace UGA/ACU <- OPQ/QPO
		sequence = sequence.replaceAll("O", "T");
		sequence = sequence.replaceAll(">", "G");
		sequence = sequence.replaceAll("<", "C");
		sequence = sequence.replaceAll("Q", "A");
		return sequence;
	}

	/**
	 * @param args:
	 * <ol>
	 *                  <li>Input scaffolds.fasta</li>
	 *                  <li>Output file.fasta</li>
	 *                  <li>Optional: annotation-file.gff</li>
	 *                  </ol>
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			System.out.println("Give three arguments:\n\t1: input scaffolds.fasta"
					+ "\n\t2: output file.fasta (will be overridden)\n\t3: (optional) annotation-file.gff");
		}

		for (int i = 0; i < args.length; i++)
			System.out.println("args[" + i + "]=" + args[i]);
		HashMap<String, List<ThreeUTR>> annotation = new HashMap<String, List<ThreeUTR>>();

		if (args.length == 3) {
			BufferedReader reader = new BufferedReader(new FileReader(new File(args[2])));
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.startsWith("#"))
					continue;
				if (line.startsWith(">"))
					break;

				String[] lineContents = line.split("\t");
				if (lineContents[2].equals("three_prime_UTR")) {
					if (!annotation.containsKey(lineContents[0])) {
						annotation.put(lineContents[0], new ArrayList<ThreeUTR>());
					}
					annotation.get(lineContents[0]).add(new ThreeUTR(Integer.parseInt(lineContents[3]) - 1,
							Integer.parseInt(lineContents[4]), lineContents[6].equals("+")));
				}
			}
			reader.close();
		}

		if (args.length > 0) {
			File infile = new File(args[0]);
			System.out.println("Reading from file: " + infile);

			String currentHeader = null;
			String currentSequence = null;

			BufferedWriter writer = new BufferedWriter(new FileWriter(args[1]));
			BufferedReader reader = new BufferedReader(new FileReader(infile));
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				if (line.startsWith(">")) {
					if (currentHeader != null) {
						System.out.println("Cleaning " + currentHeader + ", for which the annotation is: "
								+ annotation.get(currentHeader.substring(1)));

						String cleaned = replaceNonStopUGA(currentSequence, annotation.get(currentHeader.substring(1)));
						writer.write(currentHeader);
						writer.newLine();
						writer.write(cleaned);
						writer.newLine();
					}
					currentHeader = line;
					currentSequence = "";
				} else {
					currentSequence += line;
				}
			}

			System.out.println("Cleaning " + currentHeader + ", for which the annotation is: "
					+ annotation.get(currentHeader.substring(1)));

			String cleaned = replaceNonStopUGA(currentSequence, annotation.get(currentHeader.substring(1)));
			writer.write(currentHeader);
			writer.newLine();
			writer.write(cleaned);
			writer.newLine();
			writer.close();
			reader.close();
		}

		System.out.println("\nDone");
	}
}
