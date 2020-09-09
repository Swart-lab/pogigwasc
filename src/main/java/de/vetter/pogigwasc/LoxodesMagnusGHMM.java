package de.vetter.pogigwasc;

import de.vetter.pogigwasc.states.CodingState;
import de.vetter.pogigwasc.states.HMMState;
import de.vetter.pogigwasc.states.InterruptedCodonState;
import de.vetter.pogigwasc.states.IntronState;
import de.vetter.pogigwasc.states.NoncodingState;
import de.vetter.pogigwasc.states.StartRegionState;
import de.vetter.pogigwasc.states.StopRegionState;

/**
 * Class for instantiating the particular GHMM used for predicting genes in
 * <i>Loxodes magnus</i> (probably fit for other karyorelicts too). <br>
 * Takes care of setting up the states and transitions, and provides a method
 * {@link #parseToGFF(String, Parse, ModelParameters)} to convert a parse of
 * this GHMM into gff-format
 * 
 * @author David Emanuel Vetter
 *
 */
public class LoxodesMagnusGHMM extends GHMM {
	/** State-names: Name states +... for forward strand, and -... for reverse strand */
	public static final String NCS = "NCS";
	public static final String FORWARD_START = "+M";
	public static final String FORWARD_CDS = "+CDS";
	public static final String FORWARD_STOP = "+Stop";
	public static final String FORWARD_INTRON_0_0 = "+intron 0-0";
	public static final String FORWARD_INTRON_1_2 = "+intron 1-2";
	public static final String FORWARD_INTRON_2_1 = "+intron 2-1";
	public static final String FORWARD_PREINTRON_ONE = "+1 nt before intron";
	public static final String FORWARD_PREINTRON_TWO = "+2 nts before intron";
	public static final String FORWARD_POSTINTRON_ONE = "+1 nt after intron";
	public static final String FORWARD_POSTINTRON_TWO = "+2 nts after intron";

	public static final String REVERSE_START = "-M";
	public static final String REVERSE_CDS = "-CDS";
	public static final String REVERSE_STOP = "-Stop";
	public static final String REVERSE_INTRON_0_0 = "-intron 0-0";
	public static final String REVERSE_INTRON_1_2 = "-intron 1-2";
	public static final String REVERSE_INTRON_2_1 = "-intron 2-1";
	public static final String REVERSE_PREINTRON_ONE = "-1 nt before intron";
	public static final String REVERSE_PREINTRON_TWO = "-2 nts before intron";
	public static final String REVERSE_POSTINTRON_ONE = "-1 nt after intron";
	public static final String REVERSE_POSTINTRON_TWO = "-2 nts after intron";
	
	public LoxodesMagnusGHMM(final ModelParameters modelParameters) {
		super();
		/** 
		 * Setting up the model 
		 */
		addState(new NoncodingState(NCS, modelParameters));

		/** FORWARD STRAND */
		addState(new StartRegionState(FORWARD_START, true, modelParameters)); 
		addState(new CodingState(FORWARD_CDS, true, modelParameters));
		addState(new StopRegionState(FORWARD_STOP, true, modelParameters));

		/** simplest intron */
		addState(new IntronState(FORWARD_INTRON_0_0, true, modelParameters));

		/** intron preceded by one nt and followed by two nt of a codon */
		addState(new InterruptedCodonState(FORWARD_PREINTRON_ONE, true, true, 1, modelParameters));
		addState(new IntronState(FORWARD_INTRON_1_2, true, modelParameters));
		addState(new InterruptedCodonState(FORWARD_POSTINTRON_TWO, true, false, 2, modelParameters));

		/** intron preceded by two nt and followed by one nt of a codon */
		addState(new InterruptedCodonState(FORWARD_PREINTRON_TWO, true, true, 2, modelParameters));
		addState(new IntronState(FORWARD_INTRON_2_1, true, modelParameters));
		addState(new InterruptedCodonState(FORWARD_POSTINTRON_ONE, true, false, 1, modelParameters));

		/** REVERSE STRAND */
		addState(new StartRegionState(REVERSE_START, false, modelParameters));
		addState(new CodingState(REVERSE_CDS, false, modelParameters));
		addState(new StopRegionState(REVERSE_STOP, false, modelParameters));

		/** simplest intron */
		addState(new IntronState(REVERSE_INTRON_0_0, false, modelParameters));

		/** intron preceded by one nt and followed by two nt of a codon */
		addState(new InterruptedCodonState(REVERSE_PREINTRON_ONE, false, true, 1, modelParameters));
		addState(new IntronState(REVERSE_INTRON_1_2, false, modelParameters));
		addState(new InterruptedCodonState(REVERSE_POSTINTRON_TWO, false, false, 2, modelParameters));

		/** intron preceded by two nt and followed by one nt of a codon */
		addState(new InterruptedCodonState(REVERSE_PREINTRON_TWO, false, true, 2, modelParameters));
		addState(new IntronState(REVERSE_INTRON_2_1, false, modelParameters));
		addState(new InterruptedCodonState(REVERSE_POSTINTRON_ONE, false, false, 1, modelParameters));

		initialiseTransitionMatrix();
		clearTransitionMatrix();

		/** Setting the transitions */
		setTransitionProbability(1, 1, 1d); // stay in terminal
		setTransitionProbability(0, 2, 1d);

		double probabilityStayNCS = modelParameters.getProbabilityOfStayingInNCS();

		setTransitionProbability(2, 2, probabilityStayNCS);
		setTransitionProbability(2, 3, 0.4 * (1 - probabilityStayNCS)); // NCS -> +M
		setTransitionProbability(2, 15, 0.4 * (1 - probabilityStayNCS)); // NCS -> -Stop
		setTransitionProbability(2, 1, 0.2 * (1 - probabilityStayNCS)); // NCS -> terminal

		setTransitionProbability(3, 4, 1d); // +M -> +CDS

		setTransitionProbability(4, 4, modelParameters.getProbabilityOfStayingInCDS()); // +CDS -> +CDS
		setTransitionProbability(4, 5, modelParameters.getProbabilityGeneEnds()); // +CDS -> +Stop
		setTransitionProbability(4, 6, modelParameters.getProbabilityCDSToIntron()); // +CDS -> +intron 0-0
		setTransitionProbability(4, 7, modelParameters.getProbabilityCDSToIntron()); // +CDS -> + 1 nt before
		setTransitionProbability(4, 10, modelParameters.getProbabilityCDSToIntron()); // +CDS -> + 2 nts before

		setTransitionProbability(5, 2, 1d); // +stop -> NCS always
		setTransitionProbability(6, 4, 1d); // + intron 0-0 -> +CDS always

		setTransitionProbability(7, 8, 1d);
		setTransitionProbability(8, 9, 1d);
		setTransitionProbability(9, 4, 1d);

		setTransitionProbability(10, 11, 1d);
		setTransitionProbability(11, 12, 1d);
		setTransitionProbability(12, 4, 1d);

		// Reverse model:
		setTransitionProbability(15, 14, 1d); // -Stop -> -CDS

		setTransitionProbability(14, 14, modelParameters.getProbabilityOfStayingInCDS()); // stay in -CDS
		setTransitionProbability(14, 13, modelParameters.getProbabilityGeneEnds()); // -CDS -> -M
		setTransitionProbability(14, 16, modelParameters.getProbabilityCDSToIntron()); // -CDS -> -intron 0-0
		setTransitionProbability(14, 17, modelParameters.getProbabilityCDSToIntron()); // -CDS -> - 1 nt before (=
																							// to the left of) intron
		setTransitionProbability(14, 20, modelParameters.getProbabilityCDSToIntron()); // -CDS -> - 2 nts before (=
																							// to the left of) intron

		setTransitionProbability(13, 2, 1d); // -M -> NCS always

		setTransitionProbability(16, 14, 1d); // -intron 0-0 -> -CDS

		setTransitionProbability(17, 18, 1d); // - 1 nt before -> -intron 1-2
		setTransitionProbability(18, 19, 1d); // -intron 1-2 -> -2 nts after
		setTransitionProbability(19, 14, 1d); // -2 nts after -> -CDS

		setTransitionProbability(20, 21, 1d);
		setTransitionProbability(21, 22, 1d);
		setTransitionProbability(22, 14, 1d);
	}
	
	/**
	 * Translates the given parse into GFF-format, using the names and the semantic defined in this class.
	 * @param header the contig-name on which prediction was performed
	 * @param parse a parse on that contig
	 * @return The GFF-representation of that parse, ending in a newline
	 */
	public static String parseToGFF(String header, Parse parse, ModelParameters parameters) {
		StringBuffer result = new StringBuffer();
		int currentPos = 1;

		GFFFeature currentFeature = null;
		int startOfCurrentFeature = -1;

		for(int i = 0; i < parse.getNumberOfSteps(); i++) {
			HMMState state = parse.get(i).getFirst();
			switch (state.getName()) {
			case LoxodesMagnusGHMM.NCS:
				startOfCurrentFeature = currentPos + 1; // for reverse-stop
				currentFeature = null;
				break;
			case LoxodesMagnusGHMM.FORWARD_START:
			case LoxodesMagnusGHMM.FORWARD_POSTINTRON_ONE:
			case LoxodesMagnusGHMM.FORWARD_POSTINTRON_TWO:
			case LoxodesMagnusGHMM.FORWARD_CDS:
			case LoxodesMagnusGHMM.FORWARD_PREINTRON_ONE:
			case LoxodesMagnusGHMM.FORWARD_PREINTRON_TWO:

			case LoxodesMagnusGHMM.REVERSE_POSTINTRON_ONE:
			case LoxodesMagnusGHMM.REVERSE_POSTINTRON_TWO:
			case LoxodesMagnusGHMM.REVERSE_CDS:
			case LoxodesMagnusGHMM.REVERSE_PREINTRON_ONE:
			case LoxodesMagnusGHMM.REVERSE_PREINTRON_TWO:
				if (currentFeature != null && currentFeature != GFFFeature.CDS) {
					result.append(header + "\tpredicted\t" + currentFeature.getCode() + "\t"
							+ startOfCurrentFeature + "\t" + (currentPos - 1) + "\t.\t"
							+ (state.getName().startsWith("+") ? "+" : "-") + "\t.\t ");
					result.append(System.lineSeparator());

				}

				if (currentFeature != GFFFeature.CDS)
					startOfCurrentFeature = state.getName().equals(LoxodesMagnusGHMM.FORWARD_START)
							? currentPos + (parameters.getStartRegionSize() - 3)
							: currentPos;
				// TODO ^ notice this! It is for the new StartRegion-state!
				currentFeature = GFFFeature.CDS;

				break;
			case LoxodesMagnusGHMM.FORWARD_INTRON_0_0:
			case LoxodesMagnusGHMM.FORWARD_INTRON_1_2:
			case LoxodesMagnusGHMM.FORWARD_INTRON_2_1:
			case LoxodesMagnusGHMM.REVERSE_INTRON_0_0:
			case LoxodesMagnusGHMM.REVERSE_INTRON_1_2:
			case LoxodesMagnusGHMM.REVERSE_INTRON_2_1:
				// some kind of CDS must precede
				result.append(header + "\tpredicted\t" + currentFeature.getCode() + "\t"
						+ startOfCurrentFeature + "\t" + (currentPos - 1) + "\t.\t"
						+ (state.getName().startsWith("+") ? "+" : "-") + "\t.\t ");
				result.append(System.lineSeparator());
				currentFeature = GFFFeature.INTRON;
				startOfCurrentFeature = currentPos;
				break;

			case LoxodesMagnusGHMM.FORWARD_STOP:
				int firstbaseOfStop = currentPos + parameters.getStopRegionSize() - 3;
				// CDS must precede
				result.append(header + "\tpredicted\t" + currentFeature.getCode() + "\t"
						+ startOfCurrentFeature + "\t" + (firstbaseOfStop - 1) + "\t.\t+\t.\t ");
				result.append(System.lineSeparator());
				// The actual stop
				result.append(header + "\tpredicted\tstop_codon\t" + firstbaseOfStop + "\t"
						+ (firstbaseOfStop + 2) + "\t.\t+\t.\t ");
				result.append(System.lineSeparator());
				currentFeature = null; // NCS follows
				startOfCurrentFeature = currentPos + parameters.getStopRegionSize(); 
				break;
			case LoxodesMagnusGHMM.REVERSE_START:
				// Coming from CDS: -M is to be counted into the CDS
				result.append(header + "\tpredicted\t" + currentFeature.getCode() + "\t"
						+ startOfCurrentFeature + "\t" + (currentPos + 2) + "\t.\t-\t.\t ");
				result.append(System.lineSeparator());
				currentFeature = null; // NCS follows
				startOfCurrentFeature = currentPos + (parameters.getStartRegionSize() - 3); // that's where the NCS starts
				break;
			case LoxodesMagnusGHMM.REVERSE_STOP:
				result.append(header + "\tpredicted\tstop_codon\t" + currentPos + "\t" + (currentPos + 2)
						+ "\t.\t-\t.\t ");
				result.append(System.lineSeparator());
				currentFeature = GFFFeature.CDS;
				startOfCurrentFeature = currentPos + 3;
				break;
			}

			currentPos += parse.get(i).getSecond();
		}
		return result.toString();
	}

}
