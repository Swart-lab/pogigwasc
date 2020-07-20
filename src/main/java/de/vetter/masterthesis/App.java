package de.vetter.masterthesis;

import java.util.List;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import de.vetter.masterthesis.states.*;

/**
 * Hello world!
 *
 */
public class App 
{
	private static final String NCS = "NCS";
	private static final String FORWARD_START = "+M";
	private static final String FORWARD_CDS = "+CDS";
	private static final String FORWARD_STOP = "+Stop";
	private static final String FORWARD_INTRON_0_0 = "+intron 0-0";
	private static final String FORWARD_INTRON_1_2 = "+intron 1-2";
	private static final String FORWARD_INTRON_2_1 = "+intron 2-1";
	private static final String FORWARD_PREINTRON_ONE = "+1 nt before intron";
	private static final String FORWARD_PREINTRON_TWO = "+2 nts before intron";
	private static final String FORWARD_POSTINTRON_ONE = "+1 nt after intron";
	private static final String FORWARD_POSTINTRON_TWO = "+2 nts after intron";
	
	private static final String REVERSE_START = "-M";
	private static final String REVERSE_CDS = "-CDS";
	private static final String REVERSE_STOP = "-Stop";
	private static final String REVERSE_INTRON_0_0 = "-intron 0-0";
	private static final String REVERSE_INTRON_1_2 = "-intron 1-2";
	private static final String REVERSE_INTRON_2_1 = "-intron 2-1";
	private static final String REVERSE_PREINTRON_ONE = "-1 nt before intron";
	private static final String REVERSE_PREINTRON_TWO = "-2 nts before intron";
	private static final String REVERSE_POSTINTRON_ONE = "-1 nt after intron";
	private static final String REVERSE_POSTINTRON_TWO = "-2 nts after intron";
	
    public static void main( String[] args ) throws IOException
    {
    	LocalDateTime now = LocalDateTime.now();
    	if(args.length != 1 && args.length != 2) {
    		System.out.println("usage: > ___.jar [input sequences].fasta (optional 2nd arg: [output-file].gff)");
    		return;
    	}
    	
    	/** File management */
    	
    	File input = new File(args[0]);
    	BufferedReader reader = new BufferedReader(new FileReader(input));
    	
    	File output;
    	BufferedWriter writer;
    	if(args.length == 1) {
    		// cut off the .fasta (or something like .fa)
    		String genericOutputFilename = args[0].substring(0, args[0].lastIndexOf("."));
    		genericOutputFilename += "-day" + DateTimeFormatter.ofPattern("yyyy-MM-dd_HHmmss").format(now);
    		genericOutputFilename += ".predicted.gff";
    		
    		output = new File(genericOutputFilename);
    	} else {
    		output = new File(args[1]);
    	}
    	
    	writer = new BufferedWriter(new FileWriter(output));
    	
    	// TODO: allow parameter-file to be supplied as cmd-line option
    	File parameterFile = new File("resources\\de\\vetter\\masterthesis\\parameter\\parameters-examplefile.properties");
		final ModelParameters modelParameters = new ModelParameters(new FileReader(parameterFile));
        
        
		/** Setting up the model */
		
        GHMM ghmm = new GHMM();
        ghmm.addState(new NoncodingState(NCS, modelParameters));
        
        /** FORWARD STRAND */
        ghmm.addState(new StartRegionState(FORWARD_START, true, modelParameters));
        // ghmm.addState(new FixedSequenceState(FORWARD_START, "ATG")); // Name states + for forward, and - for backward strand
        ghmm.addState(new CodingState(FORWARD_CDS, true, modelParameters));
        ghmm.addState(new StopRegionState(FORWARD_STOP, true, modelParameters));
        
        /** simplest intron */
        ghmm.addState(new IntronState(FORWARD_INTRON_0_0, true, modelParameters));
        
        /** intron preceded by one nt and followed by two nt of a codon */
        ghmm.addState(new ConstantLengthState(FORWARD_PREINTRON_ONE, 1) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the first base of a codon -> use empirical 1st-position probability
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(0), 0);
			}
        });
        ghmm.addState(new IntronState(FORWARD_INTRON_1_2, true, modelParameters));
        ghmm.addState(new ConstantLengthState(FORWARD_POSTINTRON_TWO, 2) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the second and third of a codon -> use respective empirical distrs
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				return modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(0), 1)
						+ modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(1), 2);
			}
        });
        
        /** intron preceded by two nt and followed by one nt of a codon */
        ghmm.addState(new ConstantLengthState(FORWARD_PREINTRON_TWO, 2) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the first and second of a codon -> use respective empirical distrs
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				return modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(0), 0)
						+ modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(1), 1);
			}
        });
        ghmm.addState(new IntronState(FORWARD_INTRON_2_1, true, modelParameters));
        ghmm.addState(new ConstantLengthState(FORWARD_POSTINTRON_ONE, 1) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the final base of a codon -> use empirical 3rd-position probability
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(0), 2);
			}
        });
        
        
        /** REVERSE STRAND */
        ghmm.addState(new StartRegionState(REVERSE_START, false, modelParameters));
        // ghmm.addState(new FixedSequenceState(REVERSE_START, "CAT")); // Name states + for forward, and - for backward strand
        ghmm.addState(new CodingState(REVERSE_CDS, false, modelParameters));
        ghmm.addState(new StopRegionState(REVERSE_STOP, false, modelParameters));
        
        /** simplest intron */
        ghmm.addState(new IntronState(REVERSE_INTRON_0_0, false, modelParameters));
        
        /** intron preceded by one nt and followed by two nt of a codon */
        ghmm.addState(new ConstantLengthState(REVERSE_PREINTRON_ONE, 1) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the last base of a codon
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return modelParameters.getLogBaseProbabilityCDS(Utilities.reverseComplement(newEmission).charAt(0), 2);
			}
        });
        ghmm.addState(new IntronState(REVERSE_INTRON_1_2, false, modelParameters));
        ghmm.addState(new ConstantLengthState(REVERSE_POSTINTRON_TWO, 2) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the first and second of a codon
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				newEmission = Utilities.reverseComplement(newEmission);
				return modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(0), 0)
						+ modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(1), 1);
			}
        });
        
        /** intron preceded by two nt and followed by one nt of a codon */
        ghmm.addState(new ConstantLengthState(REVERSE_PREINTRON_TWO, 2) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the third and second of a codon
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				newEmission = Utilities.reverseComplement(newEmission);
				return modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(0), 1)
						+ modelParameters.getLogBaseProbabilityCDS(newEmission.charAt(1), 2);
			}
        });
        ghmm.addState(new IntronState(REVERSE_INTRON_2_1, false, modelParameters));
        ghmm.addState(new ConstantLengthState(REVERSE_POSTINTRON_ONE, 1) {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the first base of a codon
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return modelParameters.getLogBaseProbabilityCDS(Utilities.reverseComplement(newEmission).charAt(0), 0);
			}
        });        
        
        
        
        ghmm.initialiseTransitionMatrix();
        ghmm.clearTransitionMatrix();
        
        
        /** Setting the transitions */
        ghmm.setTransitionProbability(1, 1, 1d); // stay in terminal
        ghmm.setTransitionProbability(0, 2, 1d);
        
        double probabilityStayNCS = modelParameters.getProbabilityOfStayingInNCS(); //1695d/1700d;
        
        ghmm.setTransitionProbability(2, 2, probabilityStayNCS);
        ghmm.setTransitionProbability(2, 3, 0.4 * (1-probabilityStayNCS)); // NCS -> +M
        ghmm.setTransitionProbability(2, 15, 0.4 * (1-probabilityStayNCS)); // NCS -> -Stop
        ghmm.setTransitionProbability(2, 1, 0.2 * (1-probabilityStayNCS)); // NCS -> terminal
        
        ghmm.setTransitionProbability(3, 4, 1d); // +M -> +CDS
        
        //double probabilityStayCDS     = ;
        // double probabilityCDSEnd      = modelParameters.getProbabilityGeneEnds();
        // double probabilityCDSToIntron = modelParameters.getProbabilityCDSToIntron();
        
        ghmm.setTransitionProbability(4, 4, modelParameters.getProbabilityOfStayingInCDS()); // stay in +CDS (exonlength median=960, mean=1250 -> approx 1200)
        ghmm.setTransitionProbability(4, 5, modelParameters.getProbabilityGeneEnds()); // +CDS -> +Stop: from introns / gene: geometric with 0 -> empirical p ~ 61%
        ghmm.setTransitionProbability(4, 6, modelParameters.getProbabilityCDSToIntron()); // +CDS -> +intron 0-0
        ghmm.setTransitionProbability(4, 7, modelParameters.getProbabilityCDSToIntron()); // +CDS -> + 1 nt before
        ghmm.setTransitionProbability(4,10, modelParameters.getProbabilityCDSToIntron()); // +CDS -> + 2 nts before
        
        ghmm.setTransitionProbability(5, 2, 1d); // +stop -> NCS always
        ghmm.setTransitionProbability(6, 4, 1d); // + intron 0-0 -> +CDS always
        
        ghmm.setTransitionProbability(7, 8, 1d);
        ghmm.setTransitionProbability(8, 9, 1d);
        ghmm.setTransitionProbability(9, 4, 1d);
        
        ghmm.setTransitionProbability(10, 11, 1d);
        ghmm.setTransitionProbability(11, 12, 1d);
        ghmm.setTransitionProbability(12, 4, 1d);
        
        // Reverse model:
        ghmm.setTransitionProbability(15, 14, 1d); // -Stop -> -CDS
        
        ghmm.setTransitionProbability(14, 14, modelParameters.getProbabilityOfStayingInCDS()); // stay in -CDS (see above)
        ghmm.setTransitionProbability(14, 13, modelParameters.getProbabilityGeneEnds()); // -CDS -> -M: from introns / gene: geometric with 0 -> empirical p ~ 61%
        ghmm.setTransitionProbability(14, 16, modelParameters.getProbabilityCDSToIntron()); // -CDS -> -intron 0-0
        ghmm.setTransitionProbability(14, 17, modelParameters.getProbabilityCDSToIntron()); // -CDS -> - 1 nt before (= to the left of) intron
        ghmm.setTransitionProbability(14, 20, modelParameters.getProbabilityCDSToIntron()); // -CDS -> - 2 nts before (= to the left of) intron
        
        ghmm.setTransitionProbability(13, 2, 1d); // -M -> NCS always
        
        ghmm.setTransitionProbability(16, 14, 1d); // -intron 0-0 -> -CDS
        
        ghmm.setTransitionProbability(17, 18, 1d); // - 1 nt before -> -intron 1-2
        ghmm.setTransitionProbability(18, 19, 1d); // -intron 1-2 -> -2 nts after
        ghmm.setTransitionProbability(19, 14, 1d); // -2 nts after -> -CDS
        
        ghmm.setTransitionProbability(20, 21, 1d);
        ghmm.setTransitionProbability(21, 22, 1d);
        ghmm.setTransitionProbability(22, 14, 1d);

        
        System.out.println(ghmm);
        
        /** Write some book-keeping information into the output-file */
        writer.write("##gff-version 3");
        writer.newLine();
        writer.write("##Generated on: " + DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss").format(now));
        writer.newLine();
        writer.write("##Source file: " + input.getAbsolutePath()); writer.newLine();
        writer.write("##Parameter file: " + parameterFile.getAbsolutePath()); writer.newLine();
        
        
        /** Now read in the fasta file sequence for sequence, and for each sequence generate a gff-output */        
        String currentHeader = null;
		String currentSequence = null;

		for (String line = reader.readLine(); line != null; line = reader.readLine()) {
			if (line.startsWith(">")) {
				if (currentHeader != null) {
					doPredictions(ghmm, writer, currentHeader, currentSequence);
				}
				currentHeader = line.substring(1);
				currentSequence = "";
			} else {
				currentSequence += line;
			}
		}
		
		// on last sequence
		if (currentHeader != null) {
			doPredictions(ghmm, writer, currentHeader, currentSequence); // .substring(45700, 48500));
		}
		
		reader.close();
		writer.flush();
		writer.close();
		
		System.out.println("______________________\nWrote to output-file " + output.getAbsolutePath());
    }
    
    /**
     * Performs gene prediction on the given sequence and writes the result as gff into the writer.
     * @param ghmm
     * @param writer
     * @param currentHeader
     * @param currentSequence
     * @throws IOException if something goes wrong while trying to write into the out-file 
     * @throws IllegalArgumentException if the sequence contains N, not just TGAC
     */
    public static void doPredictions(GHMM ghmm, BufferedWriter writer, String currentHeader, String currentSequence) throws IOException {
    	System.out.println("\nPrediction genes in " + currentHeader + ":\n");
    	if(currentSequence.contains("NNN")) {
			throw new IllegalArgumentException(
					"The given sequence contained uncharacterised Nucleotides (N). Please provide contigs, not scaffolds");
    	}
    	
		Viterbi viterbi = new Viterbi(ghmm, currentSequence);
		viterbi.setAbbreviating(true);
		for(List<Pair<HMMState, Integer>> parse : viterbi.computeParses()) {
			System.out.println("\n\tWriting parse to file");
			int currentPos = 1;
			
			GFFFeature currentFeature = null;
			int startOfCurrentFeature = -1;
			
			for(Pair<HMMState, Integer> s : parse) {
				HMMState state = s.getFirst();
				switch (state.getName()) {
				case NCS:
					startOfCurrentFeature = currentPos + 1; // for reverse-stop
					currentFeature = null;
					break;
				case FORWARD_START:
				case FORWARD_POSTINTRON_ONE:
				case FORWARD_POSTINTRON_TWO:
				case FORWARD_CDS:
				case FORWARD_PREINTRON_ONE:
				case FORWARD_PREINTRON_TWO:

				case REVERSE_POSTINTRON_ONE:
				case REVERSE_POSTINTRON_TWO:
				case REVERSE_CDS:
				case REVERSE_PREINTRON_ONE:
				case REVERSE_PREINTRON_TWO:
					if (currentFeature != null && currentFeature != GFFFeature.CDS) {
						writer.write(currentHeader + "\tpredicted\t" + currentFeature.getCode() + "\t"
								+ startOfCurrentFeature + "\t" + (currentPos - 1) + "\t.\t"
								+ (state.getName().startsWith("+") ? "+" : "-") + "\t.\t ");
						writer.newLine();
						
					}

					if (currentFeature != GFFFeature.CDS)
						startOfCurrentFeature = state.getName().equals(FORWARD_START) ? currentPos + 3 : currentPos;
					// TODO ^ notice this! It is for the new StartRegion-state!
					currentFeature = GFFFeature.CDS;

					break;
				case FORWARD_INTRON_0_0:
				case FORWARD_INTRON_1_2:
				case FORWARD_INTRON_2_1:
				case REVERSE_INTRON_0_0:
				case REVERSE_INTRON_1_2:
				case REVERSE_INTRON_2_1:
					// some kind of CDS must precede
					writer.write(currentHeader + "\tpredicted\t" + currentFeature.getCode() + "\t"
							+ startOfCurrentFeature + "\t" + (currentPos - 1) + "\t.\t"
							+ (state.getName().startsWith("+") ? "+" : "-") + "\t.\t ");
					writer.newLine();
					currentFeature = GFFFeature.INTRON;
					startOfCurrentFeature = currentPos;
					break;

				case FORWARD_STOP:
					// CDS must precede
					writer.write(currentHeader + "\tpredicted\t" + currentFeature.getCode() + "\t"
							+ startOfCurrentFeature + "\t" + (currentPos + 20) + "\t.\t+\t.\t ");
					writer.newLine();
					// The actual stop
					writer.write(currentHeader + "\tpredicted\tstop_codon\t" + (currentPos + 21) + "\t"
							+ (currentPos + 23) + "\t.\t+\t.\t ");
					writer.newLine();
					currentFeature = null; // NCS follows
					startOfCurrentFeature = currentPos + 24; // not strictly necessary, following NCS will take care of
																// this.
					break;
				case REVERSE_START:
					// Komme von CDS! -M zählt mit in die CDS.
					writer.write(currentHeader + "\tpredicted\t" + currentFeature.getCode() + "\t"
							+ startOfCurrentFeature + "\t" + (currentPos + 2) + "\t.\t-\t.\t ");
					writer.newLine();
					currentFeature = null; // NCS follows
					startOfCurrentFeature = currentPos + 3; // that's where the NCS starts
					break;
				case REVERSE_STOP:
					writer.write(currentHeader + "\tpredicted\tstop_codon\t" + currentPos + "\t" + (currentPos + 2)
							+ "\t.\t-\t.\t ");
					writer.newLine();
					currentFeature = GFFFeature.CDS;
					startOfCurrentFeature = currentPos + 3;
					break;
				}
				
				currentPos += s.getSecond();
			}
		}
		
		System.out.println("done");
    }
}
