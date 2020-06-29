package de.vetter.masterthesis;

import java.util.List;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import de.vetter.masterthesis.states.*;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args ) throws IOException
    {
    	if(args.length != 2) {
    		System.out.println("usage: > ___.jar [input sequences].fasta [output-file].gff");
    		return;
    	}
    	
    	BufferedReader reader = new BufferedReader(new FileReader(new File(args[0])));
    	BufferedWriter writer = new BufferedWriter(new FileWriter(args[1]));
    	
    	// double a = Double.NEGATIVE_INFINITY;
        // System.out.println( "Hello World!" );
        // System.out.println("2^" + a + "=" + Math.pow(2, a));
        
        // TODO: preprocess: replace N? Definitely toUpperCase
        
        GHMM ghmm = new GHMM();
        ghmm.addState(new NoncodingState("NCS"));
        
        /** FORWARD STRAND */
        ghmm.addState(new FixedSequenceState("+M", "ATG")); // Name states + for forward, and - for backward strand
        ghmm.addState(new CodingState("+CDS", true));
        ghmm.addState(new StopRegionState("+Stop", true));
        
        /** simplest intron */
        ghmm.addState(new IntronState("+intron 0-0", true));
        
        /** intron preceded by one nt and followed by two nt of a codon */
        ghmm.addState(new HMMState("+1 nt before intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the first base of a codon -> use empirical 1st-position probability
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return Utilities.getLogBaseProbabilityCDS(newEmission.charAt(0), 0);
			}
        });
        ghmm.addState(new IntronState("+intron 1-2", true));
        ghmm.addState(new HMMState("+2 nts after intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the second and third of a codon -> use respective empirical distrs
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				return Utilities.getLogBaseProbabilityCDS(newEmission.charAt(0), 1)
						+ Utilities.getLogBaseProbabilityCDS(newEmission.charAt(1), 2);
			}
        });
        
        /** intron preceded by two nt and followed by one nt of a codon */
        ghmm.addState(new HMMState("+2 nts before intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the first and second of a codon -> use respective empirical distrs
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				return Utilities.getLogBaseProbabilityCDS(newEmission.charAt(0), 0)
						+ Utilities.getLogBaseProbabilityCDS(newEmission.charAt(1), 1);
			}
        });
        ghmm.addState(new IntronState("+intron 2-1", true));
        ghmm.addState(new HMMState("+1 nt after intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the final base of a codon -> use empirical 3rd-position probability
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return Utilities.getLogBaseProbabilityCDS(newEmission.charAt(0), 2);
			}
        });
        
        
        /** REVERSE STRAND */
        ghmm.addState(new FixedSequenceState("-M", "CAT")); // Name states + for forward, and - for backward strand
        ghmm.addState(new CodingState("-CDS", false));
        ghmm.addState(new StopRegionState("-Stop", false));
        
        /** simplest intron */
        ghmm.addState(new IntronState("-intron 0-0", false));
        
        /** intron preceded by one nt and followed by two nt of a codon */
        ghmm.addState(new HMMState("-1 nt before intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the last base of a codon
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return Utilities.getLogBaseProbabilityCDS(Utilities.reverseComplement(newEmission).charAt(0), 2);
			}
        });
        ghmm.addState(new IntronState("-intron 1-2", false));
        ghmm.addState(new HMMState("-2 nts after intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the first and second of a codon
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				newEmission = Utilities.reverseComplement(newEmission);
				return Utilities.getLogBaseProbabilityCDS(newEmission.charAt(0), 0)
						+ Utilities.getLogBaseProbabilityCDS(newEmission.charAt(1), 1);
			}
        });
        
        /** intron preceded by two nt and followed by one nt of a codon */
        ghmm.addState(new HMMState("-2 nts before intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// These are the third and second of a codon
				if(newEmission.length() != 2)
					return Double.NEGATIVE_INFINITY;
				
				newEmission = Utilities.reverseComplement(newEmission);
				return Utilities.getLogBaseProbabilityCDS(newEmission.charAt(0), 1)
						+ Utilities.getLogBaseProbabilityCDS(newEmission.charAt(1), 2);
			}
        });
        ghmm.addState(new IntronState("-intron 2-1", true));
        ghmm.addState(new HMMState("-1 nt after intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// This is the first base of a codon
				if(newEmission.length() != 1)
					return Double.NEGATIVE_INFINITY;
				
				return Utilities.getLogBaseProbabilityCDS(Utilities.reverseComplement(newEmission).charAt(0), 0);
			}
        });        
        
        
        
        ghmm.initialiseTransitionMatrix();
        ghmm.clearTransitionMatrix();
        
        
        /** Setting the transitions */
        ghmm.setTransitionProbability(1, 1, 1d); // stay in terminal
        ghmm.setTransitionProbability(0, 2, 1d);
        
        ghmm.setTransitionProbability(2, 2, 397d / 400d);
        ghmm.setTransitionProbability(2, 3, 1d / 400d); // NCS -> +M
        ghmm.setTransitionProbability(2, 15, 1d / 400d); // NCS -> -Stop
        ghmm.setTransitionProbability(2, 1, 1d / 400d); // NCS -> terminal
        
        ghmm.setTransitionProbability(3, 4, 1d); // +M -> +CDS
        
        ghmm.setTransitionProbability(4, 4, 1d - (1d/1200d)); // stay in +CDS (exonlength median=960, mean=1250 -> approx 1200)
        ghmm.setTransitionProbability(4, 5, (1d/1200d) * 0.61d); // +CDS -> +Stop: from introns / gene: geometric with 0 -> empirical p ~ 61%
        ghmm.setTransitionProbability(4, 6, (1d/1200d) * 0.13d); // +CDS -> +intron 0-0
        ghmm.setTransitionProbability(4, 7, (1d/1200d) * 0.13d); // +CDS -> + 1 nt before
        ghmm.setTransitionProbability(4, 10, (1d/1200d) * 0.13d); // +CDS -> + 2 nts before
        
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
        
        ghmm.setTransitionProbability(14, 14, 1d - (1d/1200d)); // stay in -CDS (see above)
        ghmm.setTransitionProbability(14, 13, (1d/1200d) * 0.61d); // -CDS -> -M: from introns / gene: geometric with 0 -> empirical p ~ 61%
        ghmm.setTransitionProbability(14, 16, (1d/1200d) * 0.13d); // -CDS -> -intron 0-0
        ghmm.setTransitionProbability(14, 17, (1d/1200d) * 0.13d); // -CDS -> - 1 nt before (= to the left of) intron
        ghmm.setTransitionProbability(14, 20, (1d/1200d) * 0.13d); // -CDS -> - 2 nts before (= to the left of) intron
        
        ghmm.setTransitionProbability(13, 2, 1d); // -M -> NCS always
        
        ghmm.setTransitionProbability(16, 14, 1d); // -intron 0-0 -> -CDS
        
        ghmm.setTransitionProbability(17, 18, 1d); // - 1 nt before -> -intron 1-2
        ghmm.setTransitionProbability(18, 19, 1d); // -intron 1-2 -> -2 nts after
        ghmm.setTransitionProbability(19, 14, 1d); // -2 nts after -> -CDS
        
        ghmm.setTransitionProbability(20, 21, 1d);
        ghmm.setTransitionProbability(21, 22, 1d);
        ghmm.setTransitionProbability(22, 14, 1d);

        
        System.out.println(ghmm);
        
        
        /** Now read in the fasta file sequence for sequence, and for each sequence generate a gff-output */
        String currentHeader = null;
		String currentSequence = null;

		for (String line = reader.readLine(); line != null; line = reader.readLine()) {
			if (line.startsWith(">")) {
				if (currentHeader != null) {
					doPredictions(ghmm, writer, currentHeader, currentSequence);
				}
				currentHeader = line;
				currentSequence = "";
			} else {
				currentSequence += line;
			}
		}
		
		// on last sequence
		if (currentHeader != null) {
			doPredictions(ghmm, writer, currentHeader, currentSequence.substring(9000, 12000));
		}
    }
    
    public static void doPredictions(GHMM ghmm, BufferedWriter writer, String currentHeader, String currentSequence) {
    	System.out.println("\nParses for " + currentHeader + ":\n");
		Viterbi viterbi = new Viterbi(ghmm, currentSequence);
		for(List<Pair<HMMState, Integer>> parse : viterbi.computeParses()) {
			System.out.print("\nParse:");
			int currentPos = 1;
			boolean openCDS = false;
			for(Pair<HMMState, Integer> s : parse) {
				currentPos += s.getSecond(); 
				if(s.getFirst().getName() == "NCS") {
					continue;
				} else if(s.getFirst() instanceof CodingState) {
					if(!openCDS)
						System.out.println("  " + s.getFirst().getName() + " [" + (currentPos - s.getSecond()) + ", ");
					openCDS = true;
				} else {
					if(openCDS)
						System.out.println(currentPos-1 + "]");
					openCDS = false;
					System.out.print("  " + s.getFirst().getName() + " [" + (currentPos - s.getSecond() + 1) + ", " + currentPos + "]");
				}
			}
		}
    }
}
