package de.vetter.masterthesis;

import de.vetter.masterthesis.states.*;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
    	double a = Double.NEGATIVE_INFINITY;
        System.out.println( "Hello World!" );
        System.out.println("2^" + a + "=" + Math.pow(2, a));
        
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
        ghmm.addState(new HMMState("+single nt before intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// TODO Auto-generated method stub
				// This is the first base of a codon -> use empirical 1st-position probability
				return 0;
			}
        });
        ghmm.addState(new IntronState("+intron 1-2", true));
        ghmm.addState(new HMMState("+two nts after intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// TODO Auto-generated method stub
				// These are the second and third of a codon -> use respective empirical distrs
				return 0;
			}
        });
        
        /** intron preceded by two nt and followed by one nt of a codon */
        ghmm.addState(new HMMState("+two nts before intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// TODO Auto-generated method stub
				// These are the first and second of a codon -> use respective empirical distrs
				return 0;
			}
        });
        ghmm.addState(new IntronState("+intron 2-1", true));
        ghmm.addState(new HMMState("+single nt after intron") {
			@Override
			public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
				// TODO Auto-generated method stub:
				// This is the final base of a codon -> use empirical 3rd-position probability
				return 0;
			}
        });
        
        /** REVERSE STRAND */
        // TODO
    }
}
