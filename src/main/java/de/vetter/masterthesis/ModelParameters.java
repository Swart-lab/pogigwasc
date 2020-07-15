package de.vetter.masterthesis;

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
	private double conditionalProbabilityOfChangingToStopWhenLeavingCDS; // <- geometric estimate from intron-count-distr
	
	// TODO: collect all the relevant parameters here
	// give all of them a getter
	// constructor which reads the parameter-file (which should be possibly given as a programm-argument --parameters=...
	// This then also requires a File-format -- find one/come up with one and document!
	
}
