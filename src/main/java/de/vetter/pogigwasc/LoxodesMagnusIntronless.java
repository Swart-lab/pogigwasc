package de.vetter.pogigwasc;

import de.vetter.pogigwasc.states.CodingState;
import de.vetter.pogigwasc.states.NoncodingState;
import de.vetter.pogigwasc.states.StartRegionState;
import de.vetter.pogigwasc.states.StopRegionState;

/**
 * Class for instantiating a GHMM for intron-less gene-prediction (for RNA):
 * Uses the same naming as {@link LoxodesMagnusGHMM}. Thus, parses of this model
 * can be converted to gff by
 * {@link LoxodesMagnusGHMM#parseToGFF(String, Parse, ModelParameters)}
 * 
 * @author David Emanuel Vetter
 *
 */
public class LoxodesMagnusIntronless extends GHMM {
	
	public LoxodesMagnusIntronless(ModelParameters modelParameters) {
		super();
		
		/** 
		 * Setting up the model 
		 */
		addState(new NoncodingState(LoxodesMagnusGHMM.NCS, modelParameters));

		/** FORWARD STRAND */
		addState(new StartRegionState(LoxodesMagnusGHMM.FORWARD_START, true, modelParameters)); 
		addState(new CodingState(LoxodesMagnusGHMM.FORWARD_CDS, true, modelParameters));
		addState(new StopRegionState(LoxodesMagnusGHMM.FORWARD_STOP, true, modelParameters));

		/** REVERSE STRAND */
		addState(new StartRegionState(LoxodesMagnusGHMM.REVERSE_START, false, modelParameters));
		addState(new CodingState(LoxodesMagnusGHMM.REVERSE_CDS, false, modelParameters));
		addState(new StopRegionState(LoxodesMagnusGHMM.REVERSE_STOP, false, modelParameters));
		
		initialiseTransitionMatrix();
		clearTransitionMatrix();
		
		/** Setting the transitions */
		setTransitionProbability(1, 1, 1d); // stay in terminal
		setTransitionProbability(0, 2, 1d);

		double probabilityStayNCS = modelParameters.getProbabilityOfStayingInNCS();

		setTransitionProbability(2, 2, probabilityStayNCS);
		setTransitionProbability(2, 3, 0.4 * (1 - probabilityStayNCS)); // NCS -> +M
		setTransitionProbability(2, 8, 0.4 * (1 - probabilityStayNCS)); // NCS -> -Stop
		setTransitionProbability(2, 1, 0.2 * (1 - probabilityStayNCS)); // NCS -> terminal

		setTransitionProbability(3, 4, 1d); // +M -> +CDS

		setTransitionProbability(4, 4, modelParameters.getProbabilityOfStayingInCDS()); // +CDS -> +CDS
		setTransitionProbability(4, 5, 1 - modelParameters.getProbabilityOfStayingInCDS()); // +CDS -> +Stop

		setTransitionProbability(5, 2, 1d); // +stop -> NCS always
		
		// Reverse model:
		setTransitionProbability(8, 7, 1d); // -Stop -> -CDS

		setTransitionProbability(7, 7, modelParameters.getProbabilityOfStayingInCDS()); // stay in -CDS
		setTransitionProbability(7, 6, 1 - modelParameters.getProbabilityOfStayingInCDS()); // -CDS -> -M
		setTransitionProbability(6, 2, 1d); // -M -> NCS always
		}

}
