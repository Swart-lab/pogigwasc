package de.vetter.masterthesis.states;

public class StopRegionState extends HMMState {
	private boolean strand;
	
	public StopRegionState(String name, boolean strand) { 
		super(name);
		this.strand = strand;
	}

	@Override
	public double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission) {
		// TODO: instead: have this combine the prestop 21 bases (no UGA) and the stop-UGA, with a total fixed length of 24 nt.
		
		
		return newEmission == "TGA" ? 0 : Double.NEGATIVE_INFINITY;
	}

}
