package de.vetter.masterthesis;

/**
 * Abstract class to generalise over different types of states (mainly
 * differentiated by their emmission distributions)
 * 
 * @author David Emanuel Vetter
 *
 */
public abstract class HMMState {
	
	private String name; // mainly for output/human-readability
	
	public HMMState(String name) {
		this.name = name;
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	
	/**
	 * Computes the log_2 of the probability of emitting the newEmission when in
	 * current State, having before been in previousState and having (from the
	 * start) previously emitted emissionHistory.
	 * 
	 * To be implemented by concrete instances (probably ad-hoc-inheritors)
	 * 
	 * @param previousState
	 * @param emissionHistory
	 * @param newEmission
	 * @return log_2 of e as specified by Stanke.
	 */
	public abstract double computeLogEmissionProbability(int previousState, String emissionHistory, String newEmission);

}
