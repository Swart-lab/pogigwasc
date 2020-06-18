package de.vetter.simple;

public class ThreeUTR {

	private boolean forward;
	private int start, end;
	
	/**
	 * 
	 * @param start start-INDEX (starts with 0 unlike for GFF)
	 * @param end exclusive end-INDEX (unlike GFF, no change needed though)
	 * @param forward
	 */
	public ThreeUTR(int start, int end, boolean forward) {
		this.start = start;
		this.end = end;
		this.forward = forward;
	}
	
	/**
	 * @return the leftmost index of the stop-UGA or UCA (i.e. stop is found at result:result+3)
	 */
	public int getStopPosition() {
		if(forward) {
			return start;
		} else {
			return end-3;
		}
	}
	
	public boolean isForward() {
		return forward;
	}
	
}
