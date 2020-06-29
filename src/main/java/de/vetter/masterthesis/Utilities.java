package de.vetter.masterthesis;

public class Utilities {

	// Used for introns.
	private static final double[] logFactorials = new double[] {0, 0}; // TODO
	
	// maybe rather store log factorial? yes: 21! is no longer within integer.
	
	/**
	 * T-C-A-G
	 * @param base
	 * @return
	 */
	public static int baseToIndex(char base) {
		switch(base) {
		case 'T':
		case 't':
			return 0;
		case 'C':
		case 'c':
			return 1;
		case 'A':
		case 'a': 
			return 2;
		case 'G':
		case 'g':
			return 3;
		default:
			return -1;
		}
	}
	
	
	public static double logFactorial(int k) {
		if(k < logFactorials.length) {
			return logFactorials[k];
		} else {
			double result = 0;
			for(int i = 2; i <= k; i++) {
				result += Math.log(i);
			}
			return result;
		}
	}
}
