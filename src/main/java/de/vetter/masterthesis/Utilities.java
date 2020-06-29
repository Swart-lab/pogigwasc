package de.vetter.masterthesis;

public class Utilities {

	// Used for introns. See main in this class.
	private static final double[] LOG_FACTORIALS = new double[] { 0.0, 0.0, 0.6931471805599453, 1.791759469228055,
			3.1780538303479453, 4.787491742782046, 6.579251212010101, 8.525161361065415, 10.60460290274525,
			12.80182748008147, 15.104412573075518, 17.502307845873887, 19.98721449566189, 22.552163853123425,
			25.191221182738683, 27.899271383840894, 30.671860106080675, 33.50507345013689, 36.39544520803305,
			39.339884187199495, 42.335616460753485, 45.38013889847691, 48.47118135183523, 51.60667556776438,
			54.784729398112326, 58.003605222980525, 61.26170176100201, 64.55753862700634, 67.88974313718154,
			71.25703896716801, 74.65823634883017, 78.09222355331532, 81.55795945611504, 85.05446701758153,
			88.5808275421977, 92.13617560368711, 95.71969454214322, 99.33061245478744, 102.96819861451382,
			106.63176026064347 };
	
	// maybe rather store log factorial? yes: 21! is no longer within integer.
	
	private static final double[][] BASE_PROBABILITIES_CDS = new double[][] {
			{ 0.24, 0.13, 0.33, 0.3 }, 
			{ 0.30, 0.18, 0.37, 0.15 },
			{ 0.36, 0.12, 0.38, 0.14 } };

	private static final double[] CODON_PROBABILITIES_STOPREGION = new double[] {
			0.036, 0.024, 0.027, 0.019, 0.017, 0.003, 0.014, 0.004, 
			0.026, 0.012, 0.002, 0.015, 0.011, 0.008,     0, 0.021, 
			0.022, 0.008, 0.020, 0.005, 0.017, 0.004, 0.014, 0.001, 
			0.010, 0.007, 0.018, 0.008, 0.004, 0.002, 0.001, 0.001, 
			0.030, 0.016, 0.027, 0.022, 0.015, 0.005, 0.013, 0.003, 
			0.033, 0.017, 0.061, 0.047, 0.014, 0.008, 0.038, 0.006, 
			0.028, 0.004, 0.019, 0.011, 0.026, 0.007, 0.016, 0.003, 
			0.029, 0.012, 0.048, 0.015, 0.011, 0.006, 0.027, 0.002
	};
			
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
	
	public static int codonToIndex(String codon) {
	    int result = baseToIndex(codon.charAt(0)) * 16;
	    result += baseToIndex(codon.charAt(1)) * 4;
	    result += baseToIndex(codon.charAt(2));
		return result;
	}
	
	public static String reverseComplement(String sequence) {
		StringBuffer result = new StringBuffer();
		for(char b : sequence.toCharArray()) {
			switch(b) {
			case 'T':
			case 't':
				result.append('A');
				break;
			case 'C':
			case 'c':
				result.append('G');
				break;
			case 'A':
			case 'a': 
				result.append('T');
				break;
			case 'G':
			case 'g':
				result.append('C');
				break;
			default:
				result.append('N');
				break;
			}
		}
		return result.reverse().toString();
	}
	
	public static double getLogCodonProbabilityStopRegion(String codon) {
		return Math.log(CODON_PROBABILITIES_STOPREGION[codonToIndex(codon)]);
	}
	
	
	public static double getLogBaseProbabilityCDS(char base, int positionInCodon) {
		return Math.log(BASE_PROBABILITIES_CDS[positionInCodon][baseToIndex(base)]);
	}
	
	
	public static double logFactorial(int k) {
		if(k < LOG_FACTORIALS.length) {
			return LOG_FACTORIALS[k];
		} else {
			double result = 0;
			for(int i = 2; i <= k; i++) {
				result += Math.log(i);
			}
			return result;
		}
	}
	
	/**
	 * Precomputing the factorials
	 * @param no parameters
	 */
	public static void main(String[] no) {
		System.out.println("p(UGA)="+CODON_PROBABILITIES_STOPREGION[codonToIndex("TGA")]);
		
		double sum = 0;
		for(int i = 0; i < 64; i++) {
			sum += CODON_PROBABILITIES_STOPREGION[i];
		}
		System.out.println("Summed approx. codon usage: " + sum);
		for(int k = 0; k < 40; k++)
			System.out.print(logFactorial(k) + ", ");
	}
}
