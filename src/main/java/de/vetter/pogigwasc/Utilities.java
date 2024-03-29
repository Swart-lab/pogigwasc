package de.vetter.pogigwasc;

/**
 * Collection of some general methods with no particular place. This class implements the base-index-correspondence:
 * TCAG<br>
 * (i.e. T has index 0, C index 1 etc)
 * 
 * @author David Emanuel Vetter
 */
public class Utilities {

	// Used for Poisson
	private static final double[] LOG_FACTORIALS = new double[] { 0.0, 0.0, 0.6931471805599453, 1.791759469228055,
			3.1780538303479453, 4.787491742782046, 6.579251212010101, 8.525161361065415, 10.60460290274525,
			12.80182748008147, 15.104412573075518, 17.502307845873887, 19.98721449566189, 22.552163853123425,
			25.191221182738683, 27.899271383840894, 30.671860106080675, 33.50507345013689, 36.39544520803305,
			39.339884187199495, 42.335616460753485, 45.38013889847691, 48.47118135183523, 51.60667556776438,
			54.784729398112326, 58.003605222980525, 61.26170176100201, 64.55753862700634, 67.88974313718154,
			71.25703896716801, 74.65823634883017, 78.09222355331532, 81.55795945611504, 85.05446701758153,
			88.5808275421977, 92.13617560368711, 95.71969454214322, 99.33061245478744, 102.96819861451382,
			106.63176026064347 };
			
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
	
	/**
	 * @param sequence (preferably <i>not</i> containing N
	 * @return reverse complement of that sequence -- if some other base than the
	 *         standard TCAG is found in the sequence, its complement is will be N.
	 */
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
	
	/**
	 * @param k
	 * @return ln(k!)
	 */
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

	
	public static double sumVector(double[] vector) {
		double sum = 0.0;
		for(double d : vector) {
			if(d < 0 || d > 1)
				throw new IllegalArgumentException("Probability vector contained a non-probability value " + d);
			sum += d;
		}
		return sum;
	}
	
	/**
	 * Throws an Exception of the vector is not normalised (does not sum to 1) with
	 * given tolerance (double-imprecision), and given error-message appended with
	 * the actual sum
	 * 
	 * @param vector
	 * @param tolerance
	 * @param message
	 */
	public static void assertNormalisation(double[] vector, double tolerance, String message) {
		if(Math.abs(sumVector(vector) - 1) > tolerance) {
			throw new IllegalArgumentException(message + sumVector(vector));
		}
	}
	
	
	/**
	 * Precomputing the factorials
	 * @param no parameters
	 */
	public static void main(String[] no) {
		for(int k = 0; k < 40; k++)
			System.out.print(logFactorial(k) + ", ");
	}
	
	
}
