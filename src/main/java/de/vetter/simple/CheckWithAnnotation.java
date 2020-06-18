package de.vetter.simple;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import de.vetter.masterthesis.Pair;

public class CheckWithAnnotation {
	
	private String name;
	private List<Pair<Integer, Integer>> codingSequences;
	private List<Pair<Integer, Integer>> threeUTRs;
	
	public CheckWithAnnotation(String nodeName) {
		name = nodeName;
		codingSequences = new ArrayList<Pair<Integer, Integer>>();
		threeUTRs = new ArrayList<Pair<Integer, Integer>>();
	}
	
	/**
	 * @param from sequence INDEX (starts with 0, unlike gff)
	 * @param to sequence INDEX (starts with 0, unlike gff)
	 * @param isForward
	 */
	public void addCDS(int from, int to, boolean isForward) {
		codingSequences.add(isForward ? new Pair<Integer, Integer>(from, to) : new Pair<Integer, Integer>(to, from));
	}

	/**
	 * @param from sequence INDEX (starts with 0, unlike gff)
	 * @param to sequence INDEX (starts with 0, end-exclusive!, both unlike gff)
	 * @param isForward
	 */
	public void addUTR(int from, int to, boolean isForward) {
		threeUTRs.add(isForward ? new Pair<Integer, Integer>(from, to) : new Pair<Integer, Integer>(to, from));
	}
	
	/**
	 * Go through the 3'UTRs and look for TGA at their left end (forward)/TCA at their right end (backward)
	 * @param sequence: the sequence to be checked
	 * @return
	 */
	public int countStopUGA(String sequence) {
		int count = 0;
		for(Pair<Integer, Integer> utr : threeUTRs) {
			if(utr.getFirst() < utr.getSecond()) {
				// forward
				if(sequence.substring(utr.getFirst(), utr.getFirst()+3).equals("TGA")) {
					count++;
				}
			} else {
				// reverse
				if(sequence.substring(utr.getFirst() - 3, utr.getFirst()).equals("TGA")) {
					count++;
				}
			}
		}
		return count;
	}
	
	public int getNumberOfUTRs() {
		return threeUTRs.size();
	}
	
	/**
	 * @param sequences the scaffolds.fasta (potentially uga-cleaned)
	 * @param annotation the gff which is used to assess how well the uga-cleaning worked
	 */
	public static void check(File sequences, File annotation) {
		
		
		
	}
}
