package de.vetter.simple;

import java.io.FileReader;
import java.util.List;

/**
 * Implements the logic of a GFF-file, sorting the annotations into contiguous
 * genes (starting with a CDS and ending with a 3'UTR, potentially with Introns
 * breaking up the CDS)
 * 
 * @author dvetter
 */
public class GFFFile {
	
	private String fileName;
	private List<Gene> genes;

	public GFFFile(String gffFile) {
		// TODO implement
		read(gffFile);
	}
	
	/**
	 * Reads in the given gff and constructs an appropriate representation that can
	 * be queried for statistics like the average number of in-frame UGA per gene
	 * 
	 * @param gffFile path to a gff-file to be read.
	 */
	private void read(String gffFile) {
		// TODO implement
	}
}
