package de.vetter.pogigwasc;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.junit.Before;
import org.junit.Test;

public class TestParseToGFF {

	private GHMM ghmm;
	private ModelParameters mp;
	
	@Before
	public void setup() throws IOException {
		FileReader reader = new FileReader(
				new File("resources//de//vetter//pogigwasc//parameter//parameters-examplefile.properties"));
		mp = new ModelParameters(reader);
		ghmm = new LoxodesMagnusGHMM(mp);
	}
	
	@Test
	public void testParseToGFFSingleExonForward() {
		Parse parse = new Parse();
		for (int i = 0; i < 10; i++)
			parse.add(ghmm.getState(2), 1); // Ten noncoding bases (NCS)
		
		parse.add(ghmm.getState(3), 6); // 6nt: start-region => bases 14, 15, 16 are start-ATG
		
		for (int i = 0; i < 7; i++) // seven codons of CDS
			parse.add(ghmm.getState(4), 3); 
		
		parse.add(ghmm.getState(5), 24);
		
		for (int i = 0; i < 5; i++)
			parse.add(ghmm.getState(2), 1); // five noncoding bases (NCS)
		
		String gff = LoxodesMagnusGHMM.parseToGFF("CONTIGNAME", parse, mp);
		
		String[] lines = gff.split("\n");
		String[] cdsLine = lines[0].split("\t");
		String[] stopCodonLine = lines[1].split("\t");
		
		assertEquals("CONTIGNAME", cdsLine[0]);
		assertEquals("CONTIGNAME", stopCodonLine[0]);
		
		assertEquals("CDS", cdsLine[2]);
		assertEquals(14, Integer.parseInt(cdsLine[3]));
		assertEquals(58, Integer.parseInt(cdsLine[4]));
		assertEquals("+", cdsLine[6]);
		
		assertEquals("stop_codon", stopCodonLine[2]);
		assertEquals(59, Integer.parseInt(stopCodonLine[3]));
		assertEquals(61, Integer.parseInt(stopCodonLine[4]));
		assertEquals("+", stopCodonLine[6]);
	}
	
	@Test
	public void testParseToGFFSingleExonReverse() {
		Parse parse = new Parse();
		for (int i = 0; i < 10; i++)
			parse.add(ghmm.getState(2), 1); // Ten noncoding bases (NCS)
		
		parse.add(ghmm.getState(15), 24); // stop-region (reverse) -> stop-codon is 11:T,12:C,13:A
		
		for (int i = 0; i < 7; i++) // seven codons of CDS
			parse.add(ghmm.getState(14), 3); 
		
		parse.add(ghmm.getState(13), 6); // Start region -> start-codon is 56,57,58
		
		for (int i = 0; i < 5; i++)
			parse.add(ghmm.getState(2), 1); // five noncoding bases (NCS)
		
		String gff = LoxodesMagnusGHMM.parseToGFF("CONTIGNAME", parse, mp);
		
		String[] lines = gff.split("\n");
		String[] stopCodonLine = lines[0].split("\t");
		String[] cdsLine = lines[1].split("\t");
		
		assertEquals("CONTIGNAME", cdsLine[0]);
		assertEquals("CONTIGNAME", stopCodonLine[0]);
		
		assertEquals("CDS", cdsLine[2]);
		assertEquals(14, Integer.parseInt(cdsLine[3]));
		assertEquals(58, Integer.parseInt(cdsLine[4]));
		assertEquals("-", cdsLine[6]);
		
		assertEquals("stop_codon", stopCodonLine[2]);
		assertEquals(11, Integer.parseInt(stopCodonLine[3]));
		assertEquals(13, Integer.parseInt(stopCodonLine[4]));
		assertEquals("-", stopCodonLine[6]);
	}
	
	
	@Test
	public void testParseToGFFSingleIntronForward() {
		Parse parse = new Parse();
		for (int i = 0; i < 20; i++)
			parse.add(ghmm.getState(2), 1); 
		
		parse.add(ghmm.getState(3), 6); 
		
		// initial exon
		for (int i = 0; i < 7; i++)
			parse.add(ghmm.getState(4), 3);
		
		parse.add(ghmm.getState(10), 2); // +2L
		parse.add(ghmm.getState(11), 17); // +intron 2-1
		parse.add(ghmm.getState(12), 1); // +1R
		
		// second exon
		for (int i = 0; i < 5; i++)
			parse.add(ghmm.getState(4), 3);

		// stop-region
		parse.add(ghmm.getState(5), 24);
		
		for (int i = 0; i < 15; i++)
			parse.add(ghmm.getState(2), 1); // fifteen noncoding bases (NCS)
		
		String gff = LoxodesMagnusGHMM.parseToGFF("CONTIGNAME", parse, mp);
		
		String[] lines = gff.split("\n");
		String[] exon1 = lines[0].split("\t");
		String[] intron = lines[1].split("\t");
		String[] exon2 = lines[2].split("\t");
		String[] stopCodonLine = lines[3].split("\t");
		
		assertEquals("CONTIGNAME", exon1[0]);
		assertEquals("CONTIGNAME", intron[0]);
		assertEquals("CONTIGNAME", exon2[0]);
		assertEquals("CONTIGNAME", stopCodonLine[0]);
		
		assertEquals("CDS", exon1[2]);
		assertEquals("CDS", exon2[2]);
		assertEquals(24, Integer.parseInt(exon1[3]));
		assertEquals(49, Integer.parseInt(exon1[4]));
		assertEquals("+", exon1[6]);
		
		assertEquals(50, Integer.parseInt(intron[3]));
		assertEquals(66, Integer.parseInt(intron[4]));
		assertEquals("+", intron[6]);
		
		assertEquals(67, Integer.parseInt(exon2[3]));
		assertEquals(103, Integer.parseInt(exon2[4]));
		assertEquals("+", exon2[6]);
		
		assertEquals("stop_codon", stopCodonLine[2]);
		assertEquals(104, Integer.parseInt(stopCodonLine[3]));
		assertEquals(106, Integer.parseInt(stopCodonLine[4]));
		assertEquals("+", stopCodonLine[6]);
	}

}
