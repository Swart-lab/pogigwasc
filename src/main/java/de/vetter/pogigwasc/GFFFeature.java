package de.vetter.pogigwasc;

public enum GFFFeature {
	
	CDS("CDS"),
	INTRON("intron"),
	STOP_CODON("stop_codon");
	
	private String code;
	
	GFFFeature(String code) {
		this.code = code;
	}
	
	public String getCode() {
		return code;
	}
}
