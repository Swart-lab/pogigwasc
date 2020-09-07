package de.vetter.pogigwasc;

import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import de.vetter.pogigwasc.states.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

/**
 * Main class for gene-prediction: Takes at least one, and up to three
 * parameters:
 * <ol>
 * <li>-i: the input-file (.fasta)</li>
 * <li>-o: the output file (.gff, can have arbitrary extension)</li>
 * <li>-p: the parameter-file (.properties), flagged by "-p"</li>
 * </ol>
 * 
 * If no output file is specified, one will be created in the same directory as
 * the input-file, with date and time-information in its filename.<br>
 * If no parameter-file is specified, a default parameter-file (also, and mainly
 * used for code-testing) will be used (<b>this no longer works, if the executable is removed from the sources!</b>)
 */
public class App {

	public static void main(String[] args) throws IOException, ParseException {
		LocalDateTime now = LocalDateTime.now();

		/** Setting up command-line options */
		Options commandLineOptions = new Options();
		Option help = new Option("h", "help", false, "(print this message)");
		commandLineOptions.addOption(help);
		
		Option in = new Option("i", "input", true,
				"the fasta-input-file whereon to perform gene-prediction. Supply contigs, "
						+ "not scaffolds (no unspecified nucleotides (N) allowed)");
		commandLineOptions.addOption(in);
		
		OptionGroup helpGroup = new OptionGroup();
		helpGroup.addOption(help);
		helpGroup.addOption(in);

		Option out = new Option("o", "output", true,
				"the output-file, will be in GFF3-format; if this is not specified, a .gff-file"
						+ " with a name derived from the given input-file (but also including"
						+ " information about date and time for bookkeeping) will be created "
						+ "in the same directory as the in-file");
		commandLineOptions.addOption(out);
		
		Option param = new Option("p", "parameters", true, "parameter-file in java's properties-format");
		commandLineOptions.addOption(param);
		
		CommandLine cmd = new DefaultParser().parse(commandLineOptions, args, true);
		
		/** Generally nice behaviour: if one simply runs the program without sufficient parameters: provide help */
		if(cmd.hasOption('h') || !cmd.hasOption('i')) {
			new HelpFormatter().printHelp("... -i infile.fasta [-o outfile] [-p parameterfile.properties]", commandLineOptions);
			return;
		}

		
		/** File management */

		File input = new File(cmd.getOptionValue('i'));
		BufferedReader reader = new BufferedReader(new FileReader(input));

		File output;
		BufferedWriter writer;
		if (cmd.hasOption('o')) {
			output = new File(cmd.getOptionValue('o'));
		} else {
			// cut off the .fasta (or something like .fa)
			String genericOutputFilename = cmd.getOptionValue('i').substring(0, cmd.getOptionValue('i').lastIndexOf("."));
			genericOutputFilename += "-day" + DateTimeFormatter.ofPattern("yyyy-MM-dd_HHmmss").format(now);
			genericOutputFilename += ".predicted.gff";

			output = new File(genericOutputFilename);
		}

		writer = new BufferedWriter(new FileWriter(output));

		File parameterFile;
		if (cmd.hasOption('p')) {
			parameterFile = new File(cmd.getOptionValue('p'));
		} else {
			System.out.println("Trying to use default parameter-file:\n" + "THIS IS NOT RECOMMENDED\n"
					+ "and only works, if the respective file is in the correct relative path.\n"
					+ "The file will not be particularly fit to the specific prediction!");
			parameterFile = new File(
					"resources\\de\\vetter\\masterthesis\\parameter\\parameters-examplefile.properties");
		}
		 
		final ModelParameters modelParameters = new ModelParameters(new FileReader(parameterFile));

		
		GHMM ghmm = new LoxodesMagnusGHMM(modelParameters);

		System.out.println(ghmm);

		/** Write some book-keeping information into the output-file */
		writer.write("##gff-version 3");
		writer.newLine();
		writer.write("##Generated on: " + DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss").format(now));
		writer.newLine();
		writer.write("##Source file: " + input.getAbsolutePath());
		writer.newLine();
		writer.write("##Parameter file: " + parameterFile.getAbsolutePath());
		writer.newLine();

		/**
		 * Now read in the fasta file sequence for sequence, and for each sequence
		 * generate a gff-output
		 */
		String currentHeader = null;
		String currentSequence = null;

		for (String line = reader.readLine(); line != null; line = reader.readLine()) {
			if (line.startsWith(">")) {
				if (currentHeader != null) {
					doPredictions(ghmm, writer, currentHeader, currentSequence, modelParameters);
				}
				currentHeader = line.substring(1);
				currentSequence = "";
			} else {
				currentSequence += line;
			}
		}

		// on last sequence
		if (currentHeader != null) {
			doPredictions(ghmm, writer, currentHeader, currentSequence, modelParameters);
		}

		reader.close();
		writer.flush();
		writer.close();

		System.out.println("______________________\nWrote to output-file " + output.getAbsolutePath());
	}

	/**
	 * Performs gene prediction on the given sequence and writes the result as gff
	 * into the writer.
	 * 
	 * @param ghmm
	 * @param writer
	 * @param currentHeader
	 * @param currentSequence
	 * @throws IOException              if something goes wrong while trying to
	 *                                  write into the out-file
	 * @throws IllegalArgumentException if the sequence contains N, not just TGAC
	 */
	public static void doPredictions(GHMM ghmm, BufferedWriter writer, String currentHeader, String currentSequence,
			ModelParameters parameters) throws IOException {
		
		System.out.println("\nPredicting genes in " + currentHeader + ":\n");
		if (currentSequence.contains("N")) { // or three?
			throw new IllegalArgumentException(
					"The given sequence contained uncharacterised Nucleotides (N). Please provide contigs, not scaffolds");
		}

		Viterbi viterbi = new Viterbi(ghmm, currentSequence);
		viterbi.setAbbreviating(true);
		for (Parse parse : viterbi.computeParses()) {
			System.out.println("\n\tWriting parse to file");
			writer.write(LoxodesMagnusGHMM.parseToGFF(currentHeader, parse, parameters));
		}

		System.out.println("done");
	}
}
