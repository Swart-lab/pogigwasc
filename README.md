# POGIGWASC
A gene predictor developed on _Loxodes magnus_ in the summer-semester of 2020 as part of the master-thesis of David Emanuel Vetter (Thesis-title: "Prediction of genes in genomes with ambiguous genetic codes"). Heavily based on Augustus.

Naming: POGIGWASC (IPA: [ˈpɔ.dʒɪdʒ.wəskʰ]) stands for "Prediciton of genes in genomes with ambiguous stop codons" -- the deviation from the thesis-name is due to the program only being specifically designed to handle an ambiguous UGA-codon. Its generaliseability to other cases of genetic code ambiguity is unknown.

# Running the program:
##Variant 1: Without eclipse -- just want to use the program
Prerequisites: **install Maven**

 1. download the code as a zip-file (from github) and extract the contents to some directory
 2. navigate to that directory (more precisely: navigate to the directory `masterthesis-master`, containing the **pom.xml**) and run the following maven command: `mvn package appassembler:assemble`
   _(this should print quite a bunch of information into the console; including (green) text informing about the number of tests run, and the number of failures encountered -- there should no be any Failures or errors; at the end, 'BUILD SUCCESS' should be printed to the console)_
 3. from the same directory, run `target/appassembler/bin/ghmm-predict`: this should print the standard help-text

To move the compiled program, copy/move the **entire** `appassembler`-directory (which can be renamed) -- it does not suffice to copy/move just the generated `target/appassembler/bin/ghmm-predict` file

##Variant 2: For working on the program
Prerequisites: **install Eclipse**

 1. clone the git-repository with eclipse (or directly with git)
 2. right-click the pom-xml and choose `Run As` > `2 Maven build...`, then enter e.g. `package appassembler:assemble` for the goals and run
   _This should produce an output as described above_
 3. Rightclick `de.vetter.pogigwasc/App.java` and run as java application. This should give the usual help-text (set command-line args via Run>Run Configurations)
 
# Program structure overview:
This project contains 2 separate parts:
 1. The GHMM/Gene-prediction in `src/main/java/de/vetter/pogigwasc`
 2. simplistic UGA-cleaning in `src/main/java/de/vetter/simple`
 
The second one is more a programming sketch, the first will be explained in slightly more detail:

Directly in `src/main/java/de/vetter/pogigwasc`, most relevant classes for gene-prediction are found, in particular, the main-class `App` with the `main`-method that defines what happens, when the program is run (cf. pom.xml). 
 - `ModelParameters` is used to load and query the model-parameters from an external file
 - `Viterbi` and `ViterbiSeed` implement the viterbi-algorithm
 - `GHMM` implements GHMMs
 - `Pair` and `GFFFeature` are very small and trivial (implementing a pair, and being an enum for writing GFF3-files)
 - The folder/subpackage `states/` contains the implementations of the states: These take care of emission probabilities and enumerating valid emission lengths

The Classes are documented in more detail in their respective files. Tests are found in `src/test/java/de/vetter/pogigwasc`, where running `AllTests` runs all tests.
