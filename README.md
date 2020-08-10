# masterthesis
Code for Master, to work on from MPI and home

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
 3. Rightclick `de.vetter.masterthesis/App.java` and run as java application. This should give the usual help-text (set command-line args via Run>Run Configurations)