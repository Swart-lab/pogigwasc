package de.vetter.masterthesis;


import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(value=Suite.class)
@SuiteClasses(value={TestGHMM.class, TestViterbi.class, TestPair.class, TestModelParameters.class, TestUtilities.class})
public class AllTests {

}
