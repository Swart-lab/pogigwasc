package de.vetter.masterthesis;


import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

/**
 * Test-Suite-class collecting all the test-cases.
 * 
 * @author David Emanuel Vetter
 */
@RunWith(value=Suite.class)
@SuiteClasses(value={TestGHMM.class, TestViterbi.class, TestPair.class, TestModelParameters.class, TestUtilities.class, TestParse.class, TestParseToGFF.class})
public class AllTests {

}
