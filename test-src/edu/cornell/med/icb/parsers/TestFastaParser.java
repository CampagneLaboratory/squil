package edu.cornell.med.icb.parsers;

import it.unimi.dsi.mg4j.util.MutableString;
import junit.framework.TestCase;

import java.io.IOException;
import java.io.StringReader;


/**
 * @author Fabien Campagne
 *         Date: Oct 25, 2006
 *         Time: 6:33:54 PM
 */
public class TestFastaParser extends TestCase {
    private final String rawDescription = "Some desc. line";
    private final String rawResidueCodes = "residues";

    public void testReaderSingleSequence() throws IOException {
        final String input = ">" + rawDescription +
                "\n" +
                rawResidueCodes;
        final StringReader reader = new StringReader(input);
        final FastaParser parser = new FastaParser(reader);
        assertTrue(parser.hasNext());
        final MutableString residues = new MutableString();
        final MutableString description = new MutableString();
        assertFalse(parser.next(description, residues));
        assertEquals(description, new MutableString(rawDescription));
        assertEquals(residues, new MutableString(rawResidueCodes));

    }

    public void testReaderTwoSequences() throws IOException {
        final String input = ">" + rawDescription +
                "\n" +
                rawResidueCodes + "\n" +
                ">" + rawDescription +
                "\n" + rawResidueCodes;
        final StringReader reader = new StringReader(input);
        final FastaParser parser = new FastaParser(reader);
        assertTrue(parser.hasNext());
        final MutableString residues = new MutableString();
        final MutableString description = new MutableString();
        assertTrue(parser.next(description, residues));
        assertEquals(description, new MutableString(rawDescription));
        assertEquals(residues, new MutableString(rawResidueCodes));
        assertFalse(parser.next(description, residues));
        assertEquals(description, new MutableString(rawDescription));
        assertEquals(residues, new MutableString(rawResidueCodes));
    }

    public void testAccessionCodeGuesses() throws IOException {
        final String description1 = "P08100";
        final String description2 = "P1;P08100";
        final String description3 = "P1;P08100|blah";
        final FastaParser parser = new FastaParser(new StringReader(""));
        final MutableString accessionCode = new MutableString();
        parser.guessAccessionCode(description1, accessionCode);
        assertEquals(accessionCode, new MutableString("P08100"));
        parser.guessAccessionCode(description2, accessionCode);
        assertEquals(accessionCode, new MutableString("P08100"));
        parser.guessAccessionCode(description3, accessionCode);
        assertEquals(accessionCode, new MutableString("P08100"));
    }

    public void testResidueFilter() throws IOException {
        final String rawResidues = "ABCDEFGHIKLMNPQRSTVWY-XZ";
        final String rawResidues2 = "!2378273%@*#($@ABCDEFGH369<>IKLMNPQRSTVWY-XZ";
        final String rawResidues3 = "12345678990/,<>!?";
        final FastaParser parser = new FastaParser(new StringReader(""));
        final MutableString validResidueCodes = new MutableString();
        parser.filterProteinResidues(rawResidues, validResidueCodes);
        assertEquals(validResidueCodes, new MutableString("ABCDEFGHIKLMNPQRSTVWY-XZ"));
        parser.filterProteinResidues(rawResidues2, validResidueCodes);
        assertEquals(validResidueCodes, new MutableString("ABCDEFGHIKLMNPQRSTVWY-XZ"));
        parser.filterProteinResidues(rawResidues3, validResidueCodes);
        assertEquals(validResidueCodes, new MutableString(""));

    }
}
