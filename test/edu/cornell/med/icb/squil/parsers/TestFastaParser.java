package edu.cornell.med.icb.squil.parsers;

import junit.framework.TestCase;

import java.io.StringReader;
import java.io.IOException;

import it.unimi.dsi.mg4j.util.MutableString;


/**
 * @author Fabien Campagne
 *         Date: Oct 25, 2006
 *         Time: 6:33:54 PM
 */
public class TestFastaParser extends TestCase {
    private String rawDescription = "Some desc. line";
    private String rawResidueCodes = "residues";

    public void testReaderSingleSequence() throws IOException {
        String input = ">" + rawDescription +
                "\n" +
                rawResidueCodes;
        StringReader reader = new StringReader(input);
        FastaParser parser = new FastaParser(reader);
        assertTrue(parser.hasNext());
        MutableString residues = new MutableString();
        MutableString description = new MutableString();
        assertFalse(parser.next(description, residues));
        assertEquals(description, new MutableString(rawDescription));
        assertEquals(residues, new MutableString(rawResidueCodes));

    }

    public void testReaderTwoSequences() throws IOException {
        String input = ">" + rawDescription +
                "\n" +
                rawResidueCodes + "\n" +
                ">" + rawDescription +
                "\n" + rawResidueCodes;
        StringReader reader = new StringReader(input);
        FastaParser parser = new FastaParser(reader);
        assertTrue(parser.hasNext());
        MutableString residues = new MutableString();
        MutableString description = new MutableString();
        assertTrue(parser.next(description, residues));
        assertEquals(description, new MutableString(rawDescription));
        assertEquals(residues, new MutableString(rawResidueCodes));
        assertFalse(parser.next(description, residues));
        assertEquals(description, new MutableString(rawDescription));
        assertEquals(residues, new MutableString(rawResidueCodes));
    }

    public void testAccessionCodeGuesses() throws IOException {
        String description1 = "P08100";
        String description2 = "P1;P08100";
        String description3 = "P1;P08100|blah";
        FastaParser parser = new FastaParser(new StringReader(""));
        MutableString accessionCode = new MutableString();
        parser.guessAccessionCode(description1, accessionCode);
        assertEquals(accessionCode, new MutableString("P08100"));
        parser.guessAccessionCode(description2, accessionCode);
        assertEquals(accessionCode, new MutableString("P08100"));
        parser.guessAccessionCode(description3, accessionCode);
        assertEquals(accessionCode, new MutableString("P08100"));
    }

    public void testResidueFilter() throws IOException {
        String rawResidues = "ABCDEFGHIKLMNPQRSTVWY-XZ";
        String rawResidues2 = "!2378273%@*#($@ABCDEFGH369<>IKLMNPQRSTVWY-XZ";
        String rawResidues3 = "12345678990/,<>!?";
        FastaParser parser = new FastaParser(new StringReader(""));
        MutableString validResidueCodes = new MutableString();
        parser.filterProteinResidues(rawResidues, validResidueCodes);
        assertEquals(validResidueCodes, new MutableString("ABCDEFGHIKLMNPQRSTVWY-XZ"));
        parser.filterProteinResidues(rawResidues2, validResidueCodes);
        assertEquals(validResidueCodes, new MutableString("ABCDEFGHIKLMNPQRSTVWY-XZ"));
        parser.filterProteinResidues(rawResidues3, validResidueCodes);
        assertEquals(validResidueCodes, new MutableString(""));

    }
}
