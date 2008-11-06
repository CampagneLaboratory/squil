/*
 * Copyright (C) 2006-2008 Institute for Computational Biomedicine,
 *                         Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cornell.med.icb.parsers;

import junit.framework.TestCase;

import java.io.IOException;
import java.io.StringReader;

import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.NullReader;

/**
 * Validates the functionality of the {@link edu.cornell.med.icb.parsers.FastaParser} class.
 *
 * @author Fabien Campagne
 *         Date: Oct 25, 2006
 *         Time: 6:33:54 PM
 */
public class TestFastaParser extends TestCase {
    /**
     * Sample description line from a FASTA sequence file.
     */
    private static final String RAW_DESCRIPTION = "Some desc. line";
    /**
     * Sample residue codes from a FASTA sequence file.
     */
    private static final String RAW_RESIDUE_CODES = "residues";

    /**
     * Validates that the FASTA parser can handle input containing a single sequence.
     * @throws IOException if there is a problem with the reader
     */
    public void testReaderSingleSequence() throws IOException {
        final String input = ">" + RAW_DESCRIPTION + "\n" + RAW_RESIDUE_CODES;
        final StringReader reader = new StringReader(input);
        final FastaParser parser = new FastaParser(reader);
        assertTrue(parser.hasNext());
        final MutableString residues = new MutableString();
        final MutableString description = new MutableString();
        assertFalse(parser.next(description, residues));
        assertEquals(new MutableString(RAW_DESCRIPTION), description);
        assertEquals(new MutableString(RAW_RESIDUE_CODES), residues);
        assertFalse(parser.hasNext());
    }

    /**
     * Validates that the FASTA parser can handle input containing tow single sequences.
     * @throws IOException if there is a problem with the reader
     */
    public void testReaderTwoSequences() throws IOException {
        final String input = ">" + RAW_DESCRIPTION +
                "\n" +
                RAW_RESIDUE_CODES + "\n" +
                ">" + RAW_DESCRIPTION +
                "\n" + RAW_RESIDUE_CODES;
        final StringReader reader = new StringReader(input);
        final FastaParser parser = new FastaParser(reader);
        assertTrue(parser.hasNext());
        final MutableString residues = new MutableString();
        final MutableString description = new MutableString();
        assertTrue(parser.next(description, residues));
        assertEquals(new MutableString(RAW_DESCRIPTION), description);
        assertEquals(new MutableString(RAW_RESIDUE_CODES), residues);
        assertFalse(parser.next(description, residues));
        assertEquals(new MutableString(RAW_DESCRIPTION), description);
        assertEquals(new MutableString(RAW_RESIDUE_CODES), residues);
        assertFalse(parser.hasNext());
    }

    /**
     * Validates that the FASTA parser extracts accession codes from description lines properly.
     * @throws IOException if there is a problem with the reader
     */
    public void testAccessionCodeGuesses() throws IOException {
        final String description1 = "P08100";
        final String description2 = "P1;P08100";
        final String description3 = "P1;P08100|blah";
        final MutableString accessionCode = new MutableString();

        FastaParser.guessAccessionCode(description1, accessionCode);
        assertEquals(new MutableString("P08100"), accessionCode);

        FastaParser.guessAccessionCode(description2, accessionCode);
        assertEquals(new MutableString("P08100"), accessionCode);

        FastaParser.guessAccessionCode(description3, accessionCode);
        assertEquals(new MutableString("P08100"), accessionCode);

        assertEquals("ASDF", FastaParser.guessAccessionCode("ASDF"));
        assertEquals("P1:P08100", FastaParser.guessAccessionCode("P1:P08100"));
        assertEquals("P2;P08100", FastaParser.guessAccessionCode("P2;P08100"));
        assertEquals("", FastaParser.guessAccessionCode(" P08100"));
        assertEquals("", FastaParser.guessAccessionCode("\tP08100"));
        assertEquals("", FastaParser.guessAccessionCode("|P08100"));
        assertEquals("", FastaParser.guessAccessionCode(""));
    }

    /**
     * Validates that the FASTA parser can extract and filter residue codes properly.
     * @throws IOException if there is a problem with the reader
     */
    public void testResidueFilter() throws IOException {
        final MutableString validResidueCodes = new MutableString();

        // all the residues here are valid, so nothing should change
        final String rawResidues = "ABCDEFGHIKLMNPQRSTVWY-XZ";
        FastaParser.filterProteinResidues(rawResidues, validResidueCodes);
        assertEquals(new MutableString(rawResidues), validResidueCodes);

        // some characters are invalid
        final String rawResidues2 = "!2378273%@*#($@ABCDEFGH369<>IKLMNPQRSTVWY-XZ";
        FastaParser.filterProteinResidues(rawResidues2, validResidueCodes);
        assertEquals(new MutableString(rawResidues), validResidueCodes);

        // No residues in this string are valid
        final String rawResidues3 = "12345678990/,<>!?";
        FastaParser.filterProteinResidues(rawResidues3, validResidueCodes);
        assertEquals(new MutableString(""), validResidueCodes);

        // The "." character should be replaced by a "-"
        final String rawResidues4 = "FACE.BEEF";
        assertEquals("FACE-BEEF", FastaParser.filterProteinResidues(rawResidues4));    
    }

    /**
     * Validates the functionality of the
     * {@link edu.cornell.med.icb.parsers.FastaParser#next(it.unimi.dsi.lang.MutableString,
     * it.unimi.dsi.lang.MutableString)} and
     * {@link edu.cornell.med.icb.parsers.FastaParser#hasNext()} methods.
     * @throws IOException if there is a problem with the reader
     */
    public void testNullReader() throws IOException {
        final FastaParser parser = new FastaParser();
        parser.setReader(new FastBufferedReader(NullReader.getInstance()));
        assertFalse(parser.hasNext());
        assertFalse(parser.next(null, null));
    }

    /**
     * Validates that the FASTA parser can handle input containing only residue codes.
     * @throws IOException if there is a problem with the reader
     */
    public void testReaderNoDescriptionSequence() throws IOException {
        final FastaParser parser = new FastaParser(new StringReader("\n" + RAW_RESIDUE_CODES));
        assertFalse(parser.hasNext());
    }
}
