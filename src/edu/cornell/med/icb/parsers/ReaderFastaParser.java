/*
 * Copyright (C) 2007-2009 Institute for Computational Biomedicine,
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

import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.lang.MutableString;

import java.io.IOException;
import java.io.Reader;

/**
 * A parser for fasta files which supports reading sequences one character at a time. Useful for processing
 * sequences containing entire chromosomes.
 *
 * @author Fabien Campagne
 *         Date: Nov 10, 2008
 *         Time: 6:09:57 PM
 */
public final class ReaderFastaParser {
    /**
     * The reader used to parse the FASTA sequences.
     */
    private FastBufferedReader reader;

    /**
     * Indicates whether or not the reader contains more sequences.
     */
    private boolean hasNext;

    /**
     * The list of valid protein sequence characters.
     */
    private static final String VALID_PROTEIN_RESIDUES = "ABCDEFGHIKLMNPQRSTVWY-XZ";

    /**
     * Used to store lines read from the FASTA sequence reader.
     */
    private MutableString line = new MutableString();

    /**
     * The description line previously seen from the FASTA sequence reader.
     */
    private MutableString previousDescriptionLine = new MutableString();


    /**
     * Create a parser to read sequences.
     *
     * @param fastaFileSource The reader over the FASTA formatted data.
     * @throws java.io.IOException if the sequence cannot be read using the reader
     */
    public ReaderFastaParser(final Reader fastaFileSource) throws IOException {
        this();
        setReader(fastaFileSource);
    }

    /**
     * Create a parser to read sequences.
     */
    public ReaderFastaParser() {
        super();
    }

    /**
     * Repositions this reader on a different file/data content.
     *
     * @param reader the new reader to use to parse the file
     * @throws java.io.IOException if the sequence cannot be read using the reader
     */
    public void setReader(final Reader reader) throws IOException {
        if (reader instanceof FastBufferedReader) {
            this.reader = (FastBufferedReader) reader;
        } else {
            this.reader = new FastBufferedReader(reader);
        }
        hasNext = readNextDescriptionLine(this.reader);
    }

    /**
     * Returns true if the reader has at least one more sequence.
     *
     * @return True if a call to next will return another sequence.
     */
    public boolean hasNextSequence() {
        return hasNext;
    }

    /**
     * Obtain the next sequence from the reader over the FASTA formatted content.
     * This method returns true until there is no more sequence to parse in the input. When the
     * method returns false, the content of the parameters descriptionLine and residues is
     * unspecified.
     *
     * @param descriptionLine Where the raw description line will be written.
     * @return True if hasNext() is true, False otherwise.
     * @throws java.io.IOException if there is a problem reading from the input
     */
    public boolean nextSequence(final MutableString descriptionLine) throws IOException {
        if (!hasNext) {
            return false;
        } else {
            descriptionLine.replace(previousDescriptionLine);
            return true;

        }
    }

    public Reader getBaseReader() {
        return new OneBaseAtATimeReader(reader);
    }

    /**
     * Try to extract an accession code from a FASTA description line.
     *
     * @param descriptionLine The line of text to parse for the accession code
     * @param accessionCode   The location to place the resulting accession code
     */
    public static void guessAccessionCode(final CharSequence descriptionLine,
                                          final MutableString accessionCode) {
        accessionCode.setLength(0);
        final int startIndex;
        if (descriptionLine.length() > 3 && descriptionLine.charAt(0) == 'P'
                && descriptionLine.charAt(1) == '1' && descriptionLine.charAt(2) == ';') {
            startIndex = 3;
        } else {
            startIndex = 0;
        }

        for (int i = startIndex; i < descriptionLine.length(); i++) {
            final char c = descriptionLine.charAt(i);
            if (c == ' ' || c == '\t' || c == '|') {
                break;
            }
            accessionCode.append(c);
        }
    }


    /**
     * Reade the sequence until the next description line or end of file is found.
     *
     * @param fastBufferedReader The reader object to get lines from
     * @return true if another description line was found, false otherwise
     * @throws java.io.IOException if there was a problem with the reader
     */
    private boolean readNextDescriptionLine(final FastBufferedReader fastBufferedReader)
            throws IOException {
        for (; ;) {
            // loop until a line that starts with > if found, or the end of file is reached.
            previousDescriptionLine = fastBufferedReader.readLine(previousDescriptionLine);
            if (previousDescriptionLine == null) {
                return false;
            }
            if (previousDescriptionLine.startsWith(">")) { // NOPMD charAt fails on an empty string
                previousDescriptionLine = removeBracket(previousDescriptionLine);
                return true;
            }
        }
    }

    /**
     * Removes the bracket character (">") from the description line.
     *
     * @param descriptionLine The line to remove the bracket from
     * @return The resulting line without the bracket
     */
    private MutableString removeBracket(final MutableString descriptionLine) {
        return descriptionLine.substring(1, descriptionLine.length());
    }


    private class OneBaseAtATimeReader extends Reader {
        private FastBufferedReader sequenceReader;

        public OneBaseAtATimeReader(FastBufferedReader reader) {
            sequenceReader = reader;
        }

        public int read() throws IOException {
            int c = sequenceReader.read();
            if (c == -1) {
                hasNext = false;   // no more sequences, we just found the end of file.
                return -1;
            } else {
                char character = (char) c;
                switch (character) {

                    case 'A':
                    case 'C':
                    case 'T':
                    case 'G':
                    case 'N':
                    case 'U':
                        return c;
                    case '>':
                        line.setLength(0);

                        sequenceReader.readLine(line);
                        previousDescriptionLine.replace(line);
                        hasNext = true;
                        return -1; // end of this specific sequence.
                    default:
                        return read();
                }
            }
        }

        public int read(char[] chars, int i, int i1) throws IOException {
            throw new UnsupportedOperationException();
        }

        public void close() throws IOException {
            // do nothing, we do not own the undelying reader.
        }
    }
}