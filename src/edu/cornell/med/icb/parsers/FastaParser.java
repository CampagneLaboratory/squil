/*
 * Copyright (C) 2007-2008 Institute for Computational Biomedicine,
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
 * Parse a FASTA file. In contrast to the
 * <a href=http://icb.med.cornell.edu/apidoc/crover/edu/mssm/crover/imports/FastaReader.html">
 * crover FastaParser</a>, this class uses {@link it.unimi.dsi.lang.MutableString} for efficiency,
 * and loads sequences lazily. This means that clients can start processing sequences in a file
 * before the file is completely loaded. This reader can therefore process very large files without
 * consuming more memory than is needed to process the largest sequence in the file.
 *
 * @author Fabien Campagne
 *         Date: Oct 25, 2006
 *         Time: 6:09:57 PM
 */
public final class FastaParser {
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
     * @throws IOException if the sequence cannot be read using the reader
     */
    public FastaParser(final Reader fastaFileSource) throws IOException {
        this();
        setReader(fastaFileSource);
    }

    /**
     * Create a parser to read sequences.
     */
    public FastaParser() {
        super();
    }

    /**
     * Repositions this reader on a different file/data content.
     *
     * @param reader the new reader to use to parse the file
     * @throws IOException if the sequence cannot be read using the reader
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
    public boolean hasNext() {
        return hasNext;
    }

    /**
     * Obtain the next sequence from the reader over the FASTA formatted content.
     * This method returns true until there is no more sequence to parse in the input. When the
     * method returns false, the content of the parameters descriptionLine and residues is
     * unspecified.
     *
     * @param descriptionLine Where the raw description line will be written.
     * @param residues When the raw residue lines will be written.
     * @return True if hasNext() is true, False otherwise.
     * @throws IOException if there is a problem reading from the input
     */
    public boolean next(final MutableString descriptionLine,
                        final MutableString residues) throws IOException {
        if (!hasNext) {
            return false;
        } else {
            descriptionLine.replace(previousDescriptionLine);
            // read residues:
            return readResidues(residues);
        }
    }

    /**
     * Try to extract an accession code from a FASTA description line.
     *
     * @param descriptionLine The line of text to parse for the accession code
     * @param accessionCode The location to place the resulting accession code
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
     * Filter a string to keep only protein residues.
     *
     * @param rawResidues A string that may contain any character.
     * @param filteredResidues The subset of characters that represent valid protein residue
     * codes, in the order in which they occur in the rawResidue string.
     */
    public static void filterProteinResidues(final CharSequence rawResidues,
                                                      final MutableString filteredResidues) {
        filteredResidues.setLength(0);
        for (int i = 0; i < rawResidues.length(); i++) {
            char residueCode = rawResidues.charAt(i);
            if (residueCode == '.') {
                residueCode = '-';
            } else {
                residueCode = Character.toUpperCase(residueCode);
            }

            // add the residue code only if it is valid
            if (VALID_PROTEIN_RESIDUES.indexOf(residueCode) != -1) {
                filteredResidues.append(residueCode);
            }
        }
    }

    /**
     * Reade the sequence until the next description line or end of file is found.
     * @param fastBufferedReader The reader object to get lines from
     * @return true if another description line was found, false otherwise
     * @throws IOException if there was a problem with the reader
     */
    private boolean readNextDescriptionLine(final FastBufferedReader fastBufferedReader)
            throws IOException {
        for (;;) {
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
     * @param descriptionLine The line to remove the bracket from
     * @return The resulting line without the bracket
     */
    private MutableString removeBracket(final MutableString descriptionLine) {
        return descriptionLine.substring(1, descriptionLine.length());
    }

    /**
     * Read and store the residues from the current sequence.
     * @param residues The object to store the residue sequence into.
     * @return true if there are more sequences left after getting the current sequence
     * @throws IOException if there is a problem reading the sequence
     */
    private boolean readResidues(final MutableString residues) throws IOException {
        residues.setLength(0);
        line.setLength(0);
        for (;;) {
            line = reader.readLine(line);
            if (line == null) {
                hasNext = false;
                return hasNext;
            }
            if (line.startsWith(">")) {   // NOPMD charAt fails on an empty string
                previousDescriptionLine.replace(line);
                previousDescriptionLine = removeBracket(previousDescriptionLine);
                return hasNext;
            }
            residues.append(line);
        }
    }
}
