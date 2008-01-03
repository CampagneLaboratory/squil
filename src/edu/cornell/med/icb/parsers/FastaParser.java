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

import it.unimi.dsi.mg4j.io.FastBufferedReader;
import it.unimi.dsi.mg4j.util.MutableString;

import java.io.IOException;
import java.io.Reader;

/**
 * Parse a FASTA file. In contrast to the <a href=http://icb.med.cornell.edu/apidoc/crover/edu/mssm/crover/imports/FastaReader.html">crover FastaParser</a>,
 * this class uses {@link MutableString} for efficiency, and loads sequences lazily. This means that clients can start processing sequences in a file before the
 * file is completely loaded. This reader can therefore process very large files without consuming more memory than is needed to process the largest sequence
 * in the file.
 *
 * @author Fabien Campagne
 *         Date: Oct 25, 2006
 *         Time: 6:09:57 PM
 */
public final class FastaParser {
    private FastBufferedReader reader;
    private boolean hasNext;
    private final String validProteinResidues = "ABCDEFGHIKLMNPQRSTVWY-XZ";
    private MutableString line = new MutableString();

    /**
     * Create a parser to read sequences.
     *
     * @param fastaFileSource The reader over the Fasta-formatted data.
     */
    public FastaParser(final Reader fastaFileSource) throws IOException {
        this();
        setReader(fastaFileSource);
        hasNext = readNextDescriptionLine(reader);
    }

    public FastaParser() {
        super();
        previousDescriptionLine = new MutableString();
    }

    /**
     * Repositions this reader on a different file/data content.
     *
     * @param reader
     */
    public void setReader(final Reader reader) {
        this.reader = new FastBufferedReader(reader);
    }

    /**
     * Returns true if the reader has at least one more sequence.
     *
     * @return True if a call to next will return another sequence.
     */
    public boolean hasNext() {
        return hasNext;
    }

    MutableString previousDescriptionLine;

    /**
     * Obtain the next sequence from the reader over the FASTA formatted content.
     * This method returns true until there is no more sequence to parse in the input. When the method returns false,
     * the content of the parameters descriptionLine and residues is unspecified.
     *
     * @param descriptionLine Where the raw description line will be written.
     * @param residues        When the raw residue lines will be written.
     * @return True if hasNext() is true, False otherwise.
     * @throws IOException
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
     * Try to extract an accession code from a Fasta description line.
     *
     * @param descriptionLine
     */
    public void guessAccessionCode(final CharSequence descriptionLine,
                                   final MutableString accessionCode) {
        accessionCode.setLength(0);
        int startIndex = 0;
        if (descriptionLine.length() > 3 &
                descriptionLine.charAt(0) == 'P' &&
                descriptionLine.charAt(1) == '1' &&
                descriptionLine.charAt(2) == ';'
                ) {
            startIndex = 3;
        }
        for (int i = startIndex; i < descriptionLine.length(); i++) {
            final char c = descriptionLine.charAt(i);
            if (c == ' ' || c == '\t' || c == '|') {
                return;
            }
            accessionCode.append(c);
        }
    }

    /**
     * Filter a string to keep only protein residues.
     *
     * @param rawResidues      A string that may contain any character.
     * @param filteredResidues The subset of characters that represent valid protein residue codes, in the order in which they occur in the rawResidue string.
     */
    public void filterProteinResidues(final CharSequence rawResidues,
                                      final MutableString filteredResidues) {
        filteredResidues.setLength(0);
        for (int i = 0; i < rawResidues.length(); i++) {
            char residueCode = rawResidues.charAt(i);
            if (validProteinResidues.indexOf(residueCode) == -1) {
                continue;
            }
            if (residueCode == '.') {
                residueCode = '-';
            } else {
                residueCode = Character.toUpperCase(residueCode);
            }
            filteredResidues.append(residueCode);
        }

    }

    private boolean readNextDescriptionLine(final FastBufferedReader reader) throws IOException {
        for (; ;) {
            // loop until a line that starts with > if found, or the end of file is reached.
            previousDescriptionLine = reader.readLine(previousDescriptionLine);
            if (previousDescriptionLine == null) {
                return false;
            }
            if (previousDescriptionLine.startsWith(">")) {
                previousDescriptionLine = removeBracket(previousDescriptionLine);
                return true;

            }
        }
    }

    private MutableString removeBracket(final MutableString previousDescriptionLine) {
        return previousDescriptionLine.substring(1, previousDescriptionLine.length());
    }

    private boolean readResidues(final MutableString residues) throws IOException {
        residues.setLength(0);
        line.setLength(0);
        for (; ;) {
            line = reader.readLine(line);
            if (line == null) {
                hasNext = false;
                return hasNext;
            }
            if (line.startsWith(">")) {
                previousDescriptionLine.replace(line);
                previousDescriptionLine = removeBracket(previousDescriptionLine);
                return hasNext;
            }
            residues.append(line);
        }
    }
}
