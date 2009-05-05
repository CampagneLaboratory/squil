/*
 * Copyright (C) 2008-2009 Institute for Computational Biomedicine,
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

import it.unimi.dsi.lang.MutableString;
import junit.framework.TestCase;

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;

/**
 * User: Fabien Campagne
 * Date: Nov 14, 2008
 * Time: 9:53:24 AM
 */
public class TestReaderFastaParser extends TestCase {
    public void testSetReader() throws IOException {
        final String fasta = ">1\n" +
                "ACTG\n" +
                ">2\n" +
                "AAAATT\n" +
                ">3 a\n" +
                "ACCUAT";
        final StringReader r = new StringReader(fasta);
        final ReaderFastaParser p = new ReaderFastaParser(r);
        final MutableString description = new MutableString();
        int index = 1;
        while (p.hasNextSequence()) {
            p.nextSequence(description);
            final Reader baseReader = p.getBaseReader();
            switch (index) {
                case 1:
                    assertEquals(new MutableString("1"), description);
                    assertEquals('A',(char)baseReader.read());
                    assertEquals('C',(char)baseReader.read());
                    assertEquals('T',(char)baseReader.read());
                    assertEquals('G',(char)baseReader.read());
                    assertEquals(-1,baseReader.read());
                    break;

                case 2:
                    assertEquals(new MutableString("2"), description);
                    assertEquals('A', (char)baseReader.read());
                    assertEquals('A', (char)baseReader.read());
                    assertEquals('A', (char)baseReader.read());
                    assertEquals('A', (char)baseReader.read());

                    assertEquals('T', (char)baseReader.read());
                    assertEquals('T', (char)baseReader.read());
                    assertEquals(-1, baseReader.read());
                    break;
                case 3:
                    assertEquals(new MutableString("3 a"), description);
                    assertEquals('A', (char)baseReader.read());
                    assertEquals('C', (char)baseReader.read());
                    assertEquals('C', (char)baseReader.read());
                    assertEquals('U', (char)baseReader.read());
                    assertEquals('A', (char)baseReader.read());
                    assertEquals('T', (char)baseReader.read());
                    assertEquals(-1, baseReader.read());
                    break;
            }

            index++;

        }
        assertEquals("fasta file must be recognized to have 4 sequences.", 4, index);
    }
}
