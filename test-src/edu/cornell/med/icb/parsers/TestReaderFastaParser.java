package edu.cornell.med.icb.parsers;

import it.unimi.dsi.lang.MutableString;

import java.io.StringReader;
import java.io.IOException;
import java.io.Reader;

import junit.framework.TestCase;

/**
 * User: Fabien Campagne
 * Date: Nov 14, 2008
 * Time: 9:53:24 AM
 */
public class TestReaderFastaParser extends TestCase {

    public void testSetReader() throws IOException {

        String fasta = ">1\n" +
                "ACTG\n" +
                ">2\n" +
                "AAAATT\n" +
                ">3 a\n" +
                "ACCUAT";
        StringReader r = new StringReader(fasta);
        ReaderFastaParser p = new ReaderFastaParser(r);
        MutableString description = new MutableString();
        int index = 1;
        while (p.hasNextSequence()) {

            p.nextSequence(description);
            Reader baseReader = p.getBaseReader();
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
                    assertEquals('A',(char)baseReader.read());
                    assertEquals('A',(char)baseReader.read());
                    assertEquals('A',(char)baseReader.read());
                    assertEquals('A',(char)baseReader.read());

                    assertEquals('T',(char)baseReader.read());
                    assertEquals('T',(char)baseReader.read());
                    assertEquals(-1,baseReader.read());
                      break;
                case 3:
                    assertEquals(new MutableString("3 a"), description);
                     assertEquals('A',(char)baseReader.read());
                     assertEquals('C',(char)baseReader.read());
                     assertEquals('C',(char)baseReader.read());
                     assertEquals('U',(char)baseReader.read());
                     assertEquals('A',(char)baseReader.read());
                     assertEquals('T',(char)baseReader.read());
                    assertEquals(-1,baseReader.read());
                      break;
            }

          index++;

        }
        assertEquals("fasta file must be recognized to have 4 sequences.", 4, index);
    }


}
