package edu.cornell.med.icb.plink;

import edu.mssm.crover.cli.CLI;
import edu.cornell.med.icb.io.TSVReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;

import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.FileWriter;
import java.util.List;

import org.apache.commons.lang.math.NumberUtils;
/*
 * Copyright (C) 2001-2002 Mount Sinai School of Medicine
 * Copyright (C) 2003-2008 Institute for Computational Biomedicine,
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

/**
 * @author: Fabien Campagne
 * Date: Sep 12, 2008
 * Time: 11:43:18 AM
 */
public class ConvertDbSNP2PlinkMap {
    public static void main(String args[]) {
        ConvertDbSNP2PlinkMap adapter = new ConvertDbSNP2PlinkMap();

        adapter.process(args);
    }

    ObjectSet<MutableString> allSnps;

    private void process(String[] args) {
        String dbSNPFilename = CLI.getOption(args, "-i", null); //
        if (dbSNPFilename == null) {
            System.out.println("No input specified (-i). Provide the concatenation of the chromosome specific files tab delimited files, see ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/chr_rpts/.");
        }
        String snpListFilename = CLI.getOption(args, "--snp-list", "all.snps");
        String mapFilename = CLI.getOption(args, "--map", "all.map");

        LineIterator it = null;
        try {
            it = new LineIterator(new FastBufferedReader(new FileReader(snpListFilename)));
            int mapCount = 0;

            FileWriter mapOut = new FileWriter(mapFilename);

            List<MutableString> allSnpsList = it.allLines();
            allSnps = new ObjectOpenHashSet<MutableString>();
            allSnps.addAll(allSnpsList);
            System.out.println(String.format("Read %d SNPs from %s", allSnps.size(), snpListFilename));
            TSVReader reader = new TSVReader(new FastBufferedReader(new FileReader(dbSNPFilename)));
            //    skipLines(reader, 7);
            while (reader.hasNext()) {
            
                reader.next();
                if (reader.numTokens() < 10) continue;
                String id = reader.getString();
                if ("rs#".equals(id) || id.length()==0) continue;
                int status = reader.getInt();
                if (status != 2) continue;
                skipFields(reader, 4);
                String chromosome = reader.getString(); //field 7
                skipFields(reader, 4);
                String positionString = reader.getString();      // field 12
                if (!NumberUtils.isNumber(positionString)) continue;
                int position = Integer.valueOf(positionString);
                skipFields(reader, 9);
                String mappedTo = reader.getString(); // field 22
                if (!"reference".equals(mappedTo)) continue;

                MutableString snpId = new MutableString("rs");
                snpId.append(id);
                snpId.compact();
                if (allSnps.contains(snpId)) {
                    mapOut.append(String.format("%s %s 0 %d\n", chromosome, snpId, position));
                    mapCount++;
                }

            }
            mapOut.close();
            System.out.println(String.format("Mapped %d SNPs", mapCount));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.err.println("Cannot open snp-list file (--snp-list option).");
            System.exit(10);
        } catch (IOException e) {
            e.printStackTrace();
            System.err.println("Error reading dbSNP file (-i option) " + dbSNPFilename);
            System.exit(10);
        }


    }

    private void skipFields(TSVReader reader, int i) {
        for (; i > 0; i--) {
            reader.getString();
        }
    }

    private void skipLines(TSVReader reader, int i) throws IOException {
        for (; reader.hasNext() && i > 0; i--) {

            reader.skip();

        }
    }
}

