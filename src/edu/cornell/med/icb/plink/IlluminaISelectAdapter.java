package edu.cornell.med.icb.plink;

import edu.mssm.crover.cli.CLI;
import edu.cornell.med.icb.io.TSVReader;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

import java.io.*;
import java.util.List;

import it.unimi.dsi.mg4j.util.MutableString;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import it.unimi.dsi.fastutil.io.TextIO;
import it.unimi.dsi.io.FastBufferedReader;
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
 * Converts Illumina genotyping reports to plink format. We assume case control studies where
 * family information is not available and code family for each sample as 1.
 * Generates
 * <LI>*.lgen LGEN format for each input file.
 * <LI>*.snps list of unique snps in each each input file.
 * <LI>*.samples list of unique samples in each each input file.
 * <LI>all.snps list of unique snps across all the input files.
 * <LI>all-phenotypes.cov associates each sample to a phenotype, according to which input file the sample appears into.
 * The .cov file can be used with plink <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml">--make-pheno</a> option to recode phenotypes *
 *
 * To make a basic FAM file from the results for association studies, do  awk '{print "1 "$1" 1 1 other 0"}' Control.samples > Control.fam

 *
 * @author: Fabien Campagne Date: Sep 11, 2008 Time: 6:22:00 PM
 */
public class IlluminaISelectAdapter {
    private ObjectSet<MutableString> allSnps;

    public static void main(String args[]) {
        IlluminaISelectAdapter adapter = new IlluminaISelectAdapter();

        adapter.process(args);
    }

    private void process(final String[] args) {
        allSnps = new ObjectOpenHashSet<MutableString>();
        String reportFilenames[] = CLI.getOptions(args, "-i");
        if (reportFilenames.length == 0) {
            System.out.println("No input specified (-i)");
        }
        String basenames[] = CLI.getOptions(args, "-b");
        if (basenames.length == 0) basenames = null;
        int fileIndex = 0;
        String snpsFilename = "all.snps";
        String phenotypeFilename = "all-phenotypes.cov";

        FileWriter snpsOut = null;
        FileWriter phenotypeOut = null;
        try {
            snpsOut = new FileWriter(snpsFilename);
            phenotypeOut = new FileWriter(phenotypeFilename);
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (String reportFilename : reportFilenames) {
            try {
                String basename = FilenameUtils.getBaseName(reportFilename);
                if (basenames != null) {
                    basename = basenames[fileIndex];
                }
                System.out.println("processing " + reportFilename + " basename: " + basename);
                convertReport(reportFilename, basename, phenotypeOut);
                fileIndex++;
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        try {
            dumpList(snpsOut, allSnps);
            phenotypeOut.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void convertReport(final String reportFilename, String basename, FileWriter phenotypeOut) throws IOException {

        String lgenFilename = basename + ".lgen";
        FileWriter lgenOut = new FileWriter(lgenFilename);
        String snpsFilename = basename + ".snps";
        FileWriter snpsOut = new FileWriter(snpsFilename);
        String samplesFilename = basename + ".samples";
        FileWriter samplesOut = new FileWriter(samplesFilename);

        ObjectSet<MutableString> snps = new ObjectOpenHashSet<MutableString>();
        ObjectSet<MutableString> samples = new ObjectOpenHashSet<MutableString>();
        TSVReader reader = new TSVReader(new FastBufferedReader(new FileReader(reportFilename)));
        while (reader.hasNext()) {
            reader.next();
            String token = reader.getString();
            if ("SNP Name".equals(token)) break;
        }
        while (reader.hasNext()) {
            reader.next();
            String snpId = reader.getString();
            snps.add(new MutableString(snpId));
            String sampleId = reader.getString();
            samples.add(new MutableString(sampleId));
            String allele1 = reader.getString();
            allele1 = allele1.replace('-', '0');
            String allele2 = reader.getString();
            allele2 = allele2.replace('-', '0');
            lgenOut.append(String.format("1 %s %s %s %s\n", sampleId, snpId, allele1, allele2));
        }
        lgenOut.close();
        dumpList(samplesOut, samples);
        dumpList(snpsOut, snps);
        allSnps.addAll(snps);
        for (MutableString sampleId : samples) {

            phenotypeOut.append(String.format("1 %s %s\n", sampleId, basename));
        }
        phenotypeOut.flush();
        samplesOut.close();
    }

    private void dumpList(FileWriter samplesOut, ObjectSet<MutableString> samples) throws IOException {
        for (MutableString sampleId : samples) {

            samplesOut.append(sampleId);
            samplesOut.append('\n');
        }
        samplesOut.close();
    }


}
