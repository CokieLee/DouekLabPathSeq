/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pathseqmergeqiime2taxandsalmon;

import fastqtools.DelimitedFileReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author darkosw
 */
public class PathSeqMergeQIIME2TaxAndSalmon
{

    HashMap<String, String> contigToTaxonomyHM;
    HashMap<String, HashMap<String, Double>> sampleToTaxonomyTPMHM;
    HashMap<String, HashMap<String, Double>> sampleToTaxonomyPseudocountHM;
    ArrayList<String> sampleAL;
    ArrayList<String> taxonomyAL;
    String prefix;

    public static void main(String[] args) throws IOException
    {
        String[] fileStringArray = args[0].split(",");
        String[] sampleStringArray = args[1].split(",");
        PathSeqMergeQIIME2TaxAndSalmon temp = new PathSeqMergeQIIME2TaxAndSalmon(args[3]);
        temp.readQIIME2Taxonomy(args[2]);

        for (int i = 0; i < fileStringArray.length; i++)
        {
            temp.readCounts(fileStringArray[i], sampleStringArray[i]);
        }
        temp.printCounts();
    }

    public PathSeqMergeQIIME2TaxAndSalmon(String prefix)
    {
        contigToTaxonomyHM = new HashMap(20000);
        sampleToTaxonomyTPMHM = new HashMap(1000);
        sampleToTaxonomyPseudocountHM = new HashMap(1000);
        sampleAL = new ArrayList(1000);
        this.prefix = prefix;
    }

    public void printCounts() throws FileNotFoundException
    {
        HashMap<String, Double> tpmTempHM = new HashMap(10000);
        HashMap<String, Double> pseudocountTempHM = new HashMap(10000);

        String toPrint = "taxonomy";
        String toPrintTPM, toPrintPseudocount;
        String taxonomyAssignment, sampleID;
        Double tpm, pseudocount;

        PrintWriter tpmPW = new PrintWriter(new File(prefix + "_tpm.csv"));
        PrintWriter pseudocountPW = new PrintWriter(new File(prefix + "_pseudocounts.csv"));

        for (int i = 0; i < sampleAL.size(); i++)
        {
            toPrint = toPrint + "," + sampleAL.get(i);
        }
        tpmPW.println(toPrint);
        pseudocountPW.println(toPrint);

        for (int i = 0; i < taxonomyAL.size(); i++)
        {
            taxonomyAssignment = taxonomyAL.get(i);
            toPrintTPM = taxonomyAssignment;
            toPrintPseudocount = taxonomyAssignment;
            {
                if (!toPrint.endsWith("__"))
                {
                    for (int j = 0; j < sampleAL.size(); j++)
                    {
                        sampleID = sampleAL.get(j);
                        tpmTempHM = sampleToTaxonomyTPMHM.get(sampleID);
                        pseudocountTempHM = sampleToTaxonomyPseudocountHM.get(sampleID);
                        tpm = tpmTempHM.get(taxonomyAssignment);
                        pseudocount = pseudocountTempHM.get(taxonomyAssignment);
                        toPrintTPM = toPrintTPM + "," + tpm;
                        toPrintPseudocount = toPrintPseudocount + "," + pseudocount;
                    }
                    tpmPW.println(toPrintTPM);
                    pseudocountPW.println(toPrintPseudocount);
                }
            }

        }
        tpmPW.close();
        pseudocountPW.close();
    }

    public void readCounts(String file, String sampleID) throws IOException
    {
        HashMap<String, Double> tpmTempHM = new HashMap(10000);
        HashMap<String, Double> pseudocountTempHM = new HashMap(10000);
        DelimitedFileReader dfr = new DelimitedFileReader("\t");
        dfr.readFromFile(new File(file));
        dfr.getNextLine(); //skip header
        sampleAL.add(sampleID);
        String[] tempStringArray;
        Double tpm, pseudocount;
        String contigID, taxonomyAssignment;
        while (dfr.parseDelimitedLine())
        {
            tempStringArray = dfr.getDelimitedLine();
            contigID = tempStringArray[0];
            tpm = new Double(tempStringArray[3]);
            pseudocount = new Double(tempStringArray[4]);
            taxonomyAssignment = contigToTaxonomyHM.get(contigID);
            if (tpmTempHM.containsKey(taxonomyAssignment))
            {
                tpm = tpm + tpmTempHM.get(taxonomyAssignment);
                pseudocount = pseudocount + pseudocountTempHM.get(taxonomyAssignment);
            }

            tpmTempHM.put(taxonomyAssignment, tpm);
            pseudocountTempHM.put(taxonomyAssignment, pseudocount);
        }
        System.err.println(tpmTempHM.size());
        sampleToTaxonomyTPMHM.put(sampleID, tpmTempHM);
        sampleToTaxonomyPseudocountHM.put(sampleID, pseudocountTempHM);
    }

    public void readQIIME2Taxonomy(String s) throws IOException
    {
        Set<String> taxonomySet = new HashSet<>(60000); //start large for easy resizing
        DelimitedFileReader dfr = new DelimitedFileReader(",");
        dfr.readFromFile(new File(s));
        dfr.getNextLine(); //skip header
        String[] tempStringArray;
        String contigID, taxonomyAssignment;
        while (dfr.parseDelimitedLine())
        {
            tempStringArray = dfr.getDelimitedLine();
            contigID = tempStringArray[0];
            taxonomyAssignment = tempStringArray[1];
            taxonomySet.add(taxonomyAssignment);
            contigToTaxonomyHM.put(contigID, taxonomyAssignment);
        }
        taxonomyAL = new ArrayList(taxonomySet);
        System.err.println(contigToTaxonomyHM.size() + " total contigs");
        System.err.println(taxonomyAL.size() + " total taxonomies");
    }
}
