/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pathseqkaijuconcensussplitter2;

import fastqtools.DelimitedFileReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

/**
 *
 * @author darkosw
 */
public class PathSeqKaijuConcensusSplitter2
{

    /**
     * @param args the command line arguments
     */
    HashMap<String, Double> contigToConfidenceHM;
    //HashMap<String, String[]> contigToTaxonomyHM;
    HashMap<String, HashMap<String, String>> contigToAllTaxonomyHM;
    HashMap<Integer, String> positionToTaxHM;
    HashMap<Integer, String> taxResolutionHM;
    HashMap<Integer, Integer> taxNodeHM;
    HashMap<String, String> taxonomyHM;
    HashMap<String, String> proteinFastaHM;
    ArrayList<String> validTaxAL;
    

    public static void main(String[] args) throws IOException
    {
        PathSeqKaijuConcensusSplitter2 temp = new PathSeqKaijuConcensusSplitter2();
        temp.readProteinFasta(args[2]);
        temp.readNodes(args[4]);
        temp.readTaxonomy(args[3]);
        temp.readKaijuOutput(args[0]);
        temp.readFastaAndFilter(args[1]);
    }

    public PathSeqKaijuConcensusSplitter2()
    {
        contigToConfidenceHM = new HashMap();
        //contigToTaxonomyHM = new HashMap();
        contigToAllTaxonomyHM = new HashMap();
        taxResolutionHM = new HashMap(10000000);
        taxonomyHM = new HashMap(500000000);
        taxNodeHM = new HashMap(500000000);
        proteinFastaHM = new HashMap(1000000);
        validTaxAL = new ArrayList();
        positionToTaxHM = new HashMap();
        positionToTaxHM.put(0, "kingdom");
        positionToTaxHM.put(1, "phylum");
        positionToTaxHM.put(2, "class");
        positionToTaxHM.put(3, "order");
        positionToTaxHM.put(4, "family");
        positionToTaxHM.put(5, "genus");
        positionToTaxHM.put(6, "species");
        
        validTaxAL.add("kingdom");
        validTaxAL.add("phylum");
        validTaxAL.add("class");
        validTaxAL.add("order");
        validTaxAL.add("family");
        validTaxAL.add("genus");
        validTaxAL.add("species");

    }

    public void readProteinFasta(String s) throws FileNotFoundException
    {
        Scanner sc = new Scanner(new File(s));
        String header, sequence;
        String[] tempStringArray;
        while (sc.hasNextLine())
        {
            header = sc.nextLine();
            tempStringArray = header.split(" ");
            header = tempStringArray[0];
            header = header.substring(1);
            sequence = sc.nextLine();
            if (sequence.endsWith("*"))
            {
                sequence = sequence.replaceAll("\\*", "");
            }
            proteinFastaHM.put(header, sequence);
        }
        sc.close();
        System.err.println(proteinFastaHM.size() + " predicted proteins stored");
    }

    public void readFastaAndFilter(String s) throws FileNotFoundException
    {
        HashMap<String, PrintWriter> fastaPrintWriterHM = new HashMap();
        HashMap<String, PrintWriter> tablePrintWriterHM = new HashMap();
        PrintWriter tempFastaPrintWriter, tempTablePrintWriter;
        HashMap<String, String> tempHM;
        Iterator<String> itr;
        String header, contigID, sequence, currentTax, taxAssignment;
        String[] tempStringArray;
        Scanner sc = new Scanner(new File(s));
        while (sc.hasNextLine())
        {
            header = sc.nextLine();
            sequence = sc.nextLine();
            tempStringArray = header.split(" ");
            contigID = tempStringArray[0].substring(1);
            
            if (contigToAllTaxonomyHM.containsKey(contigID))
            {

                tempHM = contigToAllTaxonomyHM.get(contigID);
                itr = tempHM.keySet().iterator();
                while (itr.hasNext())
                {
                    currentTax = itr.next();
                    if (tempHM.containsKey(currentTax))
                    {
                        taxAssignment = tempHM.get(currentTax).replaceAll(" ", "-");
                        if (fastaPrintWriterHM.containsKey(currentTax))
                        {
                            tempFastaPrintWriter = fastaPrintWriterHM.get(currentTax);
                            tempTablePrintWriter = tablePrintWriterHM.get(currentTax);
                        } else
                        {
                            tempFastaPrintWriter = new PrintWriter(new File(currentTax + "_sequences.fa"));
                            tempTablePrintWriter = new PrintWriter(new File(currentTax + "_table.csv"));
                            tempTablePrintWriter.println("contigID,taxAssignment");
                        }
                        tempFastaPrintWriter.println(header);
                        tempFastaPrintWriter.println(sequence);

                        tempTablePrintWriter.println(contigID + "," + taxAssignment);

                        fastaPrintWriterHM.put(currentTax, tempFastaPrintWriter);
                        tablePrintWriterHM.put(currentTax, tempTablePrintWriter);
                    }
                }
            } else
            {
                currentTax = "unclassified";
                if (fastaPrintWriterHM.containsKey(currentTax))
                {
                    tempFastaPrintWriter = fastaPrintWriterHM.get(currentTax);
                } else
                {
                    tempFastaPrintWriter = new PrintWriter(new File(currentTax + "_sequences.fa"));
                }
                tempFastaPrintWriter.println(header);
                tempFastaPrintWriter.println(sequence);
                fastaPrintWriterHM.put(currentTax, tempFastaPrintWriter);
            }
        }
        itr = fastaPrintWriterHM.keySet().iterator();
        while (itr.hasNext())
        {
            currentTax = itr.next();
            tempFastaPrintWriter = fastaPrintWriterHM.get(currentTax);
            tempFastaPrintWriter.close();
            if (tablePrintWriterHM.containsKey(currentTax))
            {
                tempTablePrintWriter = tablePrintWriterHM.get(currentTax);
                tempTablePrintWriter.close();
            }
        }
    }
     
    public void readKaijuOutput(String kaijuString) throws IOException
    {
        HashMap<String, String> tempHM;
        HashMap<String, Integer> resolutionCounterHM = new HashMap();
        DelimitedFileReader dfr = new DelimitedFileReader("\t");
        dfr.readFromFile(new File(kaijuString));
        String classification, taxLevel, contigID, accessionID, taxonomy, proteinAlignment;
        Double alignmentCoverage;
        Integer taxNumber, parentTaxNumber;
        String[] tempStringArray1, tempStringArray2;
        while (dfr.parseDelimitedLine())
        {
            tempStringArray1 = dfr.getDelimitedLine();
            classification = tempStringArray1[0];
            if (classification.equalsIgnoreCase("C"))
            {
                contigID = tempStringArray1[1];
                System.err.println(contigID);
                taxNumber = new Integer(tempStringArray1[2]);
                System.err.println(taxNumber);
                taxLevel = taxResolutionHM.get(taxNumber).replaceAll(" ", "_");
                System.err.println(taxLevel);
                // instead of synonyms, let's go up the node tree to find taxLevel
                //taxLevel = taxSynonymsHM.get(taxLevel);
                //System.err.println(taxLevel);
                
                if (taxNumber.compareTo(1) == 0) //there's no higher rank
                {
                    taxLevel = "unclassified";
                }
                else if(!validTaxAL.contains(taxLevel))
                {
                    System.err.println("Not a valid taxLevel");
                    do
                    {
                        parentTaxNumber = taxNodeHM.get(taxNumber);
                        taxLevel = taxResolutionHM.get(taxNumber).replaceAll(" ", "_");
                        System.err.println(taxLevel);
                        taxNumber = parentTaxNumber;
                    }
                    while(!validTaxAL.contains(taxLevel) && !taxLevel.equalsIgnoreCase("superkingdom") && taxNumber.compareTo(1) != 0);
                    if(taxLevel.equalsIgnoreCase("superkingdom") || taxNumber.compareTo(1) == 0)
                    {
                        taxLevel = "unclassified";
                    }
                }
                else
                {
                    System.err.println("Is a valid taxLevel");
                }
                System.err.println("Settled on: "+taxLevel);

                tempStringArray2 = tempStringArray1[5].split(",");
                accessionID = tempStringArray2[0];

                tempStringArray2 = tempStringArray1[6].split(",");
                proteinAlignment = tempStringArray2[0];

                alignmentCoverage = (double) proteinAlignment.length() / (double) proteinFastaHM.get(contigID).length();
                System.err.println(alignmentCoverage);
                if (alignmentCoverage.compareTo(0.70) >= 0 && !taxLevel.equalsIgnoreCase("unclassified"))
                {
                    System.err.println(alignmentCoverage+" is greater than 0.70");
                    if (taxonomyHM.containsKey(accessionID))
                    {
                        taxonomy = taxonomyHM.get(accessionID);
                        tempHM = makeTaxonomyHM(taxonomy, taxLevel);
                        contigToAllTaxonomyHM.put(contigID, tempHM);
                    } else
                    {
                        taxLevel = "unclassified";
                        System.err.println("Unclassified because taxonomy wasn't found!");
                    }
                }
                else
                {
                    taxLevel = "unclassified";
                }

                //System.err.println(accessionID);
                //System.err.println(taxNumber);
                //System.err.println("Contig: " + contigID);
                //System.err.println("Taxonomy: " + taxonomyHM.get(accessionID));
                //System.err.println("Resolution: " + taxLevel);
                if (resolutionCounterHM.containsKey(taxLevel))
                {
                    resolutionCounterHM.put(taxLevel, resolutionCounterHM.get(taxLevel) + 1);
                } else
                {
                    resolutionCounterHM.put(taxLevel, 1);
                }

            } else
            {
                taxLevel = "unclassified";
                if (resolutionCounterHM.containsKey(taxLevel))
                {
                    resolutionCounterHM.put(taxLevel, resolutionCounterHM.get(taxLevel) + 1);
                } else
                {
                    resolutionCounterHM.put(taxLevel, 1);
                }
            }
        }
        Iterator<String> itr = resolutionCounterHM.keySet().iterator();
        while (itr.hasNext())
        {
            taxLevel = itr.next();
            System.err.println(taxLevel + " : " + resolutionCounterHM.get(taxLevel));
        }
        System.err.println("All contigs with taxonomy: " + contigToAllTaxonomyHM.size());
    }

    public HashMap<String, String> makeTaxonomyHM(String tempString, String taxLevel)
    {
        HashMap<String, String> toReturn = new HashMap();
        String currentTax, taxAssignment;
        String[] tempStringArray = tempString.split("; ");

        for (int i = 0; i < tempStringArray.length; i++)
        {
            currentTax = positionToTaxHM.get(i);
            //System.err.println(tempStringArray[i]);
            //System.err.println(tempStringArray[i].substring(3));
            taxAssignment = tempStringArray[i];
            int j = i - 1;
            while (j >= 0)
            {
                taxAssignment = tempStringArray[j] + "; " + taxAssignment;
                j--;
            }
            if (taxAssignment.length() > 0)
            {
                toReturn.put(currentTax, taxAssignment);
            }
            if (taxLevel.equalsIgnoreCase(currentTax))
            {
                i = tempStringArray.length;
            }

        }
        return toReturn;
    }

    public void readNodes(String nodeString) throws IOException
    {
        DelimitedFileReader dfr = new DelimitedFileReader("\t");
        dfr.readFromFile(new File(nodeString));
        Integer taxNumber, parentTaxNumber;
        String taxLevel;
        String[] tempStringArray;
        while (dfr.parseDelimitedLine())
        {
            tempStringArray = dfr.getDelimitedLine();
            taxNumber = new Integer(tempStringArray[0]);
            parentTaxNumber = new Integer(tempStringArray[2]);
            taxLevel = tempStringArray[4];
            //System.err.println(taxNumber + " : " + taxLevel);
            taxResolutionHM.put(taxNumber, taxLevel);
            taxNodeHM.put(taxNumber, parentTaxNumber);
        }
        System.err.println(taxResolutionHM.size());
        System.err.println(taxNodeHM.size());
    }

    public void readTaxonomy(String taxonomyString) throws IOException
    {
        DelimitedFileReader dfr = new DelimitedFileReader("\t");
        dfr.readFromFile(new File(taxonomyString));
        String accessionID, taxonomy;
        String[] tempStringArray;
        while (dfr.parseDelimitedLine())
        {
            tempStringArray = dfr.getDelimitedLine();
            accessionID = tempStringArray[0];
            taxonomy = tempStringArray[1];
            taxonomyHM.put(accessionID, taxonomy);
        }
        System.err.println(taxonomyHM.size());
    }

}
