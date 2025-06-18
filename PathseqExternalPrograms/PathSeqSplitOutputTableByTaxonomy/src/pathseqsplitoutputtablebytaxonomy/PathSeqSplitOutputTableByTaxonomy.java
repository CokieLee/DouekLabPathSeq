/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pathseqsplitoutputtablebytaxonomy;

import fastqtools.DelimitedFileReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author darkosw
 */
public class PathSeqSplitOutputTableByTaxonomy
{

    /**
     * @param args the command line arguments
     */
    String taxonomyLevel;
    public static void main(String[] args) throws IOException
    {
        PathSeqSplitOutputTableByTaxonomy temp = new PathSeqSplitOutputTableByTaxonomy(args[1]);
        temp.readAndWrite(args[0]);
    }
    public PathSeqSplitOutputTableByTaxonomy(String taxonomyLevel)
    {
        this.taxonomyLevel = taxonomyLevel;
    }
    
    public void readAndWrite(String s) throws IOException
    {
        DelimitedFileReader dfr = new DelimitedFileReader(",");
        dfr.readFromFile(new File(s));
        String header = dfr.getNextLine();
        //System.err.println(header);
        HashMap<String, PrintWriter> hm = new HashMap();
        PrintWriter tempPW;
        ArrayList<String> splitAL = new ArrayList();
        String[] tempArray, tempTaxArray;
        String taxLevelAssignment, stringToTest, toPrint;
        while(dfr.parseDelimitedLine())
        {
            tempArray = dfr.getDelimitedLine();
            tempTaxArray = tempArray[0].split(";");
            for(int i=0; i<tempTaxArray.length; i++)
            {
                stringToTest = tempTaxArray[i];
                if(stringToTest.contains(taxonomyLevel+"__"))
                {
                    taxLevelAssignment = stringToTest.replaceAll(taxonomyLevel+"__", "");
                    taxLevelAssignment = taxLevelAssignment.replaceFirst("-", "");
                    //System.err.println(taxLevelAssignment);
                    i = tempTaxArray.length;
                    toPrint = tempArray[0];
                    if(hm.containsKey(taxLevelAssignment))
                    {
                        tempPW = hm.get(taxLevelAssignment);
                    }
                    else
                    {
                        tempPW = new PrintWriter(taxLevelAssignment+"_"+s);
                        tempPW.println(header);
                        splitAL.add(taxLevelAssignment);
                    }
                    for(int j=1; j<tempArray.length; j++)
                    {
                        toPrint = toPrint+","+tempArray[j];
                    }
                    toPrint = toPrint.replaceAll(";-", ";");
                    tempPW.println(toPrint);
                    hm.put(taxLevelAssignment, tempPW);
                }
            }
        }
        for(int i=0; i<splitAL.size(); i++)
        {
            hm.get(splitAL.get(i)).close();
        }
    }
}
