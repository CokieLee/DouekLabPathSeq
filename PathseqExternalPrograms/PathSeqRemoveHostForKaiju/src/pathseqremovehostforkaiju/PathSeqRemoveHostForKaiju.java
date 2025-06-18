/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pathseqremovehostforkaiju;

import fastqtools.BlastReader;
import fastqtools.BlastRecord;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Scanner;

/**
 *
 * @author darkosw
 */
public class PathSeqRemoveHostForKaiju
{

    HashMap<String, String> alignedHM;
    public static void main(String[] args) throws IOException
    {
        PathSeqRemoveHostForKaiju temp = new PathSeqRemoveHostForKaiju();
        temp.readBlastOutput(args[0]);
        temp.readAndWriteUnalignedFasta(args[1]);
    }
    public PathSeqRemoveHostForKaiju()
    {
        alignedHM = new HashMap();
    }
    public void readAndWriteUnalignedFasta(String s) throws FileNotFoundException
    {
        Scanner sc = new Scanner(new File(s));
        String fastaHeader, fastaSequence, blastQueryName;
        String[] tempStringArray;
        int aligned=0;
        int unaligned = 0;
        while(sc.hasNextLine())
        {
            fastaHeader = sc.nextLine();
            fastaSequence = sc.nextLine();
            tempStringArray = fastaHeader.split(" ");
            blastQueryName = tempStringArray[0].substring(1);
            System.err.println(blastQueryName);
            if(!alignedHM.containsKey(blastQueryName))
            {
                System.out.println(fastaHeader);
                System.out.println(fastaSequence);
                unaligned++;
            }
            else
            {
                aligned++;
            }
        }
        sc.close();
        System.err.println("aligned: "+aligned);
        System.err.println("unaligned: "+unaligned);
    }
    public void readBlastOutput(String s) throws IOException
    {
        BlastReader br = new BlastReader();
        br.readFromFile(new File(s));
        BlastRecord tempBlastRecord;
        while(br.parseBlastRecord())
        {
            tempBlastRecord = br.getBlastRecord();
            alignedHM.put(tempBlastRecord.getQueryID(), "aligned");
        }
    }
    
}
