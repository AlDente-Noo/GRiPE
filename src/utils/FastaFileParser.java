package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import objects.DNA;

/**
 * fasta file parser
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class FastaFileParser {
	public static DNA fileParser(String filename){
		DNA dna=new DNA();
		File f = new File(filename);
		if(f.exists()){
			ArrayList<Byte> strand = new ArrayList<Byte>();
			BufferedReader reader = null;
		        try
		        {
		            reader = new BufferedReader(new FileReader(filename));
		            String text = null;
		            
		            String currentName = "";
		            while ((text = reader.readLine()) != null){
			            	text=text.trim();
			            	if(!text.isEmpty()){
			            		if(text.startsWith(">")){	
			            			currentName =  text.replaceAll(">", "").trim();
			            		} else{	
			            			strand.addAll(CellUtils.getSequenceIDs(text));
			            		}
			            	}
		            }
		            
		            dna.description = currentName;
		            dna.loadSequence(strand);
		        } catch (Exception e) {
		            e.printStackTrace();
		        } finally {
		            try
		            {
		                if (reader != null){
		                    reader.close();
		                }
		            } catch (IOException e)
		            {
		                e.printStackTrace();
		            }
		        }
		}
		
	    return dna;
	}
}
