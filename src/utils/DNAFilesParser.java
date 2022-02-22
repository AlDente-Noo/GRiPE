package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import environment.Cell;
import objects.DNA;

/**
 * fasta file parser
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class DNAFilesParser {
	public static DNA fastaFileParser(String filename){
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

	// returns an error string or an empty string if succeeded
	public static String btrackFileParser(String filename, byte[] closed) {
		ArrayList<Byte> bufferClosed = new ArrayList<>();
		File f = new File(filename);
		if(f.exists()) {
			try {
				BufferedReader reader = new BufferedReader(new FileReader(filename));
				String buffer;

				while ((buffer = reader.readLine()) != null) {
					bufferClosed.add(Byte.parseByte(buffer));
				}
				reader.close();

				if (bufferClosed.size() == closed.length) {
					for (int i = 0; i < bufferClosed.size(); i++) {
						closed[i] = bufferClosed.get(i) == 1 ? Constants.BP_IS_OPEN : Constants.BP_IS_CLOSED;
					}
					return "";
				}
				return "different size of dna strand and dna availability (btrack) files";
			} catch (Exception e) {
				e.printStackTrace();
				return "failed to read the btrack file";
			}
		} else {
			return "the btrack file does not exist";
		}
	}
}
