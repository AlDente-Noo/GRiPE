package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.ArrayList;

/**
 * this class parses a CSV file and returns a matrix with the sored data and a vector with the headers
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class CSVparser {
	public ArrayList<ArrayList<String>> data;
	public HashMap<String,Integer> header;
	public char textDelimiter;
	public char cellDelimiter;
	public boolean loaded;
	public boolean foundDelimiter;
	public boolean hasVariousRowSizes;
	
	/**
	 * class constructor
	 * @param filename the file
	 * @param hasHeader  true if it has header 
	 */
	public CSVparser(String filename, boolean hasHeader, boolean processHeader){
		loaded = loadCSV(filename, hasHeader,processHeader, Constants.CSV_FILE_TEXT_DELIMITER);
	}
	
	/**
	 * class constructor
	 * @param filename the CSV file
	 * @param hasHeader true if it has header 
	 * @param textDelimiter how text is delimited, usually "
	 */
	public CSVparser(String filename, boolean hasHeader, boolean processHeader, char textDelimiter){
		loaded = loadCSV(filename, hasHeader,processHeader, textDelimiter);
	}
	
	/**
	 * loads a csv file
	 * @param filename the csv file
	 * @param hasHeader true if it has header and false otherwise
	 * @param textDelimiter the text delimiter
	 */
	private boolean loadCSV(String filename, boolean hasHeader, boolean processHeader, char textDelimiter){
		
		//if file doesn't exists skip this
		File f = new File(filename);
		if(!f.exists()){
			return false;
		}
		this.textDelimiter = textDelimiter;
		this.data = new ArrayList<ArrayList<String>>();
		this.header = new HashMap<String, Integer>();
		
		ArrayList<String> lines = new ArrayList<String>();
		BufferedReader reader = null;

		//get all lines
	    try{
            reader = new BufferedReader(new FileReader(filename));
            String text = null;
            while ((text = reader.readLine()) != null){
	            text=text.trim();
	            if(!text.isEmpty()){
	            	lines.add(text);
	            }
            }
	            
       	} catch (Exception e) {
            e.printStackTrace();
        } finally {
            try{
            	if (reader != null){
                    reader.close();
                }
            } catch (IOException e){
                e.printStackTrace();
            }
        }
        
        boolean result = false;
        this.hasVariousRowSizes = false;
        
        if(lines.size()>0){
        	
	        //get cell delimiter
        	int index = 0;
	        while(index < lines.size()){
	        	getCellDelimiter(lines.get(index));
	        	if(this.foundDelimiter){
	        		break;
	        	} else{
	        		index++;
	        	}
	        }
	        
	        // found cell delimiter then parse
	        if(this.foundDelimiter){
	        	ArrayList<String> buffer;
	            int minLineSize=-1;
	        	if(hasHeader){
	        		buffer = processLine(lines.get(0), this.cellDelimiter, this.textDelimiter);
	        		if(buffer!=null && buffer.size()>0){
	        			for(int i=0;i<buffer.size();i++){
	        				if(processHeader){
	        					buffer.set(i,this.processHeader(buffer.get(i)));
	        				}
	        				this.header.put(buffer.get(i), i);
	        			}
	        			if(minLineSize<0){
	        				minLineSize = buffer.size();
	        			} 
	        			
	        		}
	        		index++;
	        	}
	        	
	        	while(index<lines.size()){
	        		buffer = processLine(lines.get(index), this.cellDelimiter, this.textDelimiter);
	        		if(!buffer.isEmpty()){
		        		this.data.add(buffer);
		        		
	        			if(minLineSize<0){
	        				minLineSize = buffer.size();
	        			} else{
	        				if(minLineSize!=buffer.size()){
	        					this.hasVariousRowSizes=true;
	        					if(minLineSize > buffer.size()){
	        						minLineSize = buffer.size();
	        					}
	        				}
	        			}
	        		}
	        		
	        		index++;
	        	}
	        	
	        	result = true;
	        }
        }
        
        
	
        
        return result;
        
	}
	
	/**
	 * checks to find delimiter on line
	 * @param text
	 */
	private void getCellDelimiter(String text){
		char currentChar;
		boolean inCell=false;
		this.foundDelimiter = false;
		int currentPos = 0;
		while(currentPos<text.length() && !this.foundDelimiter){
			if(text.charAt(currentPos) == this.textDelimiter){
				inCell = !inCell;
			} else if(!inCell){
				currentChar = text.charAt(currentPos);
				for(char c:Constants.CSV_FILE_CELL_DELIMTER){
					if(currentChar == c){
						this.cellDelimiter = c;
						this.foundDelimiter = true;
						break;
					}
				}
			}
			currentPos++;	
		}
		
	}
	
	/**
	 * breaks a line into elements
	 * @param text the line to be parsed
	 * @param cellDelimiter the char that separates the encapsulates text cells, usually "
	 * @param textDelimiter the char that separates the cell, usually , or ;
	 * @return an array list with the cells
	 */
	private ArrayList<String> processLine(String text, char cellDelimiter, char textDelimiter){
		ArrayList<String> result = new ArrayList<String>();
		boolean inCell=false;
		int currentPos = 0, lastPos=0;
		while(currentPos<text.length()){
			if(text.charAt(currentPos) == textDelimiter){
				inCell = !inCell;
			} else if(text.charAt(currentPos) == cellDelimiter && !inCell){
				result.add(getCleanCell(text, lastPos, currentPos,textDelimiter));
				lastPos = currentPos+1;
			}
			currentPos++;	
		}
		
		//add last cell
		if(lastPos < currentPos){
			result.add(getCleanCell(text, lastPos, currentPos,textDelimiter));
		}
		
		return result;
	}
	
	/**
	 * cleans a text cell and returns the string
	 * @param text the entire line 
	 * @param start where the cell starts 
	 * @param end where the cell ends
	 * @param textDelimiter the text delimiter used
	 * @return
	 */
	private String getCleanCell(String text, int start, int end, char textDelimiter){
		String result="";
		if(end > start){
			result = text.substring(start, end).replace(textDelimiter, ' ').trim();
		}
		return result;
	}
	
	/**
	 * just uppercase letters
	 * @param text
	 * @return
	 */
	private String processHeader(String text){
		return text.toUpperCase().replaceAll("[^a-zA-Z]", "");
	}
	
}
