package objects;

import java.io.Serializable;

import utils.Constants;
import utils.Utils;

/**
 * class that specifies a DNA region
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class DNAregion  implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -1693323042845061156L;
	public long start; //position of the site on the DNA
	public long end; //position of the site on the DNA
	public String chromosome;
	public double scaleFactor;
	public int direction;	
	public double probability;
	public boolean hasProbability;
	public boolean hasDirection;
	
	/**
	 * class constructor
	 * @param chromosome the name of the chromosome where the region resides
	 * @param start the absolute start bp
	 * @param end the absolute end bp (exclusive)
	 */
	public DNAregion(String chromosome, long start, long end){
		initParams(chromosome, start, end);
	}
	
	
	/**
	 * initialise params
	 * @param chromosome the name of the chromosome where the region resides
	 * @param start the absolute start bp
	 * @param end the absolute end bp (exclusive)
	 */
	public void initParams(String chromosome, long start, long end){
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
		scaleFactor = 1.0;
		direction = 0;
		probability=1.0;
		this.hasProbability = false;
		this.hasDirection = false;
	}
	
	/**
	 * loads the data from a dna region into the current one;
	 * @param region
	 */
	public void loadRegion(DNAregion region){
		this.chromosome = region.chromosome;
		this.start = region.start;
		this.end = region.end;
		this.scaleFactor = region.scaleFactor;
		this.direction = region.direction;
		this.probability = region.probability;
		this.hasDirection = region.hasDirection;
		this.hasProbability = region.hasProbability;
	}
	
	/**
	 * constructor that loads data from a different region
	 * @param region
	 */
	public DNAregion(DNAregion region){
		loadRegion(region);
	}
	
	/**
	 * class constructor which initialises the parameters with the supplied values and attempts to parse the text
	 * @param description
	 * @param chromosome the name of the chromosome where the region resides
	 * @param start the absolute start bp
	 * @param end the absolute end bp (exclusive)
	 */
	public DNAregion(String description, String chromosome,long start, long end, boolean hasProbability, boolean hasDirection){
		initParams(chromosome,start,end);
		parseText(description, hasProbability, hasDirection);
		this.hasDirection = hasDirection;
		this.hasProbability = hasProbability;
	} 
	
	public void parseText(String description, boolean hasProbability, boolean hasDirection){
		String[] buffer, bufferLength;
		if(description!=null && !description.isEmpty()){
			if(description.contains(Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER)){
				buffer = description.split(Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER);
				
				if(buffer!=null && buffer.length>=2){
					// get chromosome name
					if(buffer[0]!=null && !buffer[0].isEmpty()){
						this.chromosome = buffer[0].trim();
					}
					// get start and end points
					if(buffer[1]!=null && !buffer[1].isEmpty()){
						bufferLength = buffer[1].split(Constants.FASTA_FILE_DESCRIPTION_INTERVAL_DELIMITER_REGEX);				
						if(bufferLength!=null && bufferLength.length==2){
							this.start = Utils.parseLong(bufferLength[0], this.start);
							this.end = Utils.parseLong(bufferLength[1], this.end);
						} 
					} 
					if(buffer.length>=3){
						if(hasProbability){
							this.probability = Utils.parseDouble(buffer[2], this.probability);
						} else if(hasDirection){
							this.direction = Utils.parseInteger(buffer[2], this.direction);
						}
					}
				} 
				
			} 
		}

	}
	
	
	/**
	 * returns that string that describes the object
	 */
	public String toString(){
		StringBuffer strBuf = new StringBuffer();
		strBuf.append(this.chromosome);
		strBuf.append(Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER);
		strBuf.append(this.start);
		strBuf.append(Constants.FASTA_FILE_DESCRIPTION_INTERVAL_DELIMITER);
		strBuf.append(this.end);
		if(this.hasProbability){
			strBuf.append(Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER);
			strBuf.append(this.probability);
		}
		if(this.hasDirection){
			strBuf.append(Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER);
			strBuf.append(this.direction);
		}
		return strBuf.toString();
	}
	
	/**
	 * checks whether the DNA region is undefined
	 * @return
	 */
	public boolean isUndefined(){
		return this.start == Constants.NONE || this.end == Constants.NONE; 
	}
	
	/**
	 * compares the current DNA region to the one provided as an argument and returns true if they are equal
	 * @param region
	 * @return
	 */
	public boolean equals(DNAregion region){
		return this.start == region.start && this.end==region.end && this.chromosome.equals(region.chromosome);
	}
	
	/**
	 * returns the size in bp of a DNA region
	 * @return
	 */
	public int size(){
		return (int) (isUndefined()?0:(int)this.end-this.start);
	}
	
	/**
	 * a specific DNA region is a defined one, not equal to the entire DNA and that has a size > 0;
	 * @param defaultRegion
	 * @return
	 */
	public boolean isSpecific(DNAregion defaultRegion){
		return !isUndefined() &&  !equals(defaultRegion) && size()>0;
		
	}
	

	/**
	 * returns true if current region includes the sub region parsed as a parameter to the file.
	 * @param region
	 * @return
	 */
	public boolean includes(DNAregion region){
		boolean result =false;
		
		if(this.start<=region.start && this.end>=region.end){
			result=true;
		}
		
		return result;
	}
	
	
	/**
	 * returns true if current region includes the sub region parsed as a parameter to the file.
	 * @param region
	 * @return
	 */
	public boolean isSubSequenceOf(DNAregion region){
		boolean result =false;
		
		if(this.start>region.start && this.end<region.end){
			result=true;
		}
		
		return result;
	}
	
}
