package objects;

import java.io.Serializable;
import java.util.HashMap;
import utils.Constants;

/**
 * This class creates an object which store informations about the nucleotides used
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class BasePairs  implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 5214691368363833137L;
	public String[] bps = {"A", "C", "G", "T", "N"}; // the letters used in the DNA strands  
	public String[] complementaryBPs = {"T", "G", "C", "A", "N"}; // the letters used in the DNA strands  

	public int numberOfBP = 4;
	
	//list of mapping of bp to a number
	public HashMap<String, Byte> bpsID;
	
	/**
	 * class constructor
	 */
	public BasePairs(){
		bpsID = buildBP(bps);
	}
	
	/**
	 * constructs a map between the chars of the nucleotides and the IDs 
	 * @param bp an array of chars that lists all nucleotides used
	 * @return the HasMap which maps a number for each nucleotide letter
	 */
	public HashMap<String, Byte> buildBP(String[] bp){
		HashMap<String, Byte> ids = new HashMap<String, Byte>();
		
		for (int i= 0; i< bp.length; i++){
			ids.put(bp[i], (byte) i);
		}
		return ids;
	}
	
	
	/**
	 * gets the complement nucleotid letter of the one supplied as a parameter
	 * @param a the letter for which the method looks for a complement
	 * @return the complement letter
	 */
	public String getComplement(String a){
		String result="";
		
		if(bpsID.containsKey(a)){
			result = complementaryBPs[bpsID.get(a)];
		}
		
		return result;
	}
	
	
	
	/**
	 * returns the ID of the complementary nucleotide
	 * @param bp the id of the base pair for which the method looks for a complement
	 * @result the id of the complementary nucleotide
	 */
	public byte getComplement(byte bp){
		byte result  = (byte) Constants.NONE;
		
		if(bp>=0 && bp<= bps.length){
			result =  bpsID.get(complementaryBPs[bp]);
		}
		
		return result;
	}
	
	/**
	 * returns the ID of the none bp
	 * @return
	 */
	public byte getANYID(){
		return bpsID.get("N");
	}
	
}
