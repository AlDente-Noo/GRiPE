package objects;

import java.io.Serializable;
import java.util.HashMap;

import utils.Constants;

public class TargetSitesGroupLogic   implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -2293416629979113990L;
	public String[] logics = {"NONE", "AND","OR","NOT"}; // the type of logic functions
	public String[] delimiters = {"", "*","+","!"}; // the type of logic functions
	public boolean[] isBinary = {false, true,true,false}; // true only if it is a binary operator
	public HashMap<String, Integer> logicIDs;
	public HashMap<String, Integer> delimiterIDs;

	public TargetSitesGroupLogic(){
		logicIDs = buildLogicIDs(logics);
		delimiterIDs = buildLogicIDs(delimiters);
	}
	
	
	/**
	 * constructs a map between the chars of the nucleotides and the IDs 
	 * @param bp an array of chars that lists all nucleotides used
	 * @return the HasMap which maps a number for each nucleotide letter
	 */
	public HashMap<String, Integer> buildLogicIDs(String[] logics){
		HashMap<String, Integer> ids = new HashMap<String, Integer>();
		
		for (int i= 0; i< logics.length; i++){
			ids.put(logics[i], i);
		}
		return ids;
	}
	
	
	/**
	 * returns the NONE ID
	 * @return
	 */
	public int getNONEid(){
		return logicIDs.get("NONE");
	}
	
	/**
	 * gets the ID of the delimiter
	 * @param str
	 * @return
	 */
	public int getDeliumiterID(String str){
		int result=Constants.NONE;
		for(int i=0;i<this.delimiters.length && result == Constants.NONE;i++){
			if(this.delimiters[i]==str){
				result=i;
			}
		}
		return result;
	}
	
	/**
	 * checks whether a string is a operator
	 * @param str
	 * @return
	 */
	public boolean isOperator(char str){
		boolean result=false;
		for(int i=0;i<this.logics.length && !result;i++){
			if(this.logics[i].equals(str)){
				result=true;
			}
		}
		return result;
	}
	
	/**
	 * checks whether a string is a operator
	 * @param str
	 * @return
	 */
	public boolean isOperator(String str){
		boolean result=false;
		for(int i=0;i<this.delimiters.length && !result;i++){
			if(this.delimiters[i].equals(str)){
				result=true;
			}
		}
		
		for(int i=0;i<this.logics.length && !result;i++){
			if(this.logics[i].equals(str)){
				result=true;
			}
		}
		
		return result;
	}
	
	/**
	 * replaces the text logic  operators to the 1 char operators 
	 * @param str
	 * @return
	 */
	public String replaceOperators(String str){
		String result = str;
		
		
		for(int i=0;i<delimiters.length;i++){
			result = result.replaceAll(this.logics[i], this.delimiters[i]);
		}
		
		return result;
	}
}
