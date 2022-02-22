package utils;

import java.io.Serializable;

public class RPNnode   implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -6844290063073910023L;
	public int argument;
	public String operator;
	public boolean isOperator;
	public boolean isBinary;
	public RPNnode Left;
	public RPNnode Right;
	public int size;
	/**
	 * this is a node in a RPN tree
	 * class constructor based on a string
	 * @param str
	 */
	public RPNnode(String str){
		if(CellUtils.tgsl.isOperator(str)){
			isOperator = true;
			isBinary = CellUtils.tgsl.isBinary[CellUtils.tgsl.delimiterIDs.get(str)];
			operator=str;
			argument=Constants.NONE;
		} else{
			isOperator = false;
			isBinary=false;
			operator="";
			argument=Utils.parseInteger(str, Constants.NONE);
		}
		size=0;
	}
	
	public String toString(){
		return this.isOperator? operator: argument+"";		
	}
	
}
