package objects;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Stack;

import utils.Constants;
import utils.RPNtree;




public class TargetSitesGroup   implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -8426350303572761221L;
	public int groupID;
	public ArrayList<Integer> targetSitesID;
	public RPNtree rpnTree;
	public String text;
	public double lastTimeUpdate;
	public double timeOccupied;
	public double firstTimeReached;
	public int timesReached;
	public boolean isOccupied;
	
	/**
	 * class contructor
	 * @param groupID the ID of the group
	 * @param isAND is true if this is a group where all sites are connected by AND logic
	 * @param isOR is true if this is a group where all sites are connected by OR logic
	 */
	public TargetSitesGroup(int groupID, String text){
		this.groupID = groupID;
		targetSitesID = new ArrayList<Integer>();
		this.lastTimeUpdate = 0;
		this.timeOccupied = 0;
		this.firstTimeReached = Constants.NONE;
		this.timesReached = 0;

		this.text = text;
		isOccupied = false;
	}
	
	/**
	 * adds a target site id to this group
	 * @param targetSiteID
	 */
	public void addTargetSite(int targetSiteID){
		this.targetSitesID.add(targetSiteID);
	}
	
	/**
	 * generate the logic expression tree for the site
	 * @param str
	 */
	public boolean generateRPN(String str){
		boolean result=true;
	
		if (str == null || str.isEmpty()){
    	 		result=false;
		} else{
			
    			char[] in = str.toCharArray();
    			Character c;
    			Stack<Character> stack = new Stack<Character>();
    			StringBuilder out = new StringBuilder();    		
    			for (int i = 0; i < in.length; i++)
    		        switch (in[i]) {
    		        case '!':
    		        case '+':
    		                while (!stack.empty() && (stack.peek() == '*' || stack.peek() == '/')){
    		                        	c = stack.pop();
    		                			out.append(' ').append(c);
    		                }
    		        case '*':
    		                out.append(' ');
    		        case '(':
    		                stack.push(in[i]);
    		        case ' ':
    		                break;
    		        case ')':
    		                while (!stack.empty() && stack.peek() != '('){
	                        	c = stack.pop();
	                			out.append(' ').append(c);
    		                }
    		                if (!stack.empty())
    		                        stack.pop();
    		                break;
    		        default:
    		                out.append(in[i]);
    		                break;
    		        }

    		    while (!stack.isEmpty()){
                	c = stack.pop();
        			out.append(' ').append(c);
                }
    			
    			Stack<String> rpn;
    		    rpn = new Stack<String>();
    		    rpn.addAll(Arrays.asList(out.toString().trim().split("[ \t]+")));
    		    rpnTree = new RPNtree(rpn);
		}
    		return result;
	}
	
	/*private boolean evalRPN(Stack<String> rpn, boolean[] occupancy) {
		String element = rpn.pop();
	    int pos;
	    boolean x,y;
	    try{
	    		pos = Integer.parseInt(element);
	    		x= occupancy[pos];
	    	} catch (Exception e){
	    		x=false;
	    		if(element.equals("+") && rpn.size()>=2){ 
	    			y = evalRPN(rpn, occupancy);
		    		x = evalRPN(rpn, occupancy);
	    			x=x || y;
	    		} else if(element.equals("*")&& rpn.size()>=2){ 
	    			y = evalRPN(rpn, occupancy);
		    		x = evalRPN(rpn, occupancy);
	    			x=x && y;
	    		} else if(element.equals("!") && rpn.size()>=1){ 
		    		x = evalRPN(rpn, occupancy);
	    			x=!x;
	    		} else{
	    			System.out.println("Unknown logic"+element);
	    		}
	   }
	   
	   return x;
		
	}
	
	
	@SuppressWarnings("unchecked")
	public boolean evalRPN(boolean[] occupancy) {
		Stack<String> buffer =(Stack<String>) rpn.clone();
		return evalRPN(buffer, occupancy);
	}*/
	
	public boolean evaluateRPNTree(boolean[] occupancy) {
		return this.rpnTree.evaluateRPNTree(occupancy);
	}
	
	/**
	 * update target site statistics
	 * @param time the current time
	 */
	public void updateTimesReachedStatistics(double time){
		if(this.firstTimeReached == Constants.NONE){
			this.firstTimeReached = time;
		}
		this.timesReached++;
	}
	
	
	/**
	 * update the ts occupancy time
	 * @param time the current time
	 */
	public void updateOccupancyStatistics(double time){
		this.timeOccupied+=(time-this.lastTimeUpdate);
	}

	/**
	 * updates the time the state of the group was last changed. 
	 * @param time
	 */
	public void updateLastTimeUpdate(double time){
		this.lastTimeUpdate = time;
	}
	
	/**
	 * generates the string
	 */
	public String toString(){
		return text;
	}
}
