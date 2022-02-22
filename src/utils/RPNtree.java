package utils;

import java.io.Serializable;
import java.util.Stack;

public class RPNtree   implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -8991196520972147021L;
	public RPNnode head;
	
	/**
	 * class constructor based on a stack of strings representing the RPN expression
	 * @param rpn
	 */
	public RPNtree(Stack<String> rpn){
		head=buildTree(rpn);
	}
	
	
	/**
	 * builds a tree recursivley
	 * @param rpn the stack of rpn
	 * @return
	 */
	private RPNnode buildTree(Stack<String> rpn){
		RPNnode node=null;
		if(rpn.size()>0){
			String element = rpn.pop();
			node = new RPNnode(element);
			if(node.isOperator){
	    			node.Left = buildTree(rpn);
	    			if(node.Left!=null){
	    				node.size++;
	    			}
	    			if(node.isBinary){
	    				node.Right = buildTree(rpn);
	    				if(node.Right!=null){
		    				node.size++;
		    			}
	    			}
			}
		}
	   return node;
	}
	
	/**
	 * evaluates recursively a RPN tree
	 * @param head
	 * @param occupancy
	 * @return
	 */
	private boolean evaluateRPNTree( RPNnode head, boolean[] occupancy){
		boolean x=false,y;
		//System.out.println(head +" -> "+head.size);
		if(head.isOperator){
			if(head.operator.equals("+") && head.size==2){ 
    				y = evaluateRPNTree(head.Right, occupancy);
    				x =  evaluateRPNTree(head.Left, occupancy);
    				//System.out.println(x+" OR "+y+" = "+(x||y));
    				x=x || y;
    			} else if(head.operator.equals("*") && head.size==2){ 
				y = evaluateRPNTree(head.Right, occupancy);
				x =  evaluateRPNTree(head.Left, occupancy);
				//System.out.println(x+" AND "+y+" = "+(x&&y));

				x=x && y;
    			} else if(head.operator.equals("!") && head.size==1){ 
				x =  evaluateRPNTree(head.Left, occupancy);
				//System.out.println("NOT"+x+" = "+(!x));

				x=!x;
    			} else{
    				//System.out.println("Unknown logic: "+head.operator);
    			}
		} else{
			x=occupancy[head.argument];
			//System.out.println(x);

		}
		
		return x;
	}
	
	/** 
	 * evaluates an RPN tree
	 * @param occupancy the vector storrying the current occupancy
	 * @return
	 */
	public boolean evaluateRPNTree(boolean[] occupancy){
		return this.evaluateRPNTree(head, occupancy);
	}
}
