package objects;

import java.io.Serializable;

/**
 * a generic parameter read from a file
 * @author n.r.zabet@gen.cam.ac.uk
 *
 * @param <E> they value type (Double, String...)
 */
public class Parameter<E>  implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3925194356792686450L;
	public E value;
	public String label;
	public String description;
	public String name;
	public String category;
	
	/**
	 * class constructor
	 * @param data the data stored by the param
	 * @param label the label used in the GUI
	 * @param description the description
	 * @param name the name used in the file
	 */
	public Parameter(String name, String label, String description, String category, E value){
		this.value = value;
		this.label = label;
		this.description = description;
		this.name = name;
		this.category = category;
	}
	
	/**
	 * empty constructor
	 */
	public Parameter(){
		this.value = null;
		this.label = "";
		this.description = "";
		this.name = "";
		this.category = "";
	}
}
