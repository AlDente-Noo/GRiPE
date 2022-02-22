package utilsGUI;

import javax.swing.JLabel;

/**
 * GUI parameter
 * @author radu
 *
 * @param <E> the type of GUI element 
 */
public class GUIParameter <E>{
	public JLabel label;
	public E component;
	
	public GUIParameter(String text, E component){
		label = new JLabel(text);
		label.setToolTipText(text);
		
		this.component = component;
	}
	
}
