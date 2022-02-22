package utilsGUI;

import java.awt.Dimension;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JPanel;

/**
 * labbeled check box
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class LabelledCheckBox extends JPanel{
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JCheckBox component;
    
    /**
     * @param lab text to go in the label
     * @param chars - the size of the text input area
     */
    public LabelledCheckBox(String lab, String toolTipText, boolean value) {
        super();
        setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));

        component = new JCheckBox(lab, value);
        component.setToolTipText(toolTipText);
        add(component);
    }
    
    /**
     * @return the label width in pixels
     */
    public int getLabelWidth() {
        return component.getPreferredSize().width;
    }
    
    /**
     * @param width - the width of the label in pixels
     */
    public void setLabelWidth(int width) {
        Dimension d = component.getPreferredSize();
        d.width = width;
        component.setPreferredSize(d);
    }
    
    /**
     * returns true whether the check box is checked or false otherwise 
     * @return
     */
    public boolean getValue(){
    	
    	return component.isSelected();
    }
    
    /**
     *  sets true whether the check box is checked or false otherwise 
     * @param text
     */
    public void setValue(boolean value){
    	component.setSelected(value);
    }
    
    /**
     * sets the enable status
     * @param e
     */
    public void setEnable(boolean e){
    	component.setEnabled(e);
    }
    
    
    /**
     * adds an listener
     * @param e
     */
    public void addActionListener(ActionListener e){
    	component.addActionListener(e);
    }
    
}
