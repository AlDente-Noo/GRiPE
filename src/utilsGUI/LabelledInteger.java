package utilsGUI;

import java.awt.Dimension;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

import javax.swing.BoxLayout;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * labelled integer box
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class LabelledInteger extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private JLabel label;
	private JFormattedTextField component;
    
    /**
     * @param lab text to go in the label
     * @param chars - the size of the text input area
     */
    public LabelledInteger(String lab, int columns, String toolTipText, int value) {
        super();
        setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
        label = new JLabel(lab+" ");
        label.setHorizontalAlignment(JLabel.LEFT);
        label.setToolTipText(lab);
        add(label);
        component =  new JFormattedTextField(NumberFormat.getIntegerInstance());
        component.setHorizontalAlignment(JFormattedTextField.RIGHT);
        component.setColumns(columns);
        component.setToolTipText(toolTipText);
        setValue(value);
        add(component);
    }
    
    /**
     * @return the label width in pixels
     */
    public int getLabelWidth() {
        return label.getPreferredSize().width;
    }
    
    /**
     * @param width - the width of the label in pixels
     */
    public void setLabelWidth(int width) {
        Dimension d = label.getPreferredSize();
        d.width = width;
        label.setPreferredSize(d);
    }
    
    /**
     * returns the typed int
     * @return
     */
    public int getValue(){
    		return ((Number)component.getValue()).intValue();
    }
    
    /**
     * sets the displayed value
     * @param text
     */
    public void setValue(int value){
    	component.setValue(new Integer(value));
    }
    
    /**
     * sets the enable status
     * @param e
     */
    public void setEditable(boolean e){
    	component.setEditable(e);
    }
    
    /**
     * adds an listener
     * @param e
     */
    public void addActionListener(ActionListener e){
    	component.addActionListener(e);
    }
    
    /**
     * adds property change listener
     */
    public void addPropertyChangeListener(String propertyName, PropertyChangeListener listener){
    	component.addPropertyChangeListener(propertyName, listener);
    }
    
}
