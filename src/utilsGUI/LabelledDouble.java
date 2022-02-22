package utilsGUI;

import java.awt.Dimension;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeListener;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import javax.swing.BoxLayout;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 *  labeled double input box 
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class LabelledDouble  extends JPanel {

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
    public LabelledDouble(String lab, int columns, String toolTipText, double value) {
        super();
        setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
        label = new JLabel(lab+" ");
        label.setHorizontalAlignment(JLabel.LEFT);
        label.setToolTipText(lab);
        add(label);
        
        NumberFormat buffer = new DecimalFormat();
        buffer.setMaximumFractionDigits(GUIconstants.DOUBLE_PECESION_MAXIMUM);
        buffer.setMinimumFractionDigits(1);
        component =  new JFormattedTextField(buffer);
        component.setColumns(columns);
        component.setToolTipText(toolTipText);
        component.setHorizontalAlignment(JFormattedTextField.RIGHT);

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
     * returns the typed text
     * @return 0 if could not parse the text
     */
    public double getValue(){
    		return ((Number)component.getValue()).doubleValue();
    }
    
    /**
     * sets the displayed value
     * @param text
     */
    public void setValue(double value){
    	component.setValue(new Double(value));
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
