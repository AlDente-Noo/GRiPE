package utilsGUI;

import java.awt.Dimension;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class LabelledComboBox  extends JPanel{
	private static final long serialVersionUID = 1L;
	private final JLabel label;
    private final JComboBox component;
    
    /**
     * @param lab text to go in the label
     */
    public LabelledComboBox(String lab, String toolTipText, String[] values, String defaultValue) {
        super();
        setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
        label = new JLabel(lab);
        label.setHorizontalAlignment(JLabel.LEFT);
        add(label);
        label.setToolTipText(lab);
        component = new JComboBox(values);
		component.setSelectedItem(defaultValue);
        component.setToolTipText(toolTipText);
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
     */
    public String getValue(){
    	return (String) component.getSelectedItem();
    }
    
    /**
     * sets the displayed value
     */
    public void setValue(String value){
    		component.setSelectedItem(value);
    }
    
    
    
    /**
     * adds an listener
     */
    public void addActionListener(ActionListener e){
    	component.addActionListener(e);
    }
}
