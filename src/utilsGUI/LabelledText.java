package utilsGUI;
import javax.swing.*;


import java.awt.*;
import java.awt.event.ActionListener;

/**
 * labelled text box
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class LabelledText extends JPanel{
	private static final long serialVersionUID = 1L;
	private final JLabel label;
    private final JTextField component;
    
    /**
     * @param lab text to go in the label
     */
    public LabelledText(String lab, int columns, String toolTipText, String value) {
        super();
        setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
        label = new JLabel(lab);
        label.setHorizontalAlignment(JLabel.LEFT);
        add(label);
        label.setToolTipText(lab);
        component = new JTextField(value);
        component.setColumns(columns);
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
    	return component.getText();
    }
    
    /**
     * sets the displayed value
     */
    public void setValue(String text){
    	component.setText(text);
    }
    
    /**
     * adds an listener
     */
    public void addActionListener(ActionListener e){
    	component.addActionListener(e);
    }
}
