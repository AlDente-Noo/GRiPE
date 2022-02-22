package utilsGUI;

import java.awt.Dimension;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * labelled file chooser
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class LabelledFileChooser extends JPanel{
	private static final long serialVersionUID = 1L;
	private JLabel label;
    private TextAreaFileChooser component;
    
    /**
     * @param lab text to go in the label
     * @param chars - the size of the text input area
     */
    public LabelledFileChooser(String lab, int columns, String toolTipText, String value, boolean isFile, boolean withClear) {
        super();
        setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
        label = new JLabel(lab+" ");
        label.setHorizontalAlignment(JLabel.LEFT);
        add(label);
        label.setToolTipText(lab);
        component = new TextAreaFileChooser(isFile, withClear);
        component.setToolTipText(toolTipText);
        add(component);
        component.setValue(value);
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
     * @return
     */
    public String getValue(){
    	return component.getValue();
    }
    
    /**
     * sets the displayed value
     * @param text
     */
    public void setValue(String text){
    	component.setValue(text);
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
