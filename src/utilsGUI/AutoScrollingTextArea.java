package utilsGUI;

import javax.swing.JTextArea;

/**
 * auto-scrolling text areas
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class AutoScrollingTextArea  extends JTextArea{

	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public AutoScrollingTextArea(int rows, int cols) {
		super(rows,cols);
	}
	
	public AutoScrollingTextArea() {
		super("");
	}
	
	public AutoScrollingTextArea(String text) {
		super(text);
	}
	
	public void append(String text) {
		int pos = this.getDocument().getLength();
		super.append(text);
		
		this.setCaretPosition(pos);
		
	}
	
}
