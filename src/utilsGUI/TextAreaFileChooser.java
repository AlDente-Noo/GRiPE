package utilsGUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JTextField;


/**
 * text area + file chooser object
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TextAreaFileChooser extends JPanel{
	private static final long serialVersionUID = 1L;
	private JTextField path;
	private JButton select;
	private JButton clear;
	private boolean isFile;
    JFileChooser chooser;

	
	
	public TextAreaFileChooser(boolean isFile, boolean withClear){
		super();
		

		

		

		
		this.setLayout(new BorderLayout());
		path = new JTextField(GUIconstants.TEXTAREA_WIDTH);
		path.setEditable(false);
		this.add(path, BorderLayout.WEST);
		

		chooser = new JFileChooser();

		InputStream in;
		ImageIcon buttonSelectImg = new ImageIcon();
	    File buttonSelectImgImgFile = new File(GUIconstants.OPEN_IMAGE);
		if(!buttonSelectImgImgFile.exists()){
			in = this.getClass().getClassLoader().getResourceAsStream(GUIconstants.OPEN_IMAGE);
			try {
				if(in.available() >0){
					buttonSelectImg = new ImageIcon(ImageIO.read(in));
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else{
			buttonSelectImg= new ImageIcon(GUIconstants.OPEN_IMAGE);
		}		
		select = new JButton(buttonSelectImg);
		select.setBackground(Color.white);
		select.setMaximumSize(new Dimension(buttonSelectImg.getIconWidth(), buttonSelectImg.getIconHeight()));
		this.isFile = isFile;
		addActionListener(select, this.isFile, true);
		this.add(select, BorderLayout.CENTER);
		
		if(withClear){
			ImageIcon buttonTrashImg = new ImageIcon();
			File buttonTrashImgFile = new File(GUIconstants.TRASH_IMAGE);
			if(!buttonTrashImgFile.exists()){
				in = this.getClass().getClassLoader().getResourceAsStream(GUIconstants.TRASH_IMAGE);
				try {
					if(in.available() >0){
						buttonTrashImg = new ImageIcon(ImageIO.read(in));
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			} else{
				buttonTrashImg= new ImageIcon(GUIconstants.TRASH_IMAGE);
			}		
			
			
			clear = new JButton(buttonTrashImg);
			clear.setBackground(Color.white);
			clear.setMaximumSize(new Dimension(buttonTrashImg.getIconWidth(), buttonTrashImg.getIconHeight()));
			addActionListener(clear, this.isFile, false);
			
			this.add(clear, BorderLayout.EAST);
		}

	
		
		
	}
		
	
	private void addActionListener(JButton button, boolean isFile, boolean isSelect){
		if(isSelect){
			if(!isFile){
				button.addActionListener(new ActionListener() {
					   public void actionPerformed(ActionEvent ae) {
				    	 chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
					     int option = chooser.showOpenDialog(TextAreaFileChooser.this);
					     if (option == JFileChooser.APPROVE_OPTION && chooser.getSelectedFile()!=null) {
					    	 	setValue(chooser.getSelectedFile().getAbsolutePath());
					     } 
					     
					   }
				});
			} else{
				button.addActionListener(new ActionListener() {
					   public void actionPerformed(ActionEvent ae) {
				    	 chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
					     int option = chooser.showOpenDialog(TextAreaFileChooser.this);
					     if (option == JFileChooser.APPROVE_OPTION && chooser.getSelectedFile()!=null) {
					    	 	setValue(chooser.getSelectedFile().getAbsolutePath());
					     } 
					     
				
					    	
					   }
				});
			}
		} else{
			button.addActionListener(new ActionListener() {
				   public void actionPerformed(ActionEvent ae) {
			    	    	 	setValue("");
				   }
			});
		}
	}

	/**
	 * sets a value of this component
	 * @param path
	 */
	public void setValue(String path){
		this.path.setText(path);
		this.path.setToolTipText(this.path.getText());
		
		File file=new File(path);
		boolean exists = file.exists();
	
		if(exists){
			if(this.isFile){
				this.chooser.setCurrentDirectory(file.getParentFile());
			} else{
				this.chooser.setCurrentDirectory(file);
			}
		}
	}
	
	
	/**
	 * returns the value of this component
	 * @return path
	 */
	public String getValue(){
		return this.path.getText();
	}
	
	/**
	 * sets the tool tip text
	 */
	public void setToolTipText(String text){
		select.setToolTipText(text);
	}
	
    /**
     * sets the enable status
     * @param e
     */
    public void setEnable(boolean e){
    	select.setEnabled(e);
    }
    
    /**
     * adds an listener
     * @param e
     */
    public void addActionListener(ActionListener e){
    	select.addActionListener(e);
    }
}
