package utilsGUI;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;

import objects.InputParameters;
import utils.Constants;

/**
 * tabbed panel with DNA parameters
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class DNAParameters  extends JPanel{

	private static final long serialVersionUID = 1L;

	private JPanel componentsStack;

	
	//DNA PARAMETERS
	public LabelledFileChooser DNA_SEQUENCE_FILE;
	
	//DNA_RANDOM PARAMETERS
	public LabelledInteger DNA_LENGTH;
	public LabelledDouble DNA_PROPORTION_OF_A;
	public LabelledDouble DNA_PROPORTION_OF_T;
	public LabelledDouble DNA_PROPORTION_OF_C;
	public LabelledDouble DNA_PROPORTION_OF_G;
	public LabelledComboBox DNA_BOUNDARY_CONDITION;
	
	public DNAParameters(InputParameters ip){
		this.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH,GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT));
		this.setLayout(new FlowLayout());
		componentsStack = new JPanel(new GridLayout(0,1, GUIconstants.GRID_HGAP, GUIconstants.GRID_WGAP));
		componentsStack.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT,GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH));
		
		
		
		JLabel label1, label2;
		label1 = new JLabel(GUIconstants.SIMULATION_AREA_DNA_LOAD_PARAMATERS);
		label2 = new JLabel(GUIconstants.SIMULATION_AREA_DNA_RANDOM_PARAMATERS);
		label1.setForeground(Color.LIGHT_GRAY);
		label2.setForeground(Color.LIGHT_GRAY);
		

		

		//DNA PARAMETERS
		DNA_SEQUENCE_FILE = new LabelledFileChooser(ip.DNA_SEQUENCE_FILE.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_SEQUENCE_FILE.description,ip.DNA_SEQUENCE_FILE.value, true, true);

		//DNA_RANDOM PARAMETERS
		DNA_LENGTH = new LabelledInteger(ip.DNA_LENGTH.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_LENGTH.description,ip.DNA_LENGTH.value.intValue());
		DNA_PROPORTION_OF_A = new LabelledDouble(ip.DNA_PROPORTION_OF_A.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_PROPORTION_OF_A.description,ip.DNA_PROPORTION_OF_A.value);	
		DNA_PROPORTION_OF_T = new LabelledDouble(ip.DNA_PROPORTION_OF_T.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_PROPORTION_OF_T.description,ip.DNA_PROPORTION_OF_T.value);	
		DNA_PROPORTION_OF_C = new LabelledDouble(ip.DNA_PROPORTION_OF_C.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_PROPORTION_OF_C.description,ip.DNA_PROPORTION_OF_C.value);	
		DNA_PROPORTION_OF_G = new LabelledDouble(ip.DNA_PROPORTION_OF_G.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_PROPORTION_OF_G.description,ip.DNA_PROPORTION_OF_G.value);	

		
		String[] boundaryConditions = {Constants.DNA_FASTA_BOUNDARY_ABSORBING, Constants.DNA_FASTA_BOUNDARY_REFLEXIVE, Constants.DNA_FASTA_BOUNDARY_PERIODIC}; 
		DNA_BOUNDARY_CONDITION = new LabelledComboBox(ip.DNA_BOUNDARY_CONDITION.label,ip.DNA_BOUNDARY_CONDITION.description,boundaryConditions ,Constants.DNA_FASTA_BOUNDARY_REFLEXIVE);
		
		
		resetLabelsWidth();
		
		
		
		
		
		

		
		
		
		//DNA PARAMETERS
		componentsStack.add(label1);
		componentsStack.add(DNA_SEQUENCE_FILE);
		
		//DNA_RANDOM PARAMETERS
		componentsStack.add(label2);
		componentsStack.add(DNA_LENGTH);
		componentsStack.add(DNA_PROPORTION_OF_A);
		componentsStack.add(DNA_PROPORTION_OF_T);
		componentsStack.add(DNA_PROPORTION_OF_C);
		componentsStack.add(DNA_PROPORTION_OF_G);
		componentsStack.add(DNA_BOUNDARY_CONDITION);
		this.add(componentsStack);
		

	}
	
	
	/**
	 * resets the labels width
	 */
	private void resetLabelsWidth(){


		
		//DNA PARAMETERS
		int max = DNA_SEQUENCE_FILE.getLabelWidth();
				
		//DNA_RANDOM PARAMETERS
		if(DNA_LENGTH.getLabelWidth() > max){
			max = DNA_LENGTH.getLabelWidth();
		}

		if(DNA_PROPORTION_OF_A.getLabelWidth() > max){
			max = DNA_PROPORTION_OF_A.getLabelWidth();
		}
		
		if(DNA_PROPORTION_OF_T.getLabelWidth() > max){
			max = DNA_PROPORTION_OF_T.getLabelWidth();
		}
		
		if(DNA_PROPORTION_OF_C.getLabelWidth() > max){
			max = DNA_PROPORTION_OF_C.getLabelWidth();
		}			
		if(DNA_PROPORTION_OF_G.getLabelWidth() > max){
			max = DNA_PROPORTION_OF_G.getLabelWidth();
		}
		if(DNA_BOUNDARY_CONDITION.getLabelWidth() > max){
			max = DNA_BOUNDARY_CONDITION.getLabelWidth();
		}
		
		
		//DNA PARAMETERS
		DNA_SEQUENCE_FILE.setLabelWidth(max);
		
		//DNA_RANDOM PARAMETERS
		DNA_LENGTH.setLabelWidth(max);	
		DNA_PROPORTION_OF_A.setLabelWidth(max);	
		DNA_PROPORTION_OF_T.setLabelWidth(max);	
		DNA_PROPORTION_OF_C.setLabelWidth(max);
		DNA_PROPORTION_OF_G.setLabelWidth(max);	
		DNA_BOUNDARY_CONDITION.setLabelWidth(max);
	}
	
	
}
