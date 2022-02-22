package utilsGUI;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;

import objects.InputParameters;

public class OutputParameters extends JPanel{
	private static final long serialVersionUID = 1L;

	private JPanel componentsStack;


	//SIMULATION-OUTPUT PARAMATERS
	public LabelledFileChooser OUTPUT_FOLDER;
	public LabelledText OUTPUT_FILENAME;
	public LabelledDouble PRINT_INTERMEDIARY_RESULTS_AFTER;	
	public LabelledCheckBox PRINT_FINAL_OCCUPANCY;
	public LabelledCheckBox DEBUG_MODE;
	public LabelledText OUTPUT_TF;
	public LabelledInteger OUTPUT_TF_POINTS;
	public LabelledCheckBox FOLLOW_TS;
	public LabelledCheckBox OUTPUT_AFFINITY_LANDSCAPE;
	public LabelledCheckBox OUTPUT_BINDING_ENERGY;
	public LabelledCheckBox OUTPUT_DNA_OCCUPANCY;
	public LabelledCheckBox DNA_OCCUPANCY_FULL_MOLECULE_SIZE;
	public LabelledCheckBox OUTPUT_SLIDING_LENGTHS;
	public LabelledInteger WIG_STEP;
	public LabelledDouble WIG_THRESHOLD;	


	
	public OutputParameters(InputParameters ip){
		this.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH,GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT));
		this.setLayout(new FlowLayout());
		componentsStack = new JPanel(new GridLayout(0,1, GUIconstants.GRID_HGAP, GUIconstants.GRID_WGAP));
		componentsStack.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT,GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH));
		
		
		
		JLabel label1,label2;
		label1 = new JLabel(GUIconstants.SIMULATION_AREA_SIMULATION_GENERAL_PARAMATERS);
		label2 = new JLabel(GUIconstants.SIMULATION_AREA_SIMULATION_OUTPUT_PARAMATERS);
		label1.setForeground(Color.LIGHT_GRAY);
		label2.setForeground(Color.LIGHT_GRAY);

		
		
		
		//simulation output params
		OUTPUT_FILENAME = new LabelledText(ip.OUTPUT_FILENAME.label, GUIconstants.TEXTAREA_WIDTH, ip.OUTPUT_FILENAME.description, ip.OUTPUT_FILENAME.value); 
		OUTPUT_FOLDER = new LabelledFileChooser(ip.OUTPUT_FOLDER.label,GUIconstants.TEXTAREA_WIDTH,ip.OUTPUT_FOLDER.description,ip.OUTPUT_FOLDER.value, false, true);
		PRINT_INTERMEDIARY_RESULTS_AFTER = new LabelledDouble(ip.PRINT_INTERMEDIARY_RESULTS_AFTER.label,GUIconstants.TEXTAREA_WIDTH,ip.PRINT_INTERMEDIARY_RESULTS_AFTER.description,ip.PRINT_INTERMEDIARY_RESULTS_AFTER.value);	
		PRINT_FINAL_OCCUPANCY = new LabelledCheckBox(ip.PRINT_FINAL_OCCUPANCY.label, ip.PRINT_FINAL_OCCUPANCY.description, ip.PRINT_FINAL_OCCUPANCY.value); 
		DEBUG_MODE = new LabelledCheckBox(ip.DEBUG_MODE.label, ip.DEBUG_MODE.description, ip.DEBUG_MODE.value); 
		OUTPUT_TF = new LabelledText(ip.OUTPUT_TF.label, GUIconstants.TEXTAREA_WIDTH, ip.OUTPUT_TF.description, ip.OUTPUT_TF.value); 
		OUTPUT_TF_POINTS = new LabelledInteger(ip.OUTPUT_TF_POINTS.label,GUIconstants.TEXTAREA_WIDTH,ip.OUTPUT_TF_POINTS.description,ip.OUTPUT_TF_POINTS.value.intValue());
		FOLLOW_TS = new LabelledCheckBox(ip.FOLLOW_TS.label, ip.FOLLOW_TS.description, ip.FOLLOW_TS.value); 
		OUTPUT_AFFINITY_LANDSCAPE = new LabelledCheckBox(ip.OUTPUT_AFFINITY_LANDSCAPE.label, ip.OUTPUT_AFFINITY_LANDSCAPE.description, ip.OUTPUT_AFFINITY_LANDSCAPE.value);
		OUTPUT_BINDING_ENERGY = new LabelledCheckBox(ip.OUTPUT_BINDING_ENERGY.label, ip.OUTPUT_BINDING_ENERGY.description, ip.OUTPUT_BINDING_ENERGY.value);
		OUTPUT_DNA_OCCUPANCY = new LabelledCheckBox(ip.OUTPUT_DNA_OCCUPANCY.label, ip.OUTPUT_DNA_OCCUPANCY.description, ip.OUTPUT_DNA_OCCUPANCY.value); 
		DNA_OCCUPANCY_FULL_MOLECULE_SIZE = new LabelledCheckBox(ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.label, ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.description, ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value); 
		OUTPUT_SLIDING_LENGTHS = new LabelledCheckBox(ip.OUTPUT_SLIDING_LENGTHS.label, ip.OUTPUT_SLIDING_LENGTHS.description, ip.OUTPUT_SLIDING_LENGTHS.value); 
		WIG_STEP = new LabelledInteger(ip.WIG_STEP.label,GUIconstants.TEXTAREA_WIDTH,ip.WIG_STEP.description,ip.WIG_STEP.value.intValue());
		WIG_THRESHOLD = new LabelledDouble(ip.WIG_THRESHOLD.label,GUIconstants.TEXTAREA_WIDTH,ip.WIG_THRESHOLD.description,ip.WIG_THRESHOLD.value);	

		
		
		resetLabelsWidth();
		
		componentsStack.add(OUTPUT_FILENAME);
		componentsStack.add(OUTPUT_FOLDER);
		componentsStack.add(PRINT_INTERMEDIARY_RESULTS_AFTER);
		componentsStack.add(PRINT_FINAL_OCCUPANCY);		
		componentsStack.add(DEBUG_MODE);
		componentsStack.add(OUTPUT_TF);
		componentsStack.add(OUTPUT_TF_POINTS);
		componentsStack.add(FOLLOW_TS);
		componentsStack.add(OUTPUT_AFFINITY_LANDSCAPE);
		componentsStack.add(OUTPUT_BINDING_ENERGY);	
		componentsStack.add(OUTPUT_DNA_OCCUPANCY);
		componentsStack.add(DNA_OCCUPANCY_FULL_MOLECULE_SIZE);
		componentsStack.add(OUTPUT_SLIDING_LENGTHS);
		componentsStack.add(WIG_STEP);
		componentsStack.add(WIG_THRESHOLD);
		
		
		this.add(componentsStack);
		
	}
	
	
	/**
	 * resets the labels width
	 */
	private void resetLabelsWidth(){
		//SIMULATION PARAMATERS
		int max = OUTPUT_FOLDER.getLabelWidth();
		
		if(OUTPUT_FILENAME.getLabelWidth() > max){
			max = OUTPUT_FILENAME.getLabelWidth();
		}	
		
		if(PRINT_INTERMEDIARY_RESULTS_AFTER.getLabelWidth() > max){
			max = PRINT_INTERMEDIARY_RESULTS_AFTER.getLabelWidth();
		}	
		
		if(PRINT_FINAL_OCCUPANCY.getLabelWidth() > max){
			max = PRINT_FINAL_OCCUPANCY.getLabelWidth();
		}	
		
		if(OUTPUT_TF.getLabelWidth() > max){
			max = OUTPUT_TF.getLabelWidth();
		}
		
		if(OUTPUT_TF_POINTS.getLabelWidth() > max){
			max = OUTPUT_TF_POINTS.getLabelWidth();
		}	
		
		if(WIG_STEP.getLabelWidth() > max){
			max = WIG_STEP.getLabelWidth();
		}	
		
		
		if(WIG_THRESHOLD.getLabelWidth() > max){
			max = WIG_THRESHOLD.getLabelWidth();
		}	

		//SIMULATION-OUTPUT PARAMATERS
		OUTPUT_FOLDER.setLabelWidth(max);
		OUTPUT_FILENAME.setLabelWidth(max);
		PRINT_INTERMEDIARY_RESULTS_AFTER.setLabelWidth(max);	
		PRINT_FINAL_OCCUPANCY.setLabelWidth(max);	
		OUTPUT_TF.setLabelWidth(max);
		OUTPUT_TF_POINTS.setLabelWidth(max);
		WIG_STEP.setLabelWidth(max);
		WIG_THRESHOLD.setLabelWidth(max);
	}
}
