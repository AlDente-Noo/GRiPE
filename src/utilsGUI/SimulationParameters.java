package utilsGUI;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;

import objects.InputParameters;

/**
 * tabbed panel with simulation parameters
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class SimulationParameters extends JPanel{
	private static final long serialVersionUID = 1L;

	//SIMULATION PARAMATERS
	public LabelledDouble STOP_TIME;
	public LabelledInteger ENSEMBLE_SIZE;
	public LabelledInteger RANDOM_SEED;
	public LabelledInteger DNA_SECTOR_SIZE;
	public LabelledInteger EVENT_LIST_SUBGROUP_SIZE;

	
	public SimulationParameters(InputParameters ip){
		this.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH,GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT));
		this.setLayout(new FlowLayout());
		JPanel componentsStack = new JPanel(new GridLayout(0, 1, GUIconstants.GRID_HGAP, GUIconstants.GRID_VGAP));
		componentsStack.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH, GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT));
		
		
		
		JLabel label1,label2;
		label1 = new JLabel(GUIconstants.SIMULATION_AREA_SIMULATION_GENERAL_PARAMATERS);
		label2 = new JLabel(GUIconstants.SIMULATION_AREA_SIMULATION_OUTPUT_PARAMATERS);
		label1.setForeground(Color.LIGHT_GRAY);
		label2.setForeground(Color.LIGHT_GRAY);

		

		//simulation params
		STOP_TIME = new LabelledDouble(ip.STOP_TIME.label,GUIconstants.TEXTAREA_WIDTH,ip.STOP_TIME.description,ip.STOP_TIME.value);	
		ENSEMBLE_SIZE = new LabelledInteger(ip.ENSEMBLE_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.ENSEMBLE_SIZE.description, ip.ENSEMBLE_SIZE.value);
		RANDOM_SEED = new LabelledInteger(ip.RANDOM_SEED.label,GUIconstants.TEXTAREA_WIDTH,ip.RANDOM_SEED.description, ip.RANDOM_SEED.value);
		DNA_SECTOR_SIZE = new LabelledInteger(ip.DNA_SECTOR_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_SECTOR_SIZE.description, ip.DNA_SECTOR_SIZE.value);
		EVENT_LIST_SUBGROUP_SIZE = new LabelledInteger(ip.EVENT_LIST_SUBGROUP_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.EVENT_LIST_SUBGROUP_SIZE.description, ip.EVENT_LIST_SUBGROUP_SIZE.value);

		resetLabelsWidth();
		
		//simulation params
		componentsStack.add(STOP_TIME);
		componentsStack.add(ENSEMBLE_SIZE);
		componentsStack.add(RANDOM_SEED);
		componentsStack.add(DNA_SECTOR_SIZE);
		componentsStack.add(EVENT_LIST_SUBGROUP_SIZE);

		this.add(componentsStack);
	}
	
	
	/**
	 * resets the labels width
	 */
	private void resetLabelsWidth(){
		//SIMULATION PARAMATERS
		int max = STOP_TIME.getLabelWidth();
				
		if(RANDOM_SEED.getLabelWidth() > max){
			max = RANDOM_SEED.getLabelWidth();
		}
		
		if(ENSEMBLE_SIZE.getLabelWidth() > max){
			max = ENSEMBLE_SIZE.getLabelWidth();
		}
		
		if(DNA_SECTOR_SIZE.getLabelWidth() > max){
			max = DNA_SECTOR_SIZE.getLabelWidth();
		}
		
		if(EVENT_LIST_SUBGROUP_SIZE.getLabelWidth() > max){
			max = EVENT_LIST_SUBGROUP_SIZE.getLabelWidth();
		}
		
		
		
		
		
		
		
		
		
		//SIMULATION PARAMATERS
		STOP_TIME.setLabelWidth(max);
		ENSEMBLE_SIZE.setLabelWidth(max);
		RANDOM_SEED.setLabelWidth(max);
		DNA_SECTOR_SIZE.setLabelWidth(max);	
		EVENT_LIST_SUBGROUP_SIZE.setLabelWidth(max);	

	}
	
	
}

