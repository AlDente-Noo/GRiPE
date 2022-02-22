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

	private JPanel componentsStack;
	
	//SIMULATION PARAMATERS
	public LabelledDouble STOP_TIME;
	public LabelledInteger ENSAMBLE_SIZE;
	public LabelledInteger RANDOM_SEED;
	public LabelledInteger COMPUTED_AFFINITY_PRECISION;
	public LabelledInteger DNA_SECTOR_SIZE;
	public LabelledInteger EVENT_LIST_SUBGROUP_SIZE;
	public LabelledCheckBox EVENT_LIST_USES_FR;

	


	
	public SimulationParameters(InputParameters ip){
		this.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH,GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT));
		this.setLayout(new FlowLayout());
		componentsStack = new JPanel(new GridLayout(0,1, GUIconstants.GRID_HGAP, GUIconstants.GRID_WGAP));
		componentsStack.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT,GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH));
		
		
		
		JLabel label1,label2;
		label1 = new JLabel(GUIconstants.SIMULATION_AREA_SIMULATION_GENERAL_PARAMATERS);
		label2 = new JLabel(GUIconstants.SIMULATION_AREA_SIMULATION_OUTPUT_PARAMATERS);
		label1.setForeground(Color.LIGHT_GRAY);
		label2.setForeground(Color.LIGHT_GRAY);

		

		//simulation params
		STOP_TIME = new LabelledDouble(ip.STOP_TIME.label,GUIconstants.TEXTAREA_WIDTH,ip.STOP_TIME.description,ip.STOP_TIME.value);	
		ENSAMBLE_SIZE = new LabelledInteger(ip.ENSEMBLE_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.ENSEMBLE_SIZE.description,ip.ENSEMBLE_SIZE.value.intValue());
		RANDOM_SEED = new LabelledInteger(ip.RANDOM_SEED.label,GUIconstants.TEXTAREA_WIDTH,ip.RANDOM_SEED.description,ip.RANDOM_SEED.value.intValue());
		COMPUTED_AFFINITY_PRECISION = new LabelledInteger(ip.COMPUTED_AFFINITY_PRECISION.label,GUIconstants.TEXTAREA_WIDTH,ip.COMPUTED_AFFINITY_PRECISION.description,ip.COMPUTED_AFFINITY_PRECISION.value.intValue());
		DNA_SECTOR_SIZE = new LabelledInteger(ip.DNA_SECTOR_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.DNA_SECTOR_SIZE.description,ip.DNA_SECTOR_SIZE.value.intValue());
		EVENT_LIST_SUBGROUP_SIZE = new LabelledInteger(ip.EVENT_LIST_SUBGROUP_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.EVENT_LIST_SUBGROUP_SIZE.description,ip.EVENT_LIST_SUBGROUP_SIZE.value.intValue());
		EVENT_LIST_USES_FR = new LabelledCheckBox(ip.EVENT_LIST_USES_FR.label, ip.EVENT_LIST_USES_FR.description, ip.EVENT_LIST_USES_FR.value); 
		
		
		
		resetLabelsWidth();
		
		//simulation params
		componentsStack.add(STOP_TIME);
		componentsStack.add(ENSAMBLE_SIZE);
		componentsStack.add(RANDOM_SEED);
		componentsStack.add(COMPUTED_AFFINITY_PRECISION);
		componentsStack.add(DNA_SECTOR_SIZE);
		componentsStack.add(EVENT_LIST_SUBGROUP_SIZE);
		componentsStack.add(EVENT_LIST_USES_FR);
		
		

		
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
		
		if(ENSAMBLE_SIZE.getLabelWidth() > max){
			max = ENSAMBLE_SIZE.getLabelWidth();
		}	
		
	
		if(COMPUTED_AFFINITY_PRECISION.getLabelWidth() > max){
			max = COMPUTED_AFFINITY_PRECISION.getLabelWidth();
		}
		
		if(DNA_SECTOR_SIZE.getLabelWidth() > max){
			max = DNA_SECTOR_SIZE.getLabelWidth();
		}
		
		if(EVENT_LIST_SUBGROUP_SIZE.getLabelWidth() > max){
			max = EVENT_LIST_SUBGROUP_SIZE.getLabelWidth();
		}
		
		
		
		
		
		
		
		
		
		//SIMULATION PARAMATERS
		STOP_TIME.setLabelWidth(max);
		ENSAMBLE_SIZE.setLabelWidth(max);	
		RANDOM_SEED.setLabelWidth(max);	
		COMPUTED_AFFINITY_PRECISION.setLabelWidth(max);	
		DNA_SECTOR_SIZE.setLabelWidth(max);	
		EVENT_LIST_SUBGROUP_SIZE.setLabelWidth(max);	

	}
	
	
}

