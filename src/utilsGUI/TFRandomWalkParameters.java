package utilsGUI;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;

import objects.InputParameters;

/**
 * tabbed panel with random walk parameters 
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFRandomWalkParameters extends JPanel{
	private static final long serialVersionUID = 1L;

	private JPanel componentsStack;
	
	//TF RANDOM WALK PARAMATERS
	public LabelledCheckBox CHECK_OCCUPANCY_ON_BINDING;
	public LabelledCheckBox CHECK_OCCUPANCY_ON_SLIDING;
	public LabelledCheckBox CHECK_OCCUPANCY_ON_REBINDING;
	
	//TF RANDOM WALK RANDOM PARAMATERS
	public LabelledDouble TF_SLIDE_LEFT_PROBABILITY;
	public LabelledDouble TF_SLIDE_RIGHT_PROBABILITY;
	public LabelledDouble TF_UNBINDING_PROBABILITY;
	public LabelledDouble TF_JUMPING_PROBABILITY;
	public LabelledDouble TF_SPECIFIC_WAITING_TIME;
	public LabelledDouble TF_COLLISION_UNBIND_PROBABILITY;
	public LabelledDouble TF_AFFINITY_LANDSCAPE_ROUGHNESS;
	public LabelledDouble TF_HOP_STD_DISPLACEMENT;

	public LabelledInteger TF_STEP_LEFT_SIZE;
	public LabelledInteger TF_STEP_RIGHT_SIZE;
	public LabelledInteger TF_UNCORRELATED_DISPLACEMENT_SIZE;
	
	public LabelledCheckBox TF_STALLS_IF_BLOCKED;

	public LabelledCheckBox TF_IS_IMMOBILE;

	
	public LabelledCheckBox IS_BIASED_RANDOM_WALK;
	public LabelledCheckBox IS_TWO_STATE_RANDOM_WALK;

	
	public TFRandomWalkParameters(InputParameters ip){
		this.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH,GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT));
		this.setLayout(new FlowLayout());
		componentsStack = new JPanel(new GridLayout(0,1, GUIconstants.GRID_HGAP, GUIconstants.GRID_WGAP));
		componentsStack.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT,GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH));
		
		
		
		JLabel label1,label2;
		label1 = new JLabel(GUIconstants.SIMULATION_AREA_TF_RANDOM_WALK_GENERAL_PARAMATERS);
		label2 = new JLabel(GUIconstants.SIMULATION_AREA_TF_RANDOM_WALK_RANDOM_PARAMATERS);
		label1.setForeground(Color.LIGHT_GRAY);
		label2.setForeground(Color.LIGHT_GRAY);


		
		//TF RANDOM WALK PARAMATERS
		CHECK_OCCUPANCY_ON_BINDING = new LabelledCheckBox(ip.CHECK_OCCUPANCY_ON_BINDING.label, ip.CHECK_OCCUPANCY_ON_BINDING.description, ip.CHECK_OCCUPANCY_ON_BINDING.value); 
		CHECK_OCCUPANCY_ON_SLIDING = new LabelledCheckBox(ip.CHECK_OCCUPANCY_ON_SLIDING.label, ip.CHECK_OCCUPANCY_ON_SLIDING.description, ip.CHECK_OCCUPANCY_ON_SLIDING.value); 
		CHECK_OCCUPANCY_ON_REBINDING = new LabelledCheckBox(ip.CHECK_OCCUPANCY_ON_REBINDING.label, ip.CHECK_OCCUPANCY_ON_REBINDING.description, ip.CHECK_OCCUPANCY_ON_REBINDING.value); 

		//TF RANDOM WALK RANDOM PARAMATERS
		TF_SLIDE_LEFT_PROBABILITY = new LabelledDouble(ip.TF_SLIDE_LEFT_PROBABILITY.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_SLIDE_LEFT_PROBABILITY.description,ip.TF_SLIDE_LEFT_PROBABILITY.value);	
		TF_SLIDE_RIGHT_PROBABILITY = new LabelledDouble(ip.TF_SLIDE_RIGHT_PROBABILITY.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_SLIDE_RIGHT_PROBABILITY.description,ip.TF_SLIDE_RIGHT_PROBABILITY.value);	
		TF_UNBINDING_PROBABILITY = new LabelledDouble(ip.TF_UNBINDING_PROBABILITY.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_UNBINDING_PROBABILITY.description,ip.TF_UNBINDING_PROBABILITY.value);	
		TF_JUMPING_PROBABILITY = new LabelledDouble(ip.TF_JUMPING_PROBABILITY.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_JUMPING_PROBABILITY.description,ip.TF_JUMPING_PROBABILITY.value);	
		TF_SPECIFIC_WAITING_TIME = new LabelledDouble(ip.TF_SPECIFIC_WAITING_TIME.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_SPECIFIC_WAITING_TIME.description,ip.TF_SPECIFIC_WAITING_TIME.value);	
		TF_COLLISION_UNBIND_PROBABILITY = new LabelledDouble(ip.TF_COLLISION_UNBIND_PROBABILITY.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_COLLISION_UNBIND_PROBABILITY.description,ip.TF_COLLISION_UNBIND_PROBABILITY.value);	
		TF_AFFINITY_LANDSCAPE_ROUGHNESS = new LabelledDouble(ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.description,ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value);	
		TF_HOP_STD_DISPLACEMENT = new LabelledDouble(ip.TF_HOP_STD_DISPLACEMENT.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_HOP_STD_DISPLACEMENT.description,ip.TF_HOP_STD_DISPLACEMENT.value);	

		TF_STEP_LEFT_SIZE = new LabelledInteger(ip.TF_STEP_LEFT_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_STEP_LEFT_SIZE.description,ip.TF_STEP_LEFT_SIZE.value.intValue());
		TF_STEP_RIGHT_SIZE = new LabelledInteger(ip.TF_STEP_RIGHT_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_STEP_RIGHT_SIZE.description,ip.TF_STEP_RIGHT_SIZE.value.intValue());
		TF_UNCORRELATED_DISPLACEMENT_SIZE = new LabelledInteger(ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.description,ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.value.intValue());

		TF_STALLS_IF_BLOCKED = new LabelledCheckBox(ip.TF_STALLS_IF_BLOCKED.label, ip.TF_STALLS_IF_BLOCKED.description, ip.TF_STALLS_IF_BLOCKED.value); 
		TF_IS_IMMOBILE = new LabelledCheckBox(ip.TF_IS_IMMOBILE.label, ip.TF_IS_IMMOBILE.description, ip.TF_IS_IMMOBILE.value); 

		IS_BIASED_RANDOM_WALK = new LabelledCheckBox(ip.IS_BIASED_RANDOM_WALK.label, ip.IS_BIASED_RANDOM_WALK.description, ip.IS_BIASED_RANDOM_WALK.value); 
		IS_TWO_STATE_RANDOM_WALK = new LabelledCheckBox(ip.IS_TWO_STATE_RANDOM_WALK.label, ip.IS_TWO_STATE_RANDOM_WALK.description, ip.IS_TWO_STATE_RANDOM_WALK.value); 
		
		
		
		resetLabelsWidth();
		
		//TF RANDOM WALK PARAMATERS
		componentsStack.add(label1);
		componentsStack.add(CHECK_OCCUPANCY_ON_BINDING);
		componentsStack.add(CHECK_OCCUPANCY_ON_SLIDING);
		componentsStack.add(CHECK_OCCUPANCY_ON_REBINDING);
		
		//TF RANDOM WALK RANDOM PARAMATERS
		componentsStack.add(label2);
		componentsStack.add(TF_SLIDE_LEFT_PROBABILITY);
		componentsStack.add(TF_SLIDE_RIGHT_PROBABILITY);
		componentsStack.add(TF_UNBINDING_PROBABILITY);
		componentsStack.add(TF_JUMPING_PROBABILITY);
		componentsStack.add(TF_SPECIFIC_WAITING_TIME);
		componentsStack.add(TF_COLLISION_UNBIND_PROBABILITY);
		componentsStack.add(TF_AFFINITY_LANDSCAPE_ROUGHNESS);
		componentsStack.add(TF_HOP_STD_DISPLACEMENT);
		
		componentsStack.add(TF_STEP_LEFT_SIZE);
		componentsStack.add(TF_STEP_RIGHT_SIZE);
		componentsStack.add(TF_UNCORRELATED_DISPLACEMENT_SIZE);

		componentsStack.add(TF_STALLS_IF_BLOCKED);
		
		componentsStack.add(TF_IS_IMMOBILE);

		
		componentsStack.add(IS_BIASED_RANDOM_WALK);
		componentsStack.add(IS_TWO_STATE_RANDOM_WALK);

		
		
		
		
		this.add(componentsStack);
		
	}
	
	
	/**
	 * resets the labels width
	 */
	private void resetLabelsWidth(){
		

		
		
		//TF RANDOM WALK RANDOM PARAMATERS
		int max = TF_SLIDE_LEFT_PROBABILITY.getLabelWidth();
				
		if(TF_SLIDE_RIGHT_PROBABILITY.getLabelWidth() > max){
			max = TF_SLIDE_RIGHT_PROBABILITY.getLabelWidth();
		}
		
		if(TF_UNBINDING_PROBABILITY.getLabelWidth() > max){
			max = TF_UNBINDING_PROBABILITY.getLabelWidth();
		}
		
		if(TF_JUMPING_PROBABILITY.getLabelWidth() > max){
			max = TF_JUMPING_PROBABILITY.getLabelWidth();
		}
		
		if(TF_SPECIFIC_WAITING_TIME.getLabelWidth() > max){
			max = TF_SPECIFIC_WAITING_TIME.getLabelWidth();
		}
		
		
		if(TF_COLLISION_UNBIND_PROBABILITY.getLabelWidth() > max){
			max = TF_COLLISION_UNBIND_PROBABILITY.getLabelWidth();
		}
		
		if(TF_AFFINITY_LANDSCAPE_ROUGHNESS.getLabelWidth() > max){
			max = TF_AFFINITY_LANDSCAPE_ROUGHNESS.getLabelWidth();
		}	
		
		if(TF_HOP_STD_DISPLACEMENT.getLabelWidth() > max){
			max = TF_HOP_STD_DISPLACEMENT.getLabelWidth();
		}
		
		if(TF_STEP_LEFT_SIZE.getLabelWidth() > max){
			max = TF_STEP_LEFT_SIZE.getLabelWidth();
		}	
		
		if(TF_STEP_RIGHT_SIZE.getLabelWidth() > max){
			max = TF_STEP_RIGHT_SIZE.getLabelWidth();
		}	
		
		
		if(TF_UNCORRELATED_DISPLACEMENT_SIZE.getLabelWidth() > max){
			max = TF_UNCORRELATED_DISPLACEMENT_SIZE.getLabelWidth();
		}		
		
		
		//TF RANDOM WALK RANDOM PARAMATERS
		TF_SLIDE_LEFT_PROBABILITY.setLabelWidth(max);
		TF_SLIDE_RIGHT_PROBABILITY.setLabelWidth(max);	
		TF_UNBINDING_PROBABILITY.setLabelWidth(max);	
		TF_JUMPING_PROBABILITY.setLabelWidth(max);	
		TF_SPECIFIC_WAITING_TIME.setLabelWidth(max);	
		TF_COLLISION_UNBIND_PROBABILITY.setLabelWidth(max);
		TF_AFFINITY_LANDSCAPE_ROUGHNESS.setLabelWidth(max);	
		TF_HOP_STD_DISPLACEMENT.setLabelWidth(max);
		
		TF_STEP_LEFT_SIZE.setLabelWidth(max);
		TF_STEP_RIGHT_SIZE.setLabelWidth(max);
		TF_UNCORRELATED_DISPLACEMENT_SIZE.setLabelWidth(max);
	}
	
}
