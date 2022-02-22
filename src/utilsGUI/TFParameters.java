package utilsGUI;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;

import objects.InputParameters;

/**
 * tabbed panel with TF parameters
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFParameters  extends JPanel{

	

		private static final long serialVersionUID = 1L;

		private JPanel componentsStack;
	
		
		//TF PARAMETERS
		public LabelledFileChooser TF_FILE;
		public LabelledFileChooser TF_COOPERATIVITY_FILE;
		public LabelledFileChooser TS_FILE;
		
		//TF_RANDOM PARAMETERS
		public LabelledInteger TF_DBD_LENGTH_MIN;
		public LabelledInteger TF_DBD_LENGTH_MAX;
		public LabelledInteger TF_SPECIES_COUNT;
		public LabelledInteger 	TF_COPY_NUMBER_MIN;
		public LabelledInteger TF_COPY_NUMBER_MAX;
		public LabelledInteger TF_SIZE_LEFT;
		public LabelledInteger TF_SIZE_RIGHT;
		public LabelledDouble TF_ES;
		public LabelledDouble TF_ASSOC_RATE;
		public LabelledDouble  TF_PREBOUND_PROPORTION;
		public LabelledCheckBox  TF_PREBOUND_TO_HIGHEST_AFFINITY;
		public LabelledCheckBox TF_READ_IN_BOTH_DIRECTIONS;
		public LabelledCheckBox SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE;

		
		public TFParameters(InputParameters ip){
			this.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH,GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT));
			this.setLayout(new FlowLayout());
			componentsStack = new JPanel(new GridLayout(0,1, GUIconstants.GRID_HGAP, GUIconstants.GRID_WGAP));
			componentsStack.setMaximumSize(new Dimension(GUIconstants.SIMULATION_PARAMETERS_SIZE_HIGHT,GUIconstants.SIMULATION_PARAMETERS_SIZE_WIDTH));
			
			
			
			JLabel label1, label2, label3;
			label1 = new JLabel(GUIconstants.SIMULATION_AREA_TF_LOAD_PARAMATERS);
			label2 = new JLabel(GUIconstants.SIMULATION_AREA_TF_RANDOM_PARAMATERS);
			label3 = new JLabel(GUIconstants.SIMULATION_AREA_TF_GENERAL_PARAMATERS);
			label1.setForeground(Color.LIGHT_GRAY);
			label2.setForeground(Color.LIGHT_GRAY);
			label3.setForeground(Color.LIGHT_GRAY);
			

			//TF PARAMETERS
			TF_FILE = new LabelledFileChooser(ip.TF_FILE.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_FILE.description,ip.TF_FILE.value, true, true);
			TF_COOPERATIVITY_FILE = new LabelledFileChooser(ip.TF_COOPERATIVITY_FILE.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_COOPERATIVITY_FILE.description,ip.TF_COOPERATIVITY_FILE.value, true, true);
			TS_FILE = new LabelledFileChooser(ip.TS_FILE.label,GUIconstants.TEXTAREA_WIDTH,ip.TS_FILE.description,ip.TS_FILE.value, true, true);

			

			TF_READ_IN_BOTH_DIRECTIONS = new LabelledCheckBox(ip.TF_READ_IN_BOTH_DIRECTIONS.label, ip.TF_READ_IN_BOTH_DIRECTIONS.description, ip.TF_READ_IN_BOTH_DIRECTIONS.value); 
			SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE = new LabelledCheckBox(ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.label, ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.description, ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value); 

			//TF_RANDOM PARAMETERS
			TF_DBD_LENGTH_MIN = new LabelledInteger(ip.TF_DBD_LENGTH_MIN.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_DBD_LENGTH_MIN.description,ip.TF_DBD_LENGTH_MIN.value.intValue());
			TF_DBD_LENGTH_MAX = new LabelledInteger(ip.TF_DBD_LENGTH_MAX.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_DBD_LENGTH_MAX.description,ip.TF_DBD_LENGTH_MAX.value.intValue());
			TF_SPECIES_COUNT = new LabelledInteger(ip.TF_SPECIES_COUNT.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_SPECIES_COUNT.description,ip.TF_SPECIES_COUNT.value.intValue());
			TF_COPY_NUMBER_MIN = new LabelledInteger(ip.TF_COPY_NUMBER_MIN.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_COPY_NUMBER_MIN.description,ip.TF_COPY_NUMBER_MIN.value.intValue());
			TF_COPY_NUMBER_MAX = new LabelledInteger(ip.TF_COPY_NUMBER_MAX.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_COPY_NUMBER_MAX.description,ip.TF_COPY_NUMBER_MAX.value.intValue());
			TF_SIZE_LEFT = new LabelledInteger(ip.TF_SIZE_LEFT.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_SIZE_LEFT.description,ip.TF_SIZE_LEFT.value.intValue());
			TF_SIZE_RIGHT = new LabelledInteger(ip.TF_SIZE_RIGHT.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_SIZE_RIGHT.description,ip.TF_SIZE_RIGHT.value.intValue());
			TF_ES = new LabelledDouble(ip.TF_ES.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_ES.description,ip.TF_ES.value);	
			TF_ASSOC_RATE = new LabelledDouble(ip.TF_ASSOC_RATE.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_ASSOC_RATE.description,ip.TF_ASSOC_RATE.value);	
			TF_PREBOUND_PROPORTION = new LabelledDouble(ip.TF_PREBOUND_PROPORTION.label,GUIconstants.TEXTAREA_WIDTH,ip.TF_PREBOUND_PROPORTION.description,ip.TF_PREBOUND_PROPORTION.value);	
			TF_PREBOUND_TO_HIGHEST_AFFINITY = new LabelledCheckBox(ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.label, ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.description, ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.value); 

			
			
			
			resetLabelsWidth();
			
			
			
			
			
			
			

			
			
			
			
			//TF PARAMETERS
			componentsStack.add(label1);
			componentsStack.add(TF_FILE);
			componentsStack.add(TF_COOPERATIVITY_FILE);
			componentsStack.add(TS_FILE);


			componentsStack.add(label3);
			componentsStack.add(TF_READ_IN_BOTH_DIRECTIONS);
			componentsStack.add(SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE);
			
			//TF_RANDOM PARAMETERS
			componentsStack.add(label2);
			componentsStack.add(TF_SPECIES_COUNT);
			componentsStack.add(TF_COPY_NUMBER_MIN);
			componentsStack.add(TF_COPY_NUMBER_MAX);
			componentsStack.add(TF_DBD_LENGTH_MIN);
			componentsStack.add(TF_DBD_LENGTH_MAX);
			componentsStack.add(TF_SIZE_LEFT);
			componentsStack.add(TF_SIZE_RIGHT);
			componentsStack.add(TF_ES);
			componentsStack.add(TF_ASSOC_RATE);
			componentsStack.add(TF_PREBOUND_PROPORTION);
			componentsStack.add(TF_PREBOUND_TO_HIGHEST_AFFINITY);
			
			
			this.add(componentsStack);
			

		}
		
		
		/**
		 * resets the labels width
		 */
		private void resetLabelsWidth(){

			
			//SIMULATION PARAMATERS
			int max = TF_FILE.getLabelWidth();
					
			if(TF_COOPERATIVITY_FILE.getLabelWidth() > max){
				max = TF_COOPERATIVITY_FILE.getLabelWidth();
			}
			
			if(TS_FILE.getLabelWidth() > max){
				max = TS_FILE.getLabelWidth();
			}
			
			//TF_RANDOM PARAMETERS
			if(TF_DBD_LENGTH_MIN.getLabelWidth() > max){
				max = TF_DBD_LENGTH_MIN.getLabelWidth();
			}
			
			if(TF_DBD_LENGTH_MAX.getLabelWidth() > max){
				max = TF_DBD_LENGTH_MAX.getLabelWidth();
			}
			
			if(TF_SPECIES_COUNT.getLabelWidth() > max){
				max = TF_SPECIES_COUNT.getLabelWidth();
			}			
			if(TF_COPY_NUMBER_MIN.getLabelWidth() > max){
				max = TF_COPY_NUMBER_MIN.getLabelWidth();
			}
			
			if(TF_COPY_NUMBER_MAX.getLabelWidth() > max){
				max = TF_COPY_NUMBER_MAX.getLabelWidth();
			}	
			
			if(TF_SIZE_LEFT.getLabelWidth() > max){
				max = TF_SIZE_LEFT.getLabelWidth();
			}
			
			if(TF_SIZE_RIGHT.getLabelWidth() > max){
				max = TF_SIZE_RIGHT.getLabelWidth();
			}	
		
			if(TF_ES.getLabelWidth() > max){
				max = TF_ES.getLabelWidth();
			}	

			
			if(TF_ASSOC_RATE.getLabelWidth() > max){
				max = TF_ASSOC_RATE.getLabelWidth();
			}	
			
			
			if(TF_PREBOUND_PROPORTION.getLabelWidth() > max){
				max = TF_PREBOUND_PROPORTION.getLabelWidth();
			}	
			

			
			//TF PARAMETERS
			TF_FILE.setLabelWidth(max);
			TF_COOPERATIVITY_FILE.setLabelWidth(max);	
			TS_FILE.setLabelWidth(max);

			//TF_RANDOM PARAMETERS
			TF_DBD_LENGTH_MIN.setLabelWidth(max);	
			TF_DBD_LENGTH_MAX.setLabelWidth(max);	
			TF_SPECIES_COUNT.setLabelWidth(max);	
			TF_COPY_NUMBER_MIN.setLabelWidth(max);
			TF_COPY_NUMBER_MAX.setLabelWidth(max);	
			TF_SIZE_LEFT.setLabelWidth(max);
			TF_SIZE_RIGHT.setLabelWidth(max);
			TF_ES.setLabelWidth(max);
			TF_ASSOC_RATE.setLabelWidth(max);
			TF_PREBOUND_PROPORTION.setLabelWidth(max);

		}
		
		
		
		
		
}
