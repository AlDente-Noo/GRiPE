package simulator;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.filechooser.FileFilter;


import objects.InputParameters;
import utils.Utils;
import utilsGUI.*;


/**
 * runnable class with GUI
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class SimulatorGUI {
	private JFrame frame;
	private JPanel setupArea;
	private JScrollPane setupAreaScroll;

	private JTabbedPane setupAreaTabbedPane;
	private JPanel simulationArea;
	private JButton simulationsStart;
	private JButton simulationsLoad;
	private JProgressBar simulationsProgress; 
	private JLabel timeLabel;
	private SimulatorThread simulatorThread;
	private boolean loaded;
	private boolean simulating;
	private boolean paused;
	private int steps;
	
	private JTextArea statusTextArea;
	private JScrollPane statusAreaScroll;

	private SimulationParameters simulationParamaters;
	private OutputParameters outputParamaters;
	private TFParameters TFParameters;
	private DNAParameters DNAParameters;
	private TFRandomWalkParameters TFRandomWalkParameters;


	
	
	private InputParameters ip;
	
	private String currentFile;
	private JFileChooser fc;
	private ImageIcon buttonPlayImg;
	private ImageIcon buttonPauseImg;

	private ImageIcon buttonLoadImg;
	private ImageIcon buttonUnloadImg;
	
	private ImageIcon logoImg;

	DateFormat df;
	
	/**
	 * class constructor
	 */
	public SimulatorGUI(String filename){
		ip = new InputParameters(filename);
		currentFile = filename;
		
		df = new SimpleDateFormat("HHH':'mm':'ss");
		df.setTimeZone(TimeZone.getTimeZone("GMT+0"));
		
		simulatorThread=null;
		loaded=false;
		simulating = false;
		paused =false;
		steps = GUIconstants.SIMULATION_PROGRESS_MAX - GUIconstants.SIMULATION_PROGRESS_MIN;	
	}
	
	
	public static void main(String args[]) {
		
		SimulatorGUI gui;
		String filename="";
		
		if(args.length >0){
			filename = args[0];
		}
		
		gui = new SimulatorGUI(filename);
		
		gui.makeFrame();
	}
	
	
	
	

	/**
	 * construct the frame
	 */
	private void makeFrame(){
		frame = new JFrame(GUIconstants.WINDOW_TITLE);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		
		Container contentPane = frame.getContentPane();
		contentPane.setLayout(new BorderLayout());
		
		fc = new JFileChooser();
		FileFilter filter1 = new ExtensionFileFilter(GUIconstants.EXTENSION_FILE_DESCRIPTION, new String[] { GUIconstants.EXTENSION_FILE});
	    fc.setFileFilter(filter1);
	    
	    
		InputStream in;
		logoImg = new ImageIcon();
	    File logoImgFile = new File(GUIconstants.LOGO_IMAGE);
		if(!logoImgFile.exists()){
			in = this.getClass().getClassLoader().getResourceAsStream(GUIconstants.LOGO_IMAGE);
			try {
				if(in.available() >0){
					logoImg = new ImageIcon(ImageIO.read(in));
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else{
			logoImg= new ImageIcon(GUIconstants.LOGO_IMAGE);
		}
		frame.setIconImage(logoImg.getImage());
		
		makeSimulationArea(contentPane);
		makeSetupArea(contentPane);
		makeStatusArea(contentPane);
		statusAreaScroll.setVisible(false);

		makeMenuBar(frame);
		
		frame.pack();
		frame.setLocationRelativeTo(null);
		frame.setResizable(true);
		frame.setVisible(true);
	}
	
	
	public void makeSimulationArea(Container contentPane){
		//simulation area
		simulationArea = new JPanel();
		simulationArea.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED));
		simulationArea.setLayout(new FlowLayout(FlowLayout.LEFT));
		simulationArea.setPreferredSize(new Dimension(GUIconstants.SETUP_AREA_WIDTH,GUIconstants.SIMULATION_AREA_HIGHT));

		// load button

		
		InputStream in;
		buttonLoadImg = new ImageIcon();
	    File buttonLoadImgFile = new File(GUIconstants.LOAD_IMAGE);
		if(!buttonLoadImgFile.exists()){
			in = this.getClass().getClassLoader().getResourceAsStream(GUIconstants.LOAD_IMAGE);
			try {
				if(in.available() >0){
					buttonLoadImg = new ImageIcon(ImageIO.read(in));
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else{
			buttonLoadImg= new ImageIcon(GUIconstants.LOAD_IMAGE);
		}
	
		buttonUnloadImg = new ImageIcon();		
	    File buttonUnloadImgFile = new File(GUIconstants.UNLOAD_IMAGE);
		if(!buttonUnloadImgFile.exists()){
			in = this.getClass().getClassLoader().getResourceAsStream(GUIconstants.UNLOAD_IMAGE);
			try {
				if(in.available() >0){
					buttonUnloadImg = new ImageIcon(ImageIO.read(in));
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else{
			buttonUnloadImg= new ImageIcon(GUIconstants.LOAD_IMAGE);
		}
		
		simulationsLoad = new JButton(buttonLoadImg);
		simulationsLoad.setText(GUIconstants.LOAD_BUTTON);
		simulationsLoad.setBackground(Color.white);
		loaded = false;
		simulationsLoad.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				loadSimulation();
			}
		});
		simulationArea.add(simulationsLoad);
		
		//Simulation start
		buttonPlayImg = new ImageIcon();
	    File buttonPlayImgFile = new File(GUIconstants.PLAY_IMAGE);
		if(!buttonPlayImgFile.exists()){
			in = this.getClass().getClassLoader().getResourceAsStream(GUIconstants.PLAY_IMAGE);
			try {
				if(in.available() >0){
					buttonPlayImg = new ImageIcon(ImageIO.read(in));
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else{
			buttonPlayImg= new ImageIcon(GUIconstants.PLAY_IMAGE);
		}
		
		buttonPauseImg = new ImageIcon();
	    File buttonPauseImgFile = new File(GUIconstants.PAUSE_IMAGE);
		if(!buttonPauseImgFile.exists()){
			in = this.getClass().getClassLoader().getResourceAsStream(GUIconstants.PAUSE_IMAGE);
			try {
				if(in.available() >0){
					buttonPauseImg = new ImageIcon(ImageIO.read(in));
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else{
			buttonPauseImg= new ImageIcon(GUIconstants.PAUSE_IMAGE);
		}
		

		
		simulationsStart = new JButton(buttonPlayImg);
		simulationsStart.setText(GUIconstants.SIMULATE_BUTTON_START);
		simulationsStart.setBackground(Color.white);
		simulationsStart.setEnabled(false);
		simulationsStart.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				startSimulation();
			}
		});
		simulating=false;
		simulationArea.add(simulationsStart);
		

		
		
		// progress bar
		simulationsProgress = new JProgressBar();
		simulationsProgress.setMinimum(GUIconstants.SIMULATION_PROGRESS_MIN);
		simulationsProgress.setMaximum(GUIconstants.SIMULATION_PROGRESS_MAX);
		simulationsProgress.setPreferredSize(new Dimension(GUIconstants.SIMULATION_PROGRESS_SIZE_WIDTH,GUIconstants.SIMULATION_PROGRESS_SIZE_HIGHT));
		simulationsProgress.setEnabled(false);
		simulationsProgress.setStringPainted(true);
		simulationsProgress.setBackground(Color.white);
		simulationsProgress.setIndeterminate(false);
		simulationArea.add(simulationsProgress);
		
		//Time label
		timeLabel = new JLabel(getTimeString(0,0)+" ");
        timeLabel.setEnabled(false);
		simulationArea.add(timeLabel);
		
		
		contentPane.add(simulationArea, BorderLayout.SOUTH);

	}
	

	public void makeSetupArea(Container contentPane){
		//setup area
		setupArea = new JPanel();		
	
		setupAreaTabbedPane = new JTabbedPane();
		setupAreaTabbedPane.setPreferredSize(new Dimension(GUIconstants.SETUP_AREA_WIDTH-GUIconstants.SCROLLBAR_SIZE,GUIconstants.SETUP_AREA_HIGHT-GUIconstants.SCROLLBAR_SIZE));

		setupArea.setLayout(new FlowLayout());
		


		simulationParamaters = new SimulationParameters(ip);
		setupAreaTabbedPane.addTab(GUIconstants.SIMULATION_AREA_SIMULATION_PARAMATERS, null, simulationParamaters, null);
		
		outputParamaters = new OutputParameters(ip);
		setupAreaTabbedPane.addTab(GUIconstants.SIMULATION_AREA_OUTPUT_PARAMATERS, null, outputParamaters, null);
		
		TFParameters = new TFParameters(ip);
		setupAreaTabbedPane.addTab(GUIconstants.SIMULATION_AREA_TF_PARAMATERS, null, TFParameters, null);

		
		DNAParameters = new DNAParameters(ip);
		setupAreaTabbedPane.addTab(GUIconstants.SIMULATION_AREA_DNA_PARAMATERS, null, DNAParameters, null);

		
		TFRandomWalkParameters = new TFRandomWalkParameters(ip);
		setupAreaTabbedPane.addTab(GUIconstants.SIMULATION_AREA_TF_RANDOM_WALK_PARAMATERS, null, TFRandomWalkParameters, null);

		
		setupArea.add(setupAreaTabbedPane);
		
		setupAreaScroll = new JScrollPane(setupArea);
		setupAreaScroll.setPreferredSize(new Dimension(GUIconstants.SETUP_AREA_WIDTH,GUIconstants.SETUP_AREA_HIGHT));
		contentPane.add(setupAreaScroll,BorderLayout.CENTER);
		
		
		//contentPane.add(setupArea,BorderLayout.CENTER);
	}
	

	/**
	 * constructs the status area
	 * @param contentPane
	 */
	public void makeStatusArea(Container contentPane){	
		statusTextArea= new AutoScrollingTextArea();
		statusTextArea.setSize(new Dimension(GUIconstants.SETUP_AREA_WIDTH,GUIconstants.SETUP_AREA_HIGHT));
		statusTextArea.setEditable(false);
		
		statusAreaScroll = new JScrollPane(statusTextArea);
		statusAreaScroll.setPreferredSize(new Dimension(GUIconstants.SETUP_AREA_WIDTH,GUIconstants.SETUP_AREA_HIGHT));
		contentPane.add(statusAreaScroll,BorderLayout.NORTH);
	}
	
	/**
	 * create the menu bar
	 * @param frame
	 */
	public void makeMenuBar(JFrame frame){
		JMenuBar menubar = new JMenuBar();
		frame.setJMenuBar(menubar);
		
		//create file menu
		JMenu fileMenu = new JMenu(GUIconstants.MENU_FILE);
		menubar.add(fileMenu);
		
		JMenuItem openItem = new JMenuItem(GUIconstants.MENU_FILE_OPEN);
		openItem.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				openFile();
			}
		});
		fileMenu.add(openItem);
		
		
		JMenuItem saveItem = new JMenuItem(GUIconstants.MENU_FILE_SAVE);
		saveItem.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				saveFile(e);
			}
		});
		fileMenu.add(saveItem);
		
		
		JMenuItem saveasItem = new JMenuItem(GUIconstants.MENU_FILE_SAVE_AS);
		saveasItem.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				saveFile(e);
			}
		});
		fileMenu.add(saveasItem);
		
		JMenuItem closeItem = new JMenuItem(GUIconstants.MENU_FILE_CLOSE);
		closeItem.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				closeFile();
			}
		});
		fileMenu.add(closeItem);
		
		//create help menu
		JMenu helpMenu = new JMenu(GUIconstants.MENU_HELP);
		menubar.add(helpMenu);
		
		JMenuItem aboutItem = new JMenuItem(GUIconstants.MENU_HELP_ABOUT);
		aboutItem.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				aboutHelp();
			}
		});
		helpMenu.add(aboutItem);
		
		
	}
	
	

	/**
	 * writes a line in the status area
	 * @param line
	 */
	public void printlnStatusArea(String line){
		statusTextArea.append(line+"\n");
	}
	
	

	/**
	 * opens a file
	 */
	private void openFile(){
		
		int returnVal = fc.showOpenDialog(setupArea);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            currentFile = file.getAbsolutePath();
            BufferedReader br;
			try {
				br = new BufferedReader(new FileReader(currentFile));
	            ip.loadParameters(br);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
            this.setInputParameters();
		}
		
		
	}

	
	/**
	 * save a file
	 */
	private boolean saveFile(ActionEvent e){
		boolean result=false;
		
		if(e!=null && e.getActionCommand().equalsIgnoreCase(GUIconstants.MENU_FILE_SAVE_AS)){
			currentFile="";
		}
		
		if(currentFile==""){
			int returnVal = fc.showSaveDialog(setupArea);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
	            File file = fc.getSelectedFile();
	            currentFile = file.getAbsolutePath();
	            if(!currentFile.endsWith(GUIconstants.EXTENSION_FILE)){
	            	currentFile+="."+GUIconstants.EXTENSION_FILE;
	            }
			}
		}
		
		if(currentFile!=""){
			//This is where a real application would save the file.
			result = true;
			saveModel(currentFile);
		}
		
		return result;

	}
	
	/**
	 * close a file
	 */
	private void closeFile(){
		System.exit(0);     
	}

	
	
	/**
	 * save a file
	 */
	private void aboutHelp(){
		  JOptionPane.showMessageDialog(frame, GUIconstants.MENU_ABOUT_DESCRIPTION, GUIconstants.MENU_HELP_ABOUT, JOptionPane.INFORMATION_MESSAGE, logoImg);
	}
	
	
	
	/**
	 * save the model
	 */
	private void saveModel(String currentFile){
		getInputParameters();
		this.currentFile = this.ip.exportParameterFile(currentFile).getAbsolutePath();
	}
	

	/**
	 * collect the new input parameters
	 */
	private void getInputParameters(){
		//SIMULATION PARAMATERS
		ip.STOP_TIME.value= simulationParamaters.STOP_TIME.getValue();	
		ip.ENSEMBLE_SIZE.value= simulationParamaters.ENSAMBLE_SIZE.getValue();	
		ip.RANDOM_SEED.value= simulationParamaters.RANDOM_SEED.getValue();	
		ip.COMPUTED_AFFINITY_PRECISION.value= simulationParamaters.COMPUTED_AFFINITY_PRECISION.getValue();	
		ip.DNA_SECTOR_SIZE.value= simulationParamaters.DNA_SECTOR_SIZE.getValue();	
		ip.EVENT_LIST_SUBGROUP_SIZE.value= simulationParamaters.EVENT_LIST_SUBGROUP_SIZE.getValue();	
		ip.EVENT_LIST_USES_FR.value= simulationParamaters.EVENT_LIST_USES_FR.getValue();
		
		//SIMULATION-OUTPUT PARAMATERS
		ip.OUTPUT_FOLDER.value= outputParamaters.OUTPUT_FOLDER.getValue();	
		ip.OUTPUT_FILENAME.value= outputParamaters.OUTPUT_FILENAME.getValue();	
		ip.PRINT_INTERMEDIARY_RESULTS_AFTER.value= outputParamaters.PRINT_INTERMEDIARY_RESULTS_AFTER.getValue();	
		ip.PRINT_FINAL_OCCUPANCY.value= outputParamaters.PRINT_FINAL_OCCUPANCY.getValue();	
		ip.DEBUG_MODE.value= outputParamaters.DEBUG_MODE.getValue();		
		ip.OUTPUT_TF.value= outputParamaters.OUTPUT_TF.getValue();	
		ip.OUTPUT_TF_POINTS.value= outputParamaters.OUTPUT_TF_POINTS.getValue();	
		ip.FOLLOW_TS.value= outputParamaters.FOLLOW_TS.getValue();
		ip.OUTPUT_AFFINITY_LANDSCAPE.value= outputParamaters.OUTPUT_AFFINITY_LANDSCAPE.getValue();
		ip.OUTPUT_BINDING_ENERGY.value= outputParamaters.OUTPUT_BINDING_ENERGY.getValue();
		ip.OUTPUT_DNA_OCCUPANCY.value= outputParamaters.OUTPUT_DNA_OCCUPANCY.getValue();	
		ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value= outputParamaters.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.getValue();	
		ip.OUTPUT_SLIDING_LENGTHS.value= outputParamaters.OUTPUT_SLIDING_LENGTHS.getValue();	
		ip.WIG_STEP.value= outputParamaters.WIG_STEP.getValue();	
		ip.WIG_THRESHOLD.value= outputParamaters.WIG_THRESHOLD.getValue();	

		//TF PARAMETERS
		ip.TF_FILE.value= TFParameters.TF_FILE.getValue();	
		ip.TF_COOPERATIVITY_FILE.value= TFParameters.TF_COOPERATIVITY_FILE.getValue();	

		//TF_RANDOM PARAMETERS
		ip.TF_DBD_LENGTH_MIN.value= TFParameters.TF_DBD_LENGTH_MIN.getValue();	
		ip.TF_DBD_LENGTH_MAX.value= TFParameters.TF_DBD_LENGTH_MAX.getValue();	
		ip.TF_SPECIES_COUNT.value= TFParameters.TF_SPECIES_COUNT.getValue();	
		ip.TF_COPY_NUMBER_MIN.value= TFParameters.TF_COPY_NUMBER_MIN.getValue();	
		ip.TF_COPY_NUMBER_MAX.value= TFParameters.TF_COPY_NUMBER_MAX.getValue();	
		ip.TF_SIZE_LEFT.value= TFParameters.TF_SIZE_LEFT.getValue();	
		ip.TF_SIZE_RIGHT.value= TFParameters.TF_SIZE_RIGHT.getValue();	
		ip.TF_ES.value= TFParameters.TF_ES.getValue();	
		ip.TF_ASSOC_RATE.value= TFParameters.TF_ASSOC_RATE.getValue();	
		ip.TF_PREBOUND_PROPORTION.value= TFParameters.TF_PREBOUND_PROPORTION.getValue();	
		ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.value= TFParameters.TF_PREBOUND_TO_HIGHEST_AFFINITY.getValue();	

		ip.TF_READ_IN_BOTH_DIRECTIONS.value= TFParameters.TF_READ_IN_BOTH_DIRECTIONS.getValue();	
		ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value= TFParameters.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.getValue();	

	
		
		
		//DNA PARAMETERS
		ip.DNA_SEQUENCE_FILE.value= DNAParameters.DNA_SEQUENCE_FILE.getValue();	

		//DNA_RANDOM PARAMETERS
		ip.DNA_LENGTH.value= DNAParameters.DNA_LENGTH.getValue();	
		ip.DNA_PROPORTION_OF_A.value= DNAParameters.DNA_PROPORTION_OF_A.getValue();	
		ip.DNA_PROPORTION_OF_T.value= DNAParameters.DNA_PROPORTION_OF_T.getValue();	
		ip.DNA_PROPORTION_OF_C.value= DNAParameters.DNA_PROPORTION_OF_C.getValue();	
		ip.DNA_PROPORTION_OF_G.value= DNAParameters.DNA_PROPORTION_OF_G.getValue();	
		ip.DNA_BOUNDARY_CONDITION.value= DNAParameters.DNA_BOUNDARY_CONDITION.getValue();	
		

		//TF RANDOM WALK PARAMATERS
		ip.CHECK_OCCUPANCY_ON_BINDING.value= TFRandomWalkParameters.CHECK_OCCUPANCY_ON_BINDING.getValue();	
		ip.CHECK_OCCUPANCY_ON_SLIDING.value= TFRandomWalkParameters.CHECK_OCCUPANCY_ON_SLIDING.getValue();	
		ip.CHECK_OCCUPANCY_ON_REBINDING.value= TFRandomWalkParameters.CHECK_OCCUPANCY_ON_REBINDING.getValue();	
	
		//TF RANDOM WALK RANDOM PARAMATERS
		ip.TF_SLIDE_LEFT_PROBABILITY.value= TFRandomWalkParameters.TF_SLIDE_LEFT_PROBABILITY.getValue();	
		ip.TF_SLIDE_RIGHT_PROBABILITY.value= TFRandomWalkParameters.TF_SLIDE_RIGHT_PROBABILITY.getValue();	
		ip.TF_UNBINDING_PROBABILITY.value= TFRandomWalkParameters.TF_UNBINDING_PROBABILITY.getValue();	
		ip.TF_JUMPING_PROBABILITY.value= TFRandomWalkParameters.TF_JUMPING_PROBABILITY.getValue();	
		ip.TF_SPECIFIC_WAITING_TIME.value= TFRandomWalkParameters.TF_SPECIFIC_WAITING_TIME.getValue();	
		ip.TF_COLLISION_UNBIND_PROBABILITY.value= TFRandomWalkParameters.TF_COLLISION_UNBIND_PROBABILITY.getValue();	
		ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value= TFRandomWalkParameters.TF_AFFINITY_LANDSCAPE_ROUGHNESS.getValue();	
		ip.TF_HOP_STD_DISPLACEMENT.value= TFRandomWalkParameters.TF_HOP_STD_DISPLACEMENT.getValue();	
		ip.TF_STEP_LEFT_SIZE.value= TFRandomWalkParameters.TF_STEP_LEFT_SIZE.getValue();	
		ip.TF_STEP_RIGHT_SIZE.value= TFRandomWalkParameters.TF_STEP_RIGHT_SIZE.getValue();	
		ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.value= TFRandomWalkParameters.TF_UNCORRELATED_DISPLACEMENT_SIZE.getValue();	
		ip.TF_STALLS_IF_BLOCKED.value= TFRandomWalkParameters.TF_STALLS_IF_BLOCKED.getValue();	
		ip.TF_IS_IMMOBILE.value= TFRandomWalkParameters.TF_IS_IMMOBILE.getValue();	
		ip.IS_BIASED_RANDOM_WALK.value= TFRandomWalkParameters.IS_BIASED_RANDOM_WALK.getValue();	
		ip.IS_TWO_STATE_RANDOM_WALK.value= TFRandomWalkParameters.IS_TWO_STATE_RANDOM_WALK.getValue();	

				
	}
	
	/**
	 * sets the new input parameters
	 */
	private void setInputParameters(){
		
		//SIMULATION PARAMATERS
		simulationParamaters.STOP_TIME.setValue(ip.STOP_TIME.value);
		simulationParamaters.ENSAMBLE_SIZE.setValue(ip.ENSEMBLE_SIZE.value);
		simulationParamaters.RANDOM_SEED.setValue(ip.RANDOM_SEED.value);
		simulationParamaters.COMPUTED_AFFINITY_PRECISION.setValue(ip.COMPUTED_AFFINITY_PRECISION.value);
		simulationParamaters.DNA_SECTOR_SIZE.setValue(ip.DNA_SECTOR_SIZE.value);
		simulationParamaters.EVENT_LIST_SUBGROUP_SIZE.setValue(ip.EVENT_LIST_SUBGROUP_SIZE.value);
		simulationParamaters.EVENT_LIST_USES_FR.setValue(ip.EVENT_LIST_USES_FR.value);
		
		//SIMULATION-OUTPUT PARAMATERS
		outputParamaters.OUTPUT_FOLDER.setValue(ip.OUTPUT_FOLDER.value);
		outputParamaters.OUTPUT_FILENAME.setValue(ip.OUTPUT_FILENAME.value);
		outputParamaters.PRINT_INTERMEDIARY_RESULTS_AFTER.setValue(ip.PRINT_INTERMEDIARY_RESULTS_AFTER.value);
		outputParamaters.PRINT_FINAL_OCCUPANCY.setValue(ip.PRINT_FINAL_OCCUPANCY.value);
		outputParamaters.DEBUG_MODE.setValue(ip.DEBUG_MODE.value);
		outputParamaters.OUTPUT_TF.setValue(ip.OUTPUT_TF.value);
		outputParamaters.OUTPUT_TF_POINTS.setValue(ip.OUTPUT_TF_POINTS.value);
		outputParamaters.FOLLOW_TS.setValue(ip.FOLLOW_TS.value);
		outputParamaters.OUTPUT_AFFINITY_LANDSCAPE.setValue(ip.OUTPUT_AFFINITY_LANDSCAPE.value);
		outputParamaters.OUTPUT_BINDING_ENERGY.setValue(ip.OUTPUT_BINDING_ENERGY.value);
		outputParamaters.OUTPUT_DNA_OCCUPANCY.setValue(ip.OUTPUT_DNA_OCCUPANCY.value);
		outputParamaters.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.setValue(ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value);
		outputParamaters.OUTPUT_SLIDING_LENGTHS.setValue(ip.OUTPUT_SLIDING_LENGTHS.value);
		outputParamaters.WIG_STEP.setValue(ip.WIG_STEP.value);
		outputParamaters.WIG_THRESHOLD.setValue(ip.WIG_THRESHOLD.value);
		
		
		
		//TF PARAMETERS
		TFParameters.TF_FILE.setValue(ip.TF_FILE.value);
		TFParameters.TF_COOPERATIVITY_FILE.setValue(ip.TF_COOPERATIVITY_FILE.value);

		//TF_RANDOM PARAMETERS
		TFParameters.TF_DBD_LENGTH_MIN.setValue(ip.TF_DBD_LENGTH_MIN.value);
		TFParameters.TF_DBD_LENGTH_MAX.setValue(ip.TF_DBD_LENGTH_MAX.value);
		TFParameters.TF_SPECIES_COUNT.setValue(ip.TF_SPECIES_COUNT.value);
		TFParameters.TF_COPY_NUMBER_MIN.setValue(ip.TF_COPY_NUMBER_MIN.value);
		TFParameters.TF_COPY_NUMBER_MAX.setValue(ip.TF_COPY_NUMBER_MAX.value);
		TFParameters.TF_SIZE_LEFT.setValue(ip.TF_SIZE_LEFT.value);
		TFParameters.TF_SIZE_RIGHT.setValue(ip.TF_SIZE_RIGHT.value);
		TFParameters.TF_ES.setValue(ip.TF_ES.value);
		TFParameters.TF_ASSOC_RATE.setValue(ip.TF_ASSOC_RATE.value);
		TFParameters.TF_PREBOUND_PROPORTION.setValue(ip.TF_PREBOUND_PROPORTION.value);
		TFParameters.TF_PREBOUND_TO_HIGHEST_AFFINITY.setValue(ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.value);
		TFParameters.TF_READ_IN_BOTH_DIRECTIONS.setValue(ip.TF_READ_IN_BOTH_DIRECTIONS.value);
		TFParameters.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.setValue(ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value);
		
		

		
		
		//DNA PARAMETERS
		DNAParameters.DNA_SEQUENCE_FILE.setValue(ip.DNA_SEQUENCE_FILE.value);
		
		//DNA_RANDOM PARAMETERS
		DNAParameters.DNA_LENGTH.setValue(ip.DNA_LENGTH.value);
		DNAParameters.DNA_PROPORTION_OF_A.setValue(ip.DNA_PROPORTION_OF_A.value);
		DNAParameters.DNA_PROPORTION_OF_T.setValue(ip.DNA_PROPORTION_OF_T.value);
		DNAParameters.DNA_PROPORTION_OF_C.setValue(ip.DNA_PROPORTION_OF_C.value);
		DNAParameters.DNA_PROPORTION_OF_G.setValue(ip.DNA_PROPORTION_OF_G.value);
		DNAParameters.DNA_BOUNDARY_CONDITION.setValue(ip.DNA_BOUNDARY_CONDITION.value);
		
		//TF RANDOM WALK PARAMATERS
		TFRandomWalkParameters.CHECK_OCCUPANCY_ON_BINDING.setValue(ip.CHECK_OCCUPANCY_ON_BINDING.value);
		TFRandomWalkParameters.CHECK_OCCUPANCY_ON_SLIDING.setValue(ip.CHECK_OCCUPANCY_ON_SLIDING.value);
		TFRandomWalkParameters.CHECK_OCCUPANCY_ON_REBINDING.setValue(ip.CHECK_OCCUPANCY_ON_REBINDING.value);
		
		
		//TF RANDOM WALK RANDOM PARAMATERS
		TFRandomWalkParameters.TF_SLIDE_LEFT_PROBABILITY.setValue(ip.TF_SLIDE_LEFT_PROBABILITY.value);
		TFRandomWalkParameters.TF_SLIDE_RIGHT_PROBABILITY.setValue(ip.TF_SLIDE_RIGHT_PROBABILITY.value);
		TFRandomWalkParameters.TF_UNBINDING_PROBABILITY.setValue(ip.TF_UNBINDING_PROBABILITY.value);
		TFRandomWalkParameters.TF_JUMPING_PROBABILITY.setValue(ip.TF_JUMPING_PROBABILITY.value);
		TFRandomWalkParameters.TF_SPECIFIC_WAITING_TIME.setValue(ip.TF_SPECIFIC_WAITING_TIME.value);
		TFRandomWalkParameters.TF_COLLISION_UNBIND_PROBABILITY.setValue(ip.TF_COLLISION_UNBIND_PROBABILITY.value);
		TFRandomWalkParameters.TF_AFFINITY_LANDSCAPE_ROUGHNESS.setValue(ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value);
		TFRandomWalkParameters.TF_HOP_STD_DISPLACEMENT.setValue(ip.TF_HOP_STD_DISPLACEMENT.value);

		TFRandomWalkParameters.TF_STEP_LEFT_SIZE.setValue(ip.TF_STEP_LEFT_SIZE.value);
		TFRandomWalkParameters.TF_STEP_RIGHT_SIZE.setValue(ip.TF_STEP_RIGHT_SIZE.value);
		TFRandomWalkParameters.TF_UNCORRELATED_DISPLACEMENT_SIZE.setValue(ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.value);
		
		TFRandomWalkParameters.TF_STALLS_IF_BLOCKED.setValue(ip.TF_STALLS_IF_BLOCKED.value);
		TFRandomWalkParameters.TF_IS_IMMOBILE.setValue(ip.TF_IS_IMMOBILE.value);
		TFRandomWalkParameters.IS_BIASED_RANDOM_WALK.setValue(ip.IS_BIASED_RANDOM_WALK.value);
		TFRandomWalkParameters.IS_TWO_STATE_RANDOM_WALK.setValue(ip.IS_TWO_STATE_RANDOM_WALK.value);

		
	}
	
	
	
	/**
	 * the load button was pressed
	 * load the simulations and first save the file
	 */
	private void loadSimulation(){
		
		if(simulatorThread!=null){
			simulatorThread = null;
		}
		
		if(!loaded){
			disableSetupArea();
			this.simulationsLoad.setEnabled(false);

			simulationsProgress.setString(GUIconstants.SIMULATION_PROGRESS_INIT);

			saveModel("");
					
			double stopTimeValue = simulationParamaters.STOP_TIME.getValue();
			
			if(stopTimeValue/this.steps > GUIconstants.SIMULATION_PROGRESS_MAX_STEP){
				this.steps = (int) Math.ceil((double)stopTimeValue/GUIconstants.SIMULATION_PROGRESS_MAX_STEP);
			} else{
				this.steps = GUIconstants.SIMULATION_PROGRESS_MAX;
			}
			
			this.simulationsProgress.setMaximum(this.steps*this.ip.ENSEMBLE_SIZE.value);
			
			statusTextArea.setText(null);
			
			simulatorThread = new SimulatorThread(currentFile, this,this.steps);
			
			if(simulatorThread.isInitialised()){
			
				loaded = true;
				
				
				this.simulationsLoad.setIcon(this.buttonUnloadImg);
				this.simulationsLoad.setText(GUIconstants.UNLOAD_BUTTON);
				this.simulationsLoad.setEnabled(true);
				this.simulationsStart.setEnabled(true);
				this.timeLabel.setEnabled(true);
			} else{
				enableSetupArea();
				this.simulationsLoad.setEnabled(true);
				  JOptionPane.showMessageDialog(frame, GUIconstants.MESSAGE_NOT_INITIALISED, GUIconstants.MENU_HELP_ABOUT, JOptionPane.INFORMATION_MESSAGE, logoImg);
			}
			simulationsProgress.setString("");
		} else{
			enableSetupArea();
			loaded = false;
			this.simulationsLoad.setIcon(this.buttonLoadImg);
			this.simulationsLoad.setText(GUIconstants.LOAD_BUTTON);
			this.simulationsStart.setEnabled(false);
			this.simulationsProgress.setValue(0);
			this.simulationsProgress.setString("");
			this.simulationsProgress.setEnabled(false);
			this.timeLabel.setEnabled(false);
		}
		simulating = false;
		paused=false;
	}
	
	/**
	 * start/pause simulation was pressed
	 */
	private void startSimulation(){
		this.simulationsStart.setEnabled(false);
		if(!simulating){

			if(simulatorThread==null){
				simulatorThread = new SimulatorThread(currentFile, this,GUIconstants.SIMULATION_PROGRESS_MAX-GUIconstants.SIMULATION_PROGRESS_MIN);
			}
			
			this.simulationsStart.setIcon(this.buttonPauseImg);
			simulationsStart.setText(GUIconstants.SIMULATE_BUTTON_PAUSE);

		
			if(!paused){
				this.simulationsLoad.setEnabled(false);
				simulationsProgress.setEnabled(true);
				this.timeLabel.setEnabled(true);
				simulatorThread.start();
			}else{
				synchronized (simulatorThread){
					simulatorThread.resumeSimulation();
				}
			}
			simulating=true;
			paused = false;
			this.simulationsLoad.setEnabled(false);
			this.simulationsStart.setEnabled(true);
		} else{
			this.simulationsStart.setIcon(this.buttonPlayImg);
			simulationsStart.setText(GUIconstants.SIMULATE_BUTTON_START);
			simulating = false;
			paused=true;
			synchronized (simulatorThread)  {
				simulatorThread.pauseSimulation();
			} 
			this.simulationsLoad.setEnabled(true);
		}
		this.simulationsStart.setEnabled(true);
	}
	
	/**
	 * generates the printed text on the time label
	 * @param elapsedTime
	 * @param estimatedTime
	 * @return
	 */
	private String getTimeString(double elapsedTime, double estimatedTime){
		return GUIconstants.TIME_ELAPSED+" "+formatTime(elapsedTime)+"  / "+GUIconstants.TIME_ESTIMATED+" "+formatTime(estimatedTime)+"";
	}
	
	
	private String formatTime(double time){
		long millis = (long) (time * 1000);
		return df.format(new Date(millis));
	}
	
	 
	 /**
	  * updates simulation progress
	  * @param i step
	  * @param elapsedTime
	  * @param estimatedTime
	  */
	 public void updateProgress(int i, double elapsedTime, double estimatedTime, int ensamble){
		 
		 	double value = Utils.roundTwoDecimals((double)(i*GUIconstants.SIMULATION_PROGRESS_MAX)/(this.steps*ensamble));
		 
			simulationsProgress.setString(value+" %");
		 
		 	
		 	
			simulationsProgress.setValue(i);
			timeLabel.setText(getTimeString(elapsedTime,estimatedTime));
	 }
	 
	 
	 
	 
	 /**
	  * disables the setup area
	  */
	 private void disableSetupArea(){
		 setupAreaTabbedPane.setVisible(false);
		 statusAreaScroll.setVisible(true); 
	 }
	 
	 /**
	  * enables the setup area
	  */
	 private void enableSetupArea(){
		 setupAreaTabbedPane.setVisible(true);
		 statusAreaScroll.setVisible(false);
	 
	 }
	 
	 /**
	  * when simulations are finished, then do...
	  */
	 public void finishedSimulations(){
		 simulating = false;
		 paused=false;
		 simulatorThread = null;
		 this.simulationsLoad.setEnabled(true);
		 this.simulationsStart.setIcon(this.buttonPlayImg);
		 simulationsStart.setText(GUIconstants.SIMULATE_BUTTON_START);
		 simulationsProgress.setString("100 %");
		 simulationsProgress.setValue(GUIconstants.SIMULATION_PROGRESS_MAX);
		 this.simulationsProgress.setEnabled(false);
		 this.timeLabel.setEnabled(false);
	 }
	 
	 /**
	  * enables the simulation start button
	  * @param enabled
	  */
	 public void setEnableSimulationStart(boolean enabled){
		 simulationsStart.setEnabled(enabled);
	 }
	
}
