package utils;

import javax.swing.*;

/**
 * general constants used in the application
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class Constants {

	// application
	public static final int NONE = -1; // constant that encodes for no value
	public static final int FIRST = 0; // constant that encodes for the index of the first element
	public static final char CSV_FILE_TEXT_DELIMITER = '\"'; // char that wrapps text in a cell
	public static final char[] CSV_FILE_CELL_DELIMTER = {',',';'}; // delimiters that separate cells
	
	// FG
	public static final int MAX_STRAND_LEN_TO_PRINT_STRAND_STATE_ARRAYS = 50;
	
	// default parameters file
	public static final String DEFAULT_PARAMS_FILE = "system.ini";
	public static final String DEFAULT_PARAMS_FILE_WIN = "system_win.ini";
	public static final String DEFAULT_PARAMS_FILE_EMPTY = "system_empty.ini";
	public static final String PARAMS_FILE_COMMENT_CHAR = "#";
	public static final String PARAMS_FILE_LINE_ENDING = ";";
	public static final String PARAMS_FILE_ASSIGNMENT_CHAR = "=";
	public static final String PARAMS_FILE_PREFIX = "params_";
	
	public static final String[] PARSER_DINUCLEOTIDE_CSV_FILE_HEADER = {"DINUCLEOTIDE", "PROPELLLERTWISTFACTOR", "SLIDEFACTOR"};

	public static final String[] PARSER_4MER_CSV_FILE_HEADER = {"POLYASEQ", "FACTOR"};

	public static final String[] PARSER_TF_CSV_FILE_HEADER = {"NAME", "DBD", "ES", "COPYNUMBER", "SIZELEFT", "SIZERIGHT", "ASSOCRATE", "INITIALDROP", "UNBINDINGPROBABILITY", "SLIDELEFTPROBABILITY", "SLIDERIGHTPROBABILITY", "JUMPINGPROBABILITY", "HOPSTDDISPLACEMENT", "SPECIFICWAITINGTIME", "STEPLEFTSIZE", "STEPRIGHTSIZE", "UNCORRELATEDDISPLACEMENTSIZE", "STALLSIFBLOCKED", "COLLISIONUNBINDPROBABILITY", "AFFINITYLANDSCAPEROUGHNESS", "PREBOUNDPROPORTION", "PREBOUNDTOHIGHESTAFFINITY", "TFISIMMOBILE", "ISBIASEDRANDOMWALK", "ISTWOSTATERANDOMWALK", "REPRESSIONRATE", "DEREPRESSIONATTENUATIONFACTOR", "REPRLENLEFT", "REPRLENRIGHT", "SPECIFICENERGYTHRESHOLD"};
	public static final String TF_CSV_FILE_TARGET_SITES_DELIMITER = ";";
	public static final String AFFINITY_CSV_FILE_DELIMITER = ",";



	
	public static final String[] PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER = {"FIRSTSPECIES", "SECONDSPECIES", "TYPE", "FIRSTSPECIESDIRECTION", "SECONDSPECIESDIRECTION", "FIRSTSPECIESDNAREGION", "SECONDSPECIESDNAREGION", "DIMMERISATIONPROBABILITY", "AFFINITYINCREASE", "ISREVERSIBLE", "ISADDITIVE", "ISFIXED"};
	
	public static final String FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER = ":";
	public static final String FASTA_FILE_DESCRIPTION_INTERVAL_DELIMITER_REGEX = "\\.\\.";
	public static final String FASTA_FILE_DESCRIPTION_INTERVAL_DELIMITER = "..";

	public static final int TF_SPECIES_FIRST_ID = 1;

	public static final int EVENT_TF_BINDING = 0;
	public static final int EVENT_TF_RANDOM_WALK_SLIDE_LEFT = 1;
	public static final int EVENT_TF_RANDOM_WALK_SLIDE_RIGHT = 2;
	public static final int EVENT_TF_RANDOM_WALK_JUMP = 3;
	public static final int EVENT_TF_RANDOM_WALK_HOP = 4;
	public static final int EVENT_TF_OLIGOMER = 5;
	public static final int EVENT_TF_REPRESSION = 6;
	public static final int EVENT_TF_DEREPRESSION = 7;

	public static final byte BP_IS_OPEN = 0;
	public static final byte BP_IS_CLOSED = 1;
	public static final byte BP_IS_REPRESSED = 2;

	public static final double MAX_PROBABILITY=1.0;
	
	public static final int MAX_STRING_LENGTH=1000;
	

	public static final int TF_COOPERATIVITY_TYPE_DIRECT = 0;
	public static final int TF_COOPERATIVITY_TYPE_DNA = 1;

	
	public static final String PARAMETR_FILE_EXTENSION = ".grp";
	public static final String STATUS_FILE_EXTENSION = ".txt";
	public static final String OCCUPANCY_FILE_EXTENSION = ".wig";
	public static final String AFFINITY_FILE_EXTENSION = ".wig";
	public static final String TF_FILE_EXTENSION = ".csv";
	public static final String TARGET_SITE_FILE_EXTENSION = ".csv";
	public static final String SLIDING_LENGTHS_FILE_EXTENSION = ".csv";
	public static final String SEQUENCE_FILE_EXTENSION = ".fasta";
	public static final String BACKUP_FILE_EXTENSION = ".bkp";
	public static final String INFO_FILE_EXTENSION = ".txt";
	
	public static final String PFM_NUCLEOTIDE_SEPARATOR = ";";
	public static final String PFM_NUCLEOTIDE_ASSIGNMENT = "=";
	public static final String PFM_NUCLEOTIDE_CONTAINER_LEFT = "(";
	public static final String PFM_NUCLEOTIDE_CONTAINER_RIGHT = ")";
	public static final String PFM_NUCLEOTIDE_CONTAINER_SEPRATOR = ",";

	public static final String DBD_TYPE_SEQ = "SEQ:";
	public static final String DBD_TYPE_PFM = "PFM:";
	public static final String DBD_TYPE_PWM = "PWM:";
	public static final String DBD_TYPE_PFM_INFO = "INFO:";
	public static final String DBD_TYPE_PFM_ENERGY = "ENERGY:";
	public static final String DBD_TYPE_SEPARATOR = ":";
	public static final double DBD_TYPE_PFM_ENERGY_CORRECTION_DEFAULT = 1;
	public static final String DBD_TYPE_LANDSCAPE = "LANDSCAPE:";
	public static final String DBD_TYPE_SEQS = "SEQS:";

	public static final boolean USE_INFORMATION_THEORY_TO_COMPUTE_BINDING_ENERGIES = true;
	public static final int MAXIMUM_NUMBER_OF_TF_MOLECULES_TO_FOLLOW=10;
	
	//event types
	public static final int NEXT_EVENT_IS_NONE = 0;
	public static final int NEXT_EVENT_IS_TF_BINDING = 1;
	public static final int NEXT_EVENT_IS_TF_RANDOM_WALK = 2;
	public static final int NEXT_EVENT_IS_TF_REPRESSION = 5;

	public static final int NEXT_EVENT_IS_NON_COGNATE_TF_BINDING = 3;
	public static final int NEXT_EVENT_IS_NON_COGNATE_TF_RANDOM_WALK = 4;	
	
	
	public static final String DNA_FASTA_SUBSEQUENCE = "subsequence";
	public static final String DNA_FASTA_COPY_NUMBER = "copy";
	public static final String DNA_FASTA_BOUNDARY = "boundary";
	public static final String DNA_FASTA_BOUNDARY_REFLEXIVE = "reflexive";
	public static final String DNA_FASTA_BOUNDARY_PERIODIC = "periodic";
	public static final String DNA_FASTA_BOUNDARY_ABSORBING = "absorbing";
	public static final String DNA_FASTA_DELIMITER = ";";
	
	public static final double DOUBLE_ZERO=1E-14;
	
	// save states after at least 4 hours 14400s
	public static final int MIN_INTERMEDIARY_STATE=1;

	public static final int NO_OF_TOP_SITES = 3;
	
	public static final boolean IS_BIAS_RANDOM_WALK=false;
	
}
