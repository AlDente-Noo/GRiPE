package environment;

import agents.DBP;
import agents.Remodeller;
import agents.TF;
import event.Event;
import event.EventList;
import event.ProteinEvent;
import event.RepressionEvent;
import objects.*;
import simulator.SimulatorGUI;
import utils.*;

import java.io.*;
import java.util.*;


/**
 * class that simulates the behaviour of a nucleus
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
@SuppressWarnings("ALL")
public class Cell implements Serializable {

    private static final long serialVersionUID = 5826964812877420521L;

    //parameters set externally
    public InputParameters ip;

    //internal parameters
    public Random randomGenerator;
    public DNA dna;
    public TFSpecies[] TFspecies;
    public int maxRepressionLeftSize, maxRepressionRightSize, maxTFSize;
    public TargetSitesAndGroups tsg;
    public DBP[] dbp; //DNA binding proteins
    public Remodeller remodeller;
    public int moleculesCopyNumber;
    public int ensemble;

    public ArrayList<ArrayList<Integer>> TFmoleculesIDs;
    public ArrayList<ArrayList<Integer>> freeTFmolecules;
    public int freeTFmoleculesTotal;
    public int TFreadingDirection;
    public double cellTime;
    public double totalStopTime;
    public double totalSimulatedTime;
    public double totalElapsedTime;
    public boolean isPartialSimulation;
    public EventList eventQueue;
    public long eventsCount;
    public int targetSitesReached;
    public boolean areTargetSitesToBeReached;

    //output files
    public String outputStatusFile;
    public File outputParamsFile;
    public String outputPath;
    public String outputDNAOccupancyFile;
    public String outputDNASequenceFile;
    public String outputAffinityLandscapeFile;
    public String outputTFFile;
    public String outputTargetSiteFollowFile;
    public String outputDNAOccupancyFinalFile;
    public String outputTargetSiteFile;
    public String outputIntermediaryBackupFile;
    public String outputIntermediaryInfoFile;
    public ArrayList<String> TargetSiteFollowLines;
    public SimulatorGUI gui;
    public boolean runUntilTSReached;
    private double lastPrintResultsAfter;
    private HashMap<Integer, String> TFstoFollow;
    private HashMap<Integer, String> TFstoFollowFiles;
    private ArrayList<Integer> TFstoFollowID;
    private boolean areTFstoFollow;
    private double nextPointTFstoFollow;
    private double stepPointTFstoFollow;

    //private boolean simulationFinished;
    //public BufferedWriter statusBuffer;

    /**
     * constructor that initialises a nucleus
     *
     * @throws FileNotFoundException
     */
    public Cell(String parametersFile, SimulatorGUI gui, boolean newParamsFile) throws FileNotFoundException {

        this.gui = gui;

        // FG: parse input parameters (params.grp)
        ip = new InputParameters(parametersFile);

        String paramsFile = parametersFile;
        if (newParamsFile) {
            paramsFile = "";
        }

        this.outputParamsFile = ip.exportParameterFile(paramsFile);
        this.outputPath = ip.outputDirectory;

        generateOutputFilenames(this.outputParamsFile);

        //initialise internal parameters
        // FG: here TF files are parsed
        initialiseInternalParameters();

        printInitInfo();
        printPreprocesedInfo();
    }


    /**
     * generate output filename
     *
     * @param paramsFilename the filename of the params
     */
    private void generateOutputFilenames(File paramsFilename) {
        this.outputStatusFile =
                paramsFilename.getName().replaceAll("params", "status").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.STATUS_FILE_EXTENSION);
        this.outputDNAOccupancyFile =
                paramsFilename.getName().replaceAll("params", "occupancy").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.OCCUPANCY_FILE_EXTENSION);
        this.outputTargetSiteFile =
                paramsFilename.getName().replaceAll("params", "target_site").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.TARGET_SITE_FILE_EXTENSION);
        this.outputTargetSiteFollowFile =
                paramsFilename.getName().replaceAll("params", "target_site_follow").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.TF_FILE_EXTENSION);
        this.outputDNASequenceFile =
                paramsFilename.getName().replaceAll("params", "DNA_seq").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.SEQUENCE_FILE_EXTENSION);
        this.outputAffinityLandscapeFile =
                paramsFilename.getName().replaceAll("params", "affinity_landscape").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.AFFINITY_FILE_EXTENSION);
        this.outputTFFile =
                paramsFilename.getName().replaceAll("params", "TF_species").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.TF_FILE_EXTENSION);
        this.outputDNAOccupancyFinalFile =
                paramsFilename.getName().replaceAll("params", "DNA_occupancy_final").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.OCCUPANCY_FILE_EXTENSION);
        this.outputIntermediaryBackupFile =
                "backup_" + paramsFilename.getName().replaceAll("_params", "").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.BACKUP_FILE_EXTENSION);
        this.outputIntermediaryInfoFile =
                "backup_" + paramsFilename.getName().replaceAll("_params", "").replaceAll(Constants.PARAMETR_FILE_EXTENSION, Constants.INFO_FILE_EXTENSION);
    }

    /**
     * initialises the internal parameters
     */
    private void initialiseInternalParameters() {

        //create the random number generator
        createRandomNumberGenerator();

        this.cellTime = 0;
        this.totalSimulatedTime = 0;
        this.totalStopTime = this.ip.STOP_TIME.value;
        this.totalElapsedTime = 0;
        this.isPartialSimulation = false;
        this.lastPrintResultsAfter = 0;
        this.eventsCount = 0;
        this.TFreadingDirection = 1;
        if (this.ip.TF_READ_IN_BOTH_DIRECTIONS.value) {
            this.TFreadingDirection = 2;
        }

        // FG: here files are parsed
        generateObjects();

        //TS to be reached;
        areTargetSitesToBeReached = dna.areTargetSites;
        this.targetSitesReached = 0;
        initTargetSitesToFollow();

        TFstoFollow = new HashMap<Integer, String>();
        TFstoFollowFiles = new HashMap<Integer, String>();
        TFstoFollowID = new ArrayList<Integer>();
        nextPointTFstoFollow = this.totalStopTime;
        stepPointTFstoFollow = this.totalStopTime;

        getTFsToFollow();

        eventQueue = new EventList(this);

        this.bindMolecules();
        this.ensemble = 0;

        this.runUntilTSReached = false;
    }


    /**
     * initialises the internal parameters
     */
    private void resetInternalParameters() {

        //create the random number generator
        createRandomNumberGenerator();

        this.cellTime = 0;
        eventQueue = new EventList(this);
        this.unbindMoleculesAndDerepressDNA();
        this.bindMolecules();

        for (DBP value : this.dbp) {
            assert !value.isRepressingDNA();
        }

        //reset target site reached
        for (int i = 0; i < TFspecies.length; i++) {
            for (int j = 0; j < dna.strand.length; j++) {
                for (int k = 0; k < this.TFreadingDirection; k++) {
                    this.dna.firstReached[i][j][k] = Constants.NONE;
                }
            }
        }

        //TS to be reached;
        areTargetSitesToBeReached = dna.areTargetSites;
        this.targetSitesReached = 0;
        this.runUntilTSReached = false;
    }


    /**
     * initiates the random number generator
     */
    private void createRandomNumberGenerator() {
        randomGenerator = new Random();
        if (this.ip.RANDOM_SEED.value > 0) {
            randomGenerator = new Random(this.ip.RANDOM_SEED.value);
        }
    }


    /**
     * generates the TFs to follows based on a string
     */
    private void getTFsToFollow() {
        areTFstoFollow = false;
        double avgStep = 1.0, bufferStep;

        if (this.ip.ENSEMBLE_SIZE.value == 1 && !this.ip.OUTPUT_TF.value.isEmpty() && this.ip.OUTPUT_TF_POINTS.value > 0) {
            String buffer, filename;
            StringBuilder strBuf = new StringBuilder();
            buffer = this.ip.OUTPUT_TF.value + ",";

            StringTokenizer st = new StringTokenizer(buffer, ",");
            int id;
            while (st.hasMoreTokens()) {
                buffer = st.nextToken().trim();
                id = this.getTFspeciesID(buffer);
                if (id != Constants.NONE) {
                    if (TFmoleculesIDs.get(id).size() > 0) {
                        filename = this.outputTFFile.replaceAll("TF_species", "TF_" + buffer);
                        strBuf.delete(0, strBuf.length());
                        strBuf.append("time, freeMolecules, freeMoleculesProportion, freePositions, " +
                                "freePositionsProportion");

                        if (TFmoleculesIDs.get(id).size() <= Constants.MAXIMUM_NUMBER_OF_TF_MOLECULES_TO_FOLLOW) {
                            for (int i = 0; i < TFmoleculesIDs.get(id).size(); i++) {
                                strBuf.append(", TF");
                                strBuf.append(TFmoleculesIDs.get(id).get(i));
                            }

                            bufferStep = CellUtils.computeAvgMoveRate(TFspecies[id].specificWaitingTime, 0);
                            if (bufferStep < avgStep) {
                                avgStep = bufferStep;
                            }
                        }
                        BufferedWriter bufferFile = null;
                        try {
                            //Construct the BufferedWriter object
                            if (this.outputPath.isEmpty()) {
                                bufferFile = new BufferedWriter(new FileWriter(filename));
                            } else {
                                bufferFile = new BufferedWriter(new FileWriter(new File(this.outputPath, filename)));
                            }
                            bufferFile.write(strBuf.toString());
                            bufferFile.newLine();
                            this.TFstoFollowID.add(id);
                            this.TFstoFollow.put(id, buffer);
                            this.TFstoFollowFiles.put(id, filename);
                            if (!this.areTFstoFollow) {
                                this.areTFstoFollow = true;
                            }
                            bufferFile.flush();
                            bufferFile.close();
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }
                    } else {
                        this.printDebugInfo("There is no molecule of TF " + buffer);
                    }
                } else {
                    this.printDebugInfo("TF " + buffer + " to follow was not find");
                }
            }
        }

        if (areTFstoFollow) {
            if (this.totalStopTime > 0) {
                this.nextPointTFstoFollow = this.totalStopTime / (double) this.ip.OUTPUT_TF_POINTS.value;
                stepPointTFstoFollow = this.nextPointTFstoFollow;
            } else {
                this.nextPointTFstoFollow = avgStep;
                stepPointTFstoFollow = avgStep;
            }
        }
    }


    /**
     * generates the Target Sites to follow file based on a string
     */
    private void initTargetSitesToFollow() {
        if (this.ip.ENSEMBLE_SIZE.value == 1 && this.ip.FOLLOW_TS.value) {
            BufferedWriter bufferFile = null;
            try {
                //Construct the BufferedWriter object
                if (this.outputPath.isEmpty()) {
                    bufferFile = new BufferedWriter(new FileWriter(this.outputTargetSiteFollowFile));
                } else {
                    bufferFile = new BufferedWriter(new FileWriter(new File(this.outputPath,
                            outputTargetSiteFollowFile)));
                }
                bufferFile.write("time, " + this.tsg.getTargetSiteGroupsString());
                bufferFile.newLine();

                bufferFile.write("0, " + this.tsg.getTargetSiteGroupsOccupancyString());
                bufferFile.newLine();

                bufferFile.flush();
                bufferFile.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            TargetSiteFollowLines = new ArrayList<String>();
        } else if (this.ip.FOLLOW_TS.value) {
            this.printDebugInfo("Ensemble size is bigger than 1, TS cannot be followed.");
        }
    }


    /**
     * either load the sequence from a file or generate it randomly
     */
    private void createDNAstrand() {
        DNA bufferDNA = DNAFilesParser.fastaFileParser(this.ip.DNA_SEQUENCE_FILE.value);
        if (bufferDNA != null && bufferDNA.strand != null && bufferDNA.strand.length > 0) {
            dna = new DNA(bufferDNA);
        } else {
            //generate the DNA strand randomly
            dna = new DNA(randomGenerator, this.ip.DNA_LENGTH.value, this.ip.DNA_PROPORTION_OF_A.value,
                    this.ip.DNA_PROPORTION_OF_T.value, this.ip.DNA_PROPORTION_OF_C.value,
                    this.ip.DNA_PROPORTION_OF_G.value, this.ip.DNA_BOUNDARY_CONDITION.value);
        }
        String str = DNAFilesParser.btrackFileParser(this.ip.DNA_AVAILABILITY_FILE.value, dna.closed);
        if (!str.equals("")) {
            printDebugInfo("Warning: " + str + "; the dnase accessibility (btrack) file is ignored");
        }
    }


    /**
     * either load them from a file or create them randomly
     */
    private void createTFSpecies() {
        moleculesCopyNumber = 0;
        maxRepressionRightSize = 0;
        maxRepressionLeftSize = 0;
        maxTFSize = 0;

        // load TFs from file
        TFfileParser TFparser = new TFfileParser(this, this.ip.TF_FILE.value,
                Utils.generateNextInteger(randomGenerator, this.ip.TF_COPY_NUMBER_MIN.value,
                        this.ip.TF_COPY_NUMBER_MAX.value), this.ip.TF_ES.value, this.ip.TF_SIZE_LEFT.value,
                this.ip.TF_SIZE_RIGHT.value, this.ip.TF_ASSOC_RATE.value, this.dna.strand.length,
                this.ip.TF_UNBINDING_PROBABILITY.value, this.ip.TF_SLIDE_LEFT_PROBABILITY.value,
                this.ip.TF_SLIDE_RIGHT_PROBABILITY.value, this.ip.TF_JUMPING_PROBABILITY.value,
                this.ip.TF_HOP_STD_DISPLACEMENT.value, this.ip.TF_SPECIFIC_WAITING_TIME.value,
                this.ip.TF_STEP_LEFT_SIZE.value, this.ip.TF_STEP_RIGHT_SIZE.value,
                this.ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.value,
                this.ip.TF_STALLS_IF_BLOCKED.value, this.ip.TF_COLLISION_UNBIND_PROBABILITY.value,
                this.ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value, this.ip.TF_PREBOUND_PROPORTION.value,
                this.ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.value, this.ip.TF_IS_IMMOBILE.value, dna.subsequence,
                this.ip.IS_BIASED_RANDOM_WALK.value, this.ip.IS_TWO_STATE_RANDOM_WALK.value, this.ip.TF_SPECIFIC_ENERGY_THRESHOLD.value,
                this.ip.TF_REPRESSION_RATE.value, this.ip.TF_DEREPRESSION_ATTENUATION_FACTOR.value, this.ip.TF_REPR_LEN_LEFT.value,
                this.ip.TF_REPR_LEN_RIGHT.value, this.ip.TF_TAU.value);
        if (TFparser.parsed) {
            // load TF species
            if (TFparser.data != null && !TFparser.data.isEmpty()) {
                TFspecies = new TFSpecies[TFparser.data.size()];
                for (int i = 0; i < TFparser.data.size(); i++) {
                    TFspecies[i] = TFparser.data.get(i);
                    moleculesCopyNumber += TFspecies[i].copyNumber;
                    if (maxRepressionLeftSize < TFspecies[i].repressionLeftSize) {
                        maxRepressionLeftSize = TFspecies[i].repressionLeftSize;
                    }
                    if (maxRepressionRightSize < TFspecies[i].repressionRightSize) {
                        maxRepressionRightSize = TFspecies[i].repressionRightSize;
                    }
                    if (maxTFSize < TFspecies[i].sizeTotal) {
                        maxTFSize = TFspecies[i].sizeTotal;
                    }
                }
            }

            //cooperativity data
            if (!this.ip.TF_COOPERATIVITY_FILE.value.isEmpty()) {
                TFcooperativityFileParser TFcooperativityParser = new TFcooperativityFileParser(this,
                        this.ip.TF_COOPERATIVITY_FILE.value, Constants.NONE, Constants.NONE, new DNAregion("", 0,
                        this.dna.strand.length), Constants.MAX_PROBABILITY, dna.region);
                if (TFcooperativityParser.parsed) {
                    // load TF species
                    if (TFcooperativityParser.data != null && !TFcooperativityParser.data.isEmpty()) {

                        TFcooperativity buffer;
                        this.printDebugInfo("Cooperativity data is loaded from the file " + this.ip.TF_COOPERATIVITY_FILE.value);

                        for (int i = 0; i < TFcooperativityParser.data.size(); i++) {
                            //rescale region 1 and region 0 to remove the size of the TF
                            TFcooperativityParser.data.get(i).region0.end =
                                    Math.max(TFcooperativityParser.data.get(i).region0.start + 1,
                                            TFcooperativityParser.data.get(i).region0.end - TFspecies[TFcooperativityParser.data.get(i).species0ID].sizeTotal);
                            TFcooperativityParser.data.get(i).region1.end =
                                    Math.max(TFcooperativityParser.data.get(i).region1.start + 1,
                                            TFcooperativityParser.data.get(i).region1.end - TFspecies[TFcooperativityParser.data.get(i).species1ID].sizeTotal);

                            TFspecies[TFcooperativityParser.data.get(i).species0ID].addCooperativity(TFcooperativityParser.data.get(i), this.dna.strand.length, this.TFreadingDirection);
                            if (TFcooperativityParser.data.get(i).type == Constants.TF_COOPERATIVITY_TYPE_DIRECT && TFcooperativityParser.data.get(i).species0ID != TFcooperativityParser.data.get(i).species1ID) {
                                buffer = new TFcooperativity(TFcooperativityParser.data.get(i));
                                buffer.swap();
                                TFspecies[TFcooperativityParser.data.get(i).species1ID].addCooperativity(buffer,
                                        this.dna.strand.length, this.TFreadingDirection);
                            }
                        }
                    } else {
                        this.printDebugInfo("Could not find any valid cooperativity data in file "
                                + "'" + this.ip.TF_COOPERATIVITY_FILE.value + "'");
                    }
                }
            }

            //target site file
            if (!this.ip.TS_FILE.value.isEmpty()) {
                TSfileParser TSfileParser = new TSfileParser(this.ip.TS_FILE.value, this);
                tsg = TSfileParser.tsg;
            }
        } else {
            this.printDebugInfo("The TF cooperativity file '" + this.ip.TF_COOPERATIVITY_FILE.value + "' is not parsed");
        }

        // if could not load them then generate TF species randomly
        if (TFspecies == null || TFspecies.length <= 0) {
            int size = this.ip.TF_SPECIES_COUNT.value;
            if (TFspecies != null) {
                size += TFspecies.length;
            }
            if (this.ip.TF_COPY_NUMBER_MAX.value > 0) {
                tsg = new TargetSitesAndGroups();
                TFspecies = new TFSpecies[size];
                int dbdLength, copyNumber, step, pos;
                step = dna.strand.length / (TFspecies.length + 1);
                int i = 0;
                while (i < TFspecies.length) {
                    dbdLength = Utils.generateNextInteger(randomGenerator, this.ip.TF_DBD_LENGTH_MIN.value,
                            this.ip.TF_DBD_LENGTH_MAX.value);
                    copyNumber = Utils.generateNextInteger(randomGenerator, this.ip.TF_COPY_NUMBER_MIN.value,
                            this.ip.TF_COPY_NUMBER_MAX.value);
                    if (copyNumber > 0) {
                        pos = step * (i + 1);
                        TFspecies[i] = new TFSpecies(i, dna.strand, pos, dbdLength, copyNumber, this.ip.TF_ES.value,
                                this.ip.TF_SIZE_LEFT.value, this.ip.TF_SIZE_RIGHT.value,
                                this.ip.TF_ASSOC_RATE.value, new DNAregion("", 0, this.dna.strand.length), true,
                                this.ip.TF_UNBINDING_PROBABILITY.value, this.ip.TF_SLIDE_LEFT_PROBABILITY.value,
                                this.ip.TF_SLIDE_RIGHT_PROBABILITY.value, this.ip.TF_JUMPING_PROBABILITY.value,
                                this.ip.TF_HOP_STD_DISPLACEMENT.value, this.ip.TF_SPECIFIC_WAITING_TIME.value,
                                this.ip.TF_STEP_LEFT_SIZE.value, this.ip.TF_STEP_RIGHT_SIZE.value,
                                this.ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.value,
                                this.ip.TF_STALLS_IF_BLOCKED.value, this.ip.TF_COLLISION_UNBIND_PROBABILITY.value,
                                this.ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value, this.ip.TF_PREBOUND_PROPORTION.value,
                                this.ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.value, this.ip.TF_IS_IMMOBILE.value,
                                this.ip.IS_BIASED_RANDOM_WALK.value,
                                this.ip.IS_TWO_STATE_RANDOM_WALK.value, this.ip.TF_SPECIFIC_ENERGY_THRESHOLD.value,
                                this.ip.TF_REPRESSION_RATE.value, this.ip.TF_DEREPRESSION_ATTENUATION_FACTOR.value,
                                this.ip.TF_REPR_LEN_LEFT.value, this.ip.TF_REPR_LEN_RIGHT.value, this.ip.TF_TAU.value);
                        moleculesCopyNumber += TFspecies[i].copyNumber;
                        tsg.addGroup(this, pos, i);
                        i++;
                    }
                }
            } else {
                //this.printDebugInfo("TF species cannot have even at least 1 molecule. None were created!");
                this.printDebugInfo("Molecules of this TF species have not been created");
            }
        }
    }

    /**
     * create all molecules TFs
     */
    private void createMolecules() {
        dbp = new DBP[this.moleculesCopyNumber];
        int id = 0;

        TFmoleculesIDs = new ArrayList<ArrayList<Integer>>();
        freeTFmolecules = new ArrayList<ArrayList<Integer>>();
        freeTFmoleculesTotal = 0;
        for (int i = 0; i < TFspecies.length; i++) {
            TFmoleculesIDs.add(new ArrayList<Integer>());
            freeTFmolecules.add(new ArrayList<Integer>());
            for (int j = 0; j < TFspecies[i].copyNumber; j++) {
                dbp[id] = new TF(id, Constants.NONE, Constants.NONE, 0, i, TFspecies[i].sizeTotal, Constants.NONE,
                        TFspecies[i].hasDNAbasedCooperativity, TFspecies[i].hasDirectCooperativity);
                TFmoleculesIDs.get(i).add(id);
                freeTFmolecules.get(i).add(id);
                freeTFmoleculesTotal++;
                id++;
            }
        }
    }

    private void createRemodeller() {
        this.remodeller = new Remodeller(ip.DNA_DEREPRESSION_RATE.value);
    }

    /**
     * pre-bind a bulk of molecules on the DNA
     */
    public void bindMolecules() {
        this.createRandomNumberGenerator();

        int moleculesToBind, bound, newPosition;
        int[] newLocation;
        //prebind cognate TF
        for (int i = 0; i < TFspecies.length; i++) {
            if (TFspecies[i].preboundProportion > 0) {
                moleculesToBind = (int) Math.round(TFspecies[i].copyNumber * TFspecies[i].preboundProportion);
                if (moleculesToBind > 0) {
                    if (TFspecies[i].preboundToHighestAffinity) {
                        bound = 0;
                        for (int j = 0; j < moleculesToBind && bound != Constants.NONE; j++) {
                            newLocation = this.getStrongestAvailableSite(i);
                            newPosition = newLocation[0];
                            if (newPosition != Constants.NONE) {
                                bound = this.bindTFMoleculeToPosition(i, this.cellTime, newPosition, newLocation[1]);
                                if (bound != Constants.NONE && !this.TFspecies[this.dbp[bound].speciesID].isImmobile) {
                                    eventQueue.scheduleNextTFRandomWalkEvent(this, bound, this.cellTime);
                                }
                            }
                        }
                    } else {
                        bound = 0;
                        for (int j = 0; j < moleculesToBind && bound != Constants.NONE; j++) {
                            bound = this.bindTFMolecule(i, this.cellTime);
                            if (bound != Constants.NONE && !this.TFspecies[this.dbp[bound].speciesID].isImmobile) {
                                eventQueue.scheduleNextTFRandomWalkEvent(this, bound, this.cellTime);
                            }
                        }
                    }
                }
            }
        }

        if (this.ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value) {
            this.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(this);
        }
    }

    /**
     * unbind all molecules on the DNA (FG: and derepresses DNA)
     */
    public void unbindMoleculesAndDerepressDNA() {
        // FG: unbind all molecules
        for (DBP value : this.dbp) {
            if (value.getPosition() != Constants.NONE) {
                value.unbindMolecule(this, this.cellTime);
            }
        }
        // FG: execute all derepression events in the queue
        RepressionEvent re = (RepressionEvent) this.eventQueue.TFRepressionEventQueue.pop();
        while (re != null && !re.isEmpty()) {
            if (re.nextAction == Constants.EVENT_TF_DEREPRESSION) {
                re.time = this.cellTime;
                this.remodeller.act(this, re);
            }
            re = (RepressionEvent) this.eventQueue.TFRepressionEventQueue.pop();
        }
        this.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(this);
    }


    /**
     * attempts to bind a TF to the DNA
     *
     * @param speciesID the ID of the species of the molecule that binds
     * @param time      the time
     */
    public int bindTFMolecule(int speciesID, double time) {
        int bound = Constants.NONE;
        int moleculeID = 0;
        int newPosition = Constants.NONE;

        if (this.freeTFmolecules.get(speciesID).size() > 0 && dna.effectiveTFavailabilitySum[speciesID] > 0) {
            moleculeID = this.freeTFmolecules.get(speciesID).get(this.freeTFmolecules.get(speciesID).size() - 1);
            while (bound == Constants.NONE && dna.effectiveTFavailabilitySum[speciesID] > 0) {
                newPosition = getTFbindingPosition(moleculeID, speciesID);
                bound = this.dbp[moleculeID].bindMolecule(this, time, newPosition);
            }
        }
        return bound;
    }


    /**
     * attempts to bind a TF to a specific position on the DNA
     *
     * @param speciesID the ID of the species of the molecule that binds
     * @param time      the time
     */
    public int bindTFMoleculeToPosition(int speciesID, double time, int newPosition, int direction) {
        int bound = Constants.NONE;

        int moleculeID = 0;

        if (this.freeTFmolecules.get(speciesID).size() > 0 && dna.effectiveTFavailabilitySum[speciesID] > 0) {
            moleculeID = this.freeTFmolecules.get(speciesID).get(this.freeTFmolecules.get(speciesID).size() - 1);
            while (bound == Constants.NONE && dna.effectiveTFavailabilitySum[speciesID] > 0) {
                bound = this.dbp[moleculeID].bindMolecule(this, time, newPosition, direction);
            }
        }
        return bound;

    }

    /**
     * computes a binding position for a TF on the DNA
     *
     * @param moleculeID
     * @param speciesID
     * @return
     */
    public int getTFbindingPosition(int moleculeID, int speciesID) {
        int newPosition;
        double random = randomGenerator.nextDouble();

        if (this.dbp[moleculeID].wasBound() || TFspecies[speciesID].relInitialDrop == null
                || (TFspecies[speciesID].relInitialDrop.start <= 0
                && TFspecies[speciesID].relInitialDrop.end >= this.dna.strand.length - 1)
                || (this.randomGenerator.nextDouble() > TFspecies[speciesID].initialDrop.probability))
        {
            newPosition = Gillespie.getNextPositionToBind(random * dna.effectiveTFavailabilitySum[speciesID],
                    dna.effectiveTFavailability[speciesID], dna.effectiveTFsectorsAvailabilitySum[speciesID],
                    dna.DNAsectorSize);
        } else {
            int sum = 0;
            boolean[] buffer = new boolean[TFspecies[speciesID].relInitialDrop.size()];
            for (int i = 0; i < TFspecies[speciesID].relInitialDrop.size(); i++) {
                buffer[i] = dna.effectiveTFavailability[speciesID][i + (int) TFspecies[speciesID].relInitialDrop.start];
                if (buffer[i]) {
                    sum++;
                }
            }
            if (sum > 0) {
                newPosition = Gillespie.getNextReaction(random * sum, buffer);
                newPosition += (int) TFspecies[speciesID].relInitialDrop.start;
            } else {
                newPosition = Gillespie.getNextPositionToBind(random * dna.effectiveTFavailabilitySum[speciesID],
                        dna.effectiveTFavailability[speciesID], dna.effectiveTFsectorsAvailabilitySum[speciesID],
                        dna.DNAsectorSize);
            }
        }

        return newPosition;
    }


    /**
     * returns the available site with the highest score.
     *
     * @param speciesID the TF species
     * @return the position or Constants.NONE if none
     */
    public int[] getStrongestAvailableSite(int speciesID) {
        int direction = Constants.NONE;
        int newPosition = Constants.NONE;
        int maxDirection;
        for (int i = 0; i < dna.strand.length; i++) {
            if (dna.effectiveTFavailability[speciesID][i]) {
                maxDirection = 1;
                if (dna.TFavgMoveRate[speciesID][i].length == 1 || dna.TFavgMoveRate[speciesID][i][0] >= dna.TFavgMoveRate[speciesID][i][1]) {
                    maxDirection = 0;
                }
                if (newPosition == Constants.NONE || (dna.TFavgMoveRate[speciesID][newPosition][direction] > dna.TFavgMoveRate[speciesID][i][maxDirection])) {
                    direction = maxDirection;
                    newPosition = i;
                }

            }
        }

        int[] result = new int[2];
        result[0] = newPosition;
        result[1] = direction;

        return result;
    }


    /**
     * get the ID of the species mentioned through the name
     *
     * @param name
     * @return
     */
    public int getTFspeciesID(String name) {
        int result = Constants.NONE;
        name = name.trim();
        for (int i = 0; i < TFspecies.length && result == Constants.NONE; i++) {
            if (name.equals(TFspecies[i].name)) {
                result = TFspecies[i].id;
            }
        }
        return result;
    }

    /**
     * generate TFs and DNA strand random
     */
    private void generateObjects() {

        //DNA
        createDNAstrand();

        //TF species
        createTFSpecies();

        createRemodeller();

        //compute TF affinity landscape
        dna.computeTFaffinityLandscape(this, this.randomGenerator, TFspecies, //this.ip.TF_SPECIFIC_WAITING_TIME.value,
                this.TFreadingDirection, this.ip.DNA_SECTOR_SIZE.value, this.ip.PRINT_FINAL_OCCUPANCY.value);

        createMolecules();
        //System.out.println("created the molecules");
    }


    /**
     * print init information to status file
     */
    private void printInitInfo() {
        printDebugInfo("The model parameters were saved to the file: " + this.outputParamsFile, true);

        if (this.ip.RANDOM_SEED.value > 0) {
            this.printDebugInfo("The specific seed is provided for the random numbers");
        } else {
            this.printDebugInfo("The specific seed is not provided for the random numbers");
        }
        printDebugInfo("Cell simulation time: " + this.ip.STOP_TIME.value + " seconds");
        printDebugInfo("=====================================================");
    }


    /**
     * prints the information for all TFs that are followed
     *
     * @throws FileNotFoundException
     */
    private void printTFtoFollowInformation(double time) throws FileNotFoundException {
        for (int id : this.TFstoFollowID) {
            printTFtoFollowInformation(id, time);
        }
    }


    /**
     * prints the information of an TF
     *
     * @param TFspeciesID the id of the TFSpecies
     * @throws FileNotFoundException
     */
    private void printTFtoFollowInformation(int TFspeciesID, double time) throws FileNotFoundException {
        StringBuffer strBuf = new StringBuffer();
        int pos;
        strBuf.append(time);
        strBuf.append(", ");
        strBuf.append(this.freeTFmolecules.get(TFspeciesID).size());
        strBuf.append(", ");
        strBuf.append(((double) this.freeTFmolecules.get(TFspeciesID).size() / this.TFspecies[TFspeciesID].copyNumber));
        strBuf.append(", ");
        strBuf.append(dna.effectiveTFavailabilitySum[TFspeciesID]);
        strBuf.append(", ");
        strBuf.append(((double) dna.effectiveTFavailabilitySum[TFspeciesID] / dna.effectiveTFavailabilityMaxSum[TFspeciesID]));
        if (TFmoleculesIDs.get(TFspeciesID).size() <= Constants.MAXIMUM_NUMBER_OF_TF_MOLECULES_TO_FOLLOW) {
            for (int i = 0; i < TFmoleculesIDs.get(TFspeciesID).size(); i++) {
                pos = dbp[TFmoleculesIDs.get(TFspeciesID).get(i)].getPosition();
                if (pos != Constants.NONE) {
                    strBuf.append(", ");
                    strBuf.append(pos + dna.subsequence.start);
                } else {
                    strBuf.append(", ");
                    strBuf.append(Constants.NONE);
                }
            }
        }
        strBuf.append("\n");
        printTFtoFollowToFile(TFspeciesID, strBuf.toString());
    }


    private void printTFtoFollowToFile(int TFspeciesID, String str) throws FileNotFoundException {

        if (this.TFstoFollowFiles.containsKey(TFspeciesID)) {

            BufferedWriter bufferFile = null;
            try {
                //Construct the BufferedWriter object
                if (this.outputPath.isEmpty()) {
                    bufferFile = new BufferedWriter(new FileWriter(this.TFstoFollowFiles.get(TFspeciesID), true));
                } else {
                    bufferFile = new BufferedWriter(new FileWriter(new File(this.outputPath,
                            this.TFstoFollowFiles.get(TFspeciesID)), true));
                }
                bufferFile.write(str);
                //bufferFile.newLine();
                bufferFile.flush();
                bufferFile.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }


    /**
     * prints where the output was saved
     */
    private void printOutputFiles() {

        this.printDebugInfo("The model parameters were saved in the file: " + this.outputParamsFile);
        this.printDebugInfo("The status mesages were saved in the file: " + this.outputStatusFile);
        if (this.dna.isRandom || this.dna.useSubSequence) {
            this.printDebugInfo("The DNA sequence was saved in the file: " + this.outputDNASequenceFile);
        }
        if (this.ip.OUTPUT_AFFINITY_LANDSCAPE.value) {
            this.printDebugInfo("The DNA affinity landscape was saved in the file: " + this.outputAffinityLandscapeFile);
        }
        if (this.ip.OUTPUT_DNA_OCCUPANCY.value) {
            this.printDebugInfo("The DNA occupancy was saved in the file: " + this.outputDNAOccupancyFile);
        }
        this.printDebugInfo("The target site information was saved in the file: " + this.outputTargetSiteFile);
        if (this.ip.FOLLOW_TS.value) {
            this.printDebugInfo("The TF information was saved in the file: " + this.outputTFFile);
        }

        if (this.areTFstoFollow) {
            for (int i : this.TFstoFollowID) {
                this.printDebugInfo("The information about the following TFs (" + this.TFstoFollow.get(i)
                        + ") was saved in file: " + this.TFstoFollowFiles.get(i));
            }
        }
    }

    /**
     * prints to the debug file or screen the debug message
     *
     * @param str the string to print
     */
    public void printDebugInfo(String str) {
        printDebugInfo(str, true);
    }

    /**
     * prints to the debug file or screen the debug message
     *
     * @param str
     * @param append
     */
    public void printDebugInfo(String str, boolean append) {
        BufferedWriter statusBuffer = null;
        try {
            //Construct the BufferedWriter object
            if (this.outputPath.isEmpty()) {
                statusBuffer = new BufferedWriter(new FileWriter(this.outputStatusFile, true));
            } else {
                statusBuffer = new BufferedWriter(new FileWriter(
                        new File(this.outputPath, this.outputStatusFile), true));
            }

            statusBuffer.write(str);
            statusBuffer.newLine();
            statusBuffer.flush();
            statusBuffer.close();

        } catch (IOException ex) {
            ex.printStackTrace();
        }

        //if there is a gui then write the line also to the gui
        if (gui != null) {
            gui.printlnStatusArea(str);
        }

    }

    /**
     * prints to the debug file or screen the debug message
     */
    private void printTargetSiteToFollowInfo() {

        if (this.ip.ENSEMBLE_SIZE.value == 1 && this.ip.FOLLOW_TS.value) {

            BufferedWriter statusBuffer = null;
            try {
                //Construct the BufferedWriter object
                if (this.outputPath.isEmpty()) {
                    statusBuffer = new BufferedWriter(new FileWriter(this.outputTargetSiteFollowFile, true));
                } else {
                    statusBuffer = new BufferedWriter(new FileWriter(
                            new File(this.outputPath, this.outputTargetSiteFollowFile), true));
                }
                for (int i = 0; i < (this.TargetSiteFollowLines.size() - 1); i++) {
                    statusBuffer.write(this.TargetSiteFollowLines.get(i));
                    statusBuffer.newLine();
                }
                statusBuffer.flush();
                statusBuffer.close();
                this.TargetSiteFollowLines.clear();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    /**
     * run the simulations for a certain amount of time only.
     *
     * @param stopTime
     * @throws FileNotFoundException
     */
    public double runInterval(double stopTime, double totalElapsedTime) throws FileNotFoundException {
        this.ip.STOP_TIME.value = stopTime;
        this.totalElapsedTime = totalElapsedTime;

        //restart random number generator
        if (this.cellTime == 0) {
            createRandomNumberGenerator();
        }

        this.isPartialSimulation = (this.ip.STOP_TIME.value != this.totalStopTime)
                || (this.ip.STOP_TIME.value == this.totalStopTime && this.cellTime != 0);

        if (this.isPartialSimulation && !this.isInDebugMode() && this.cellTime == 0 && this.ensemble == 0) {
            printDebugInfo("#sample, cellTime, elapsedTimeSec");
        }

        return run();
    }

    /**
     * run the simulations for a number of time steps
     */
    public double run() throws FileNotFoundException {

        boolean hasNextEvent = true;
        long curTime = System.currentTimeMillis();
        double nextPercent = 10;
        double elapsedTimeSec = 0;

        while (!isEndOfPartialSimulation(this.cellTime) && hasNextEvent) {

            //if no TF binding event is scheduled then schedule one
            if (eventQueue.TFBindingEventQueue.isEmpty() && this.freeTFmoleculesTotal > 0) {
                eventQueue.scheduleNextTFBindingEvent(this, this.cellTime);
            }

            hasNextEvent = executeNextEvent(true);

            if (this.areTFstoFollow && this.nextPointTFstoFollow <= cellTime) {
                printTFtoFollowInformation(cellTime);
                nextPointTFstoFollow += stepPointTFstoFollow;
            }

            // print intermediary simulation time points
            if (!this.isPartialSimulation && !this.isInDebugMode()
                    && nextPercent < 100 * this.cellTime / this.totalStopTime)
            {
                // compute real time of simulation
                elapsedTimeSec = this.computeElapsedTime(curTime);
                this.printDebugInfo(nextPercent + "% finished in " + elapsedTimeSec + " seconds");
                nextPercent = nextPercent + 10;
            }

            //print intermediary steady state results
            if (this.ip.PRINT_INTERMEDIARY_RESULTS_AFTER.value > 0 && !isEndOfSimulation(this.cellTime)
                    && this.cellTime >= this.lastPrintResultsAfter + this.ip.PRINT_INTERMEDIARY_RESULTS_AFTER.value
                    && this.totalStopTime > this.lastPrintResultsAfter + 2 * this.ip.PRINT_INTERMEDIARY_RESULTS_AFTER.value)
            {
                this.lastPrintResultsAfter += this.ip.PRINT_INTERMEDIARY_RESULTS_AFTER.value;
                printSteadyStates(this.lastPrintResultsAfter);
            }

        }
        elapsedTimeSec = this.computeElapsedTime(curTime);

        printFinalInfo(elapsedTimeSec);

        //print steady state info
        if (isEndOfSimulation(this.cellTime)) {
            performEndSampleActions(curTime);
            if (ensemble >= this.ip.ENSEMBLE_SIZE.value) {
                performEndActions(elapsedTimeSec);
                printFinalDebugInfo(elapsedTimeSec);
            }
        }

        return elapsedTimeSec;
    }

    public boolean isEndOfSimulation(double time) {
        double doubleZero = Constants.DOUBLE_ZERO * this.totalStopTime;
        return (this.totalStopTime - time < doubleZero);
    }

    private boolean isEndOfPartialSimulation(double time) {
        double doubleZero = Constants.DOUBLE_ZERO * this.totalStopTime;
        return (this.ip.STOP_TIME.value - time < doubleZero);
    }

    /**
     * performs the actions after each sample finished
     *
     * @param curTime
     */
    private void performEndSampleActions(double curTime) {
        //update the number of simulated samples
        ensemble++;

        //update final occupancy
        if (this.ip.PRINT_FINAL_OCCUPANCY.value) {
            dna.updateFinalPosition(this);
        }

        //update total simulate time
        this.totalSimulatedTime += this.cellTime;

        //update boud time
        this.updateDNABoundTime();

        //reset simulation if not last simulation
        if (this.ensemble < this.ip.ENSEMBLE_SIZE.value) {
            this.resetInternalParameters();
        }
    }

    /**
     * performs the final actions ot the end of the simulations
     *
     * @param elapsedTimeSec
     * @throws FileNotFoundException
     */
    private void performEndActions(double elapsedTimeSec) throws FileNotFoundException {
        //compute the time TFs were bound to DNA.
        this.computeDNABoundTime();

        //record last sliding lengths
        this.recordLastSlidingLengths();

        //if the simulator follows the occupancy of target sites write the last values;
        if (this.ip.FOLLOW_TS.value) {
            this.TargetSiteFollowLines.add(this.cellTime + ", " + this.tsg.getTargetSiteGroupsOccupancyString());
            this.printTargetSiteToFollowInfo();
        }

        //count collisions
        this.dna.collisionsCountTotal = Utils.computeSum(this.dna.collisionsCount);
    }

    /**
     * method that executes the next event in the list
     *
     * @return
     */
    private boolean executeNextEvent(boolean fixStopTime) {
        Event e;
        boolean result = true;
        int proteinID;

        //if no TF binding event is scheduled then schedule one
        if (eventQueue.TFBindingEventQueue.isEmpty() && this.freeTFmoleculesTotal > 0) {
            eventQueue.scheduleNextTFBindingEvent(this, this.cellTime);
        }

        e = eventQueue.getNextEvent();

        if (e != null) {
            eventsCount++;
            // there is at least one event in the queue
            if (e instanceof ProteinEvent) {
                // next event is a random walk event
                ProteinEvent pe = (ProteinEvent) e;
                proteinID = pe.proteinID;

                if (!fixStopTime || e.time <= this.totalStopTime) {
                    this.dbp[proteinID].act(this, pe);
                    this.cellTime = e.time;
                } else {
                    // fixStopTime == True && e.time > this.totalStopTime
                    this.cellTime = this.totalStopTime;
                }
                this.eventQueue.scheduleNextTFOnDNAEvent(this, pe.proteinID, e.time);

            } else if (e instanceof RepressionEvent) {
                // an event caused by a repressor
                RepressionEvent re = (RepressionEvent) e;
                if (!fixStopTime || e.time <= this.totalStopTime) {
                    this.remodeller.act(this, re);
                    this.cellTime = e.time;
                } else {
                    // fixStopTime == True && e.time > this.totalStopTime
                    this.cellTime = this.totalStopTime;
                }
                // prevent scheduling event if this is derepression event scheduled due to repressor unbinding
                if (re.scheduleNextEvent) {
                    this.eventQueue.scheduleNextTFOnDNAEvent(this, re.proteinID, e.time);
                }
            }
        } else {
            if (!canTFMoleculeBind()) {
                printDebugInfo("Hm... no more events at time " + (this.cellTime + this.totalSimulatedTime));
                result = false;
            }
        }

        // FG: print the dna state arrays, if they are not too big
        int strandLength = this.dna.strand.length;
        if (this.isInDebugMode() && strandLength < Constants.MAX_STRAND_LEN_TO_PRINT_STRAND_STATE_ARRAYS) {
            StringBuilder arrayPos      = new StringBuilder("position: ");
            StringBuilder arrayClosed   = new StringBuilder("closed:   ");
            StringBuilder arrayOccupied = new StringBuilder("occupied: ");
            StringBuilder arrayTFavail;
            for (int i = 0; i < strandLength; i++) {
                arrayPos.append(String.format("%2d,", i));
                arrayClosed.append(String.format("%2d,", this.dna.closed[i]));
                arrayOccupied.append(String.format("%2d,", this.dna.occupied[i]));
            }
            printDebugInfo(arrayPos.toString());
            printDebugInfo(arrayClosed.toString());
            printDebugInfo(arrayOccupied.toString());
            int numTF = this.TFspecies.length;
            for (int i_tf = 0; i_tf < numTF; i_tf++) {
                arrayTFavail = new StringBuilder("effAv " + this.TFspecies[i_tf].name + ": ");
                for (int i = 0; i < strandLength; i++) {
                    arrayTFavail.append(String.format("%2d,", (this.dna.effectiveTFavailability[i_tf][i] ? 1 : 0)));
                }
                printDebugInfo(arrayTFavail.toString());
            }
        }

        return result;
    }


    /**
     * run the simulations  until the target sites are reached
     */
    public double runUntilTSReached() throws FileNotFoundException {

        this.runUntilTSReached = true;

        long curTime = System.currentTimeMillis();
        double nextTimeStep = 1;
        double elapsedTimeSec = 0;
        boolean hasNextEvent = true;
        while (this.ensemble < this.ip.ENSEMBLE_SIZE.value && hasNextEvent) {
            while (this.areTargetSitesToBeReached && hasNextEvent) {
                hasNextEvent = this.executeNextEvent(false);
                if (this.areTFstoFollow && this.nextPointTFstoFollow <= cellTime) {
                    printTFtoFollowInformation(cellTime);
                    nextPointTFstoFollow += stepPointTFstoFollow;
                }
                // print intermediary simulation time points
                if (nextTimeStep < this.cellTime) {
                    // compute real time of simulation
                    elapsedTimeSec = this.computeElapsedTime(curTime);
                    nextTimeStep++;
                }
            }

            this.printDebugInfo("sample " + this.ensemble + ", time " + this.cellTime
                    + " simulated in " + elapsedTimeSec + " seconds");
            elapsedTimeSec = this.computeElapsedTime(curTime);
            this.computeDNABoundTime();
            performEndSampleActions(elapsedTimeSec);
        }
        this.dna.collisionsCountTotal = Utils.computeSum(this.dna.collisionsCount);
        printFinalInfo(elapsedTimeSec);

        if (this.cellTime >= this.totalStopTime && this.ensemble >= this.ip.ENSEMBLE_SIZE.value) {
            printFinalDebugInfo(elapsedTimeSec);
        }

        return elapsedTimeSec;
    }


    /**
     * compute elapsed time
     */
    public double computeElapsedTime(double curTime) {
        long finalTime = System.currentTimeMillis();
        return (finalTime - curTime) / 1000;
    }


    /**
     * checks whether the system is in debug mode
     */
    public boolean isInDebugMode() {
        return this.ip.DEBUG_MODE.value;
    }

    /**
     * returns the final stop time of the simulations
     */
    public double getTotalStopTime() {
        return this.totalStopTime;
    }


    /**
     * prints initial information
     */
    public void printPreprocesedInfo() {
        int start = 0;
        int end = dna.strand.length;

        //very slow writing
        try {
            if (dna.isRandom || dna.useSubSequence) {
                dna.printDNAstrand(this.outputPath, this.outputDNASequenceFile);
            }
            if (this.ip.OUTPUT_AFFINITY_LANDSCAPE.value) {
                dna.printAffinities(this.outputPath, this.outputAffinityLandscapeFile, start, end,
                        this.TFspecies, this.ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value,
                        this.ip.WIG_STEP.value, this.ip.WIG_THRESHOLD.value, this.ip.OUTPUT_BINDING_ENERGY.value);
            }

            printTFspecies(0, true, this.outputTFFile);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * prints the objects (Target sites/TFs/Ocupancy) to the files
     *
     * @throws FileNotFoundException
     */
    public void printSteadyStates(double time) {
        int start = 0;
        int end = dna.strand.length;
        String filename;

        try {
            //occupancy
            if (this.ip.OUTPUT_DNA_OCCUPANCY.value) {
                filename = outputDNAOccupancyFile;
                if (!isEndOfSimulation(time)) {
                    filename = filename.replaceAll("occupancy", "occupancy_" + time + "s");
                }
                dna.printDNAoccupancy(this.outputPath, filename, start, end,
                        this.ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value, this.ip.WIG_STEP.value,
                        this.ip.WIG_THRESHOLD.value);
            }

            //print current occupancy
            if (this.ip.PRINT_FINAL_OCCUPANCY.value) {
                filename = outputDNAOccupancyFinalFile;
                if (!isEndOfSimulation(time)) {
                    filename = filename.replaceAll("occupancy_final", "occupancy_final_" + time + "s");
                }
                dna.printFinalPosition(this.outputPath, filename, start, end,
                        this.ip.DNA_OCCUPANCY_FULL_MOLECULE_SIZE.value, this.ip.WIG_STEP.value,
                        this.ip.WIG_THRESHOLD.value);
            }

            //print the sliding lengths
            if (this.ip.OUTPUT_SLIDING_LENGTHS.value) {
                for (int i = 0; i < TFspecies.length; i++) {
                    filename = this.outputParamsFile.getName().replaceAll("params", TFspecies[i].name +
                            "_sliding_lengths").replaceAll(Constants.PARAMETR_FILE_EXTENSION,
                            Constants.SLIDING_LENGTHS_FILE_EXTENSION);
                    if (!isEndOfSimulation(time)) {
                        filename = filename.replaceAll("sliding_lengths", "sliding_lengths_" + time + "s");
                    }
                    TFspecies[i].printSlidingLengths(this.outputPath, filename);
                    this.printDebugInfo("sliding lengths for " + TFspecies[i].name + " were printed in file " + filename);
                    filename = this.outputParamsFile.getName().replaceAll("params", TFspecies[i].name +
                            "_observed_sliding_lengths").replaceAll(Constants.PARAMETR_FILE_EXTENSION,
                            Constants.SLIDING_LENGTHS_FILE_EXTENSION);
                    if (!isEndOfSimulation(time)) {
                        filename = filename.replaceAll("observed_sliding_lengths",
                                "observed_sliding_lengths_" + time + "s");
                    }
                    TFspecies[i].printObservedSlidingLengths(this.outputPath, filename);
                    this.printDebugInfo("observed sliding lengths for " + TFspecies[i].name + " were printed in file "
                            + filename);
                }
            }

            // FG: print the repressed lengths
            if (this.ip.OUTPUT_REPRESSED_LENGTHS.value) {
                filename = this.outputParamsFile.getName().replaceAll("params",
                        "repressed_lengths").replaceAll(Constants.PARAMETR_FILE_EXTENSION,
                        Constants.SLIDING_LENGTHS_FILE_EXTENSION);
                if (!isEndOfSimulation(time)) {
                    filename = filename.replaceAll("repressed_lengths", "repressed_lengths_" + time + "s");
                }
                this.printRepressedLengths(this.outputPath, filename);
                this.printDebugInfo("repressed lengths were printed in file " + filename);
            }

            //Target sites
            filename = this.outputTargetSiteFile;
            if (!isEndOfSimulation(time)) {
                filename = filename.replaceAll("target_site", "target_site_" + time + "s");
            }
            printTargetSitesInformation();
            printTFspecies(time, false, filename);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /** FG
     * Print the repressed data to the file
     */
    public void printRepressedLengths(String path, String filename) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path + filename));
            // header
            out.write("\"time\", \"repressedLength\", \"repressedRepScore\", \"repressedActScore\" \n");
            // data
            for (int i = 0; i < this.dna.repressionData.size(); i++) {
                RepressionData repData = this.dna.repressionData.get(i);
                out.write(String.format("%f, %d, %g, %g", repData.getTime(), repData.getRepressedLength(),
                        repData.getRepressedRepScore(), repData.getRepressedActScore()));
                out.newLine();
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void printTFspecies(double time, boolean reduced, String filename) throws IOException {
        //TF species
        if (TFspecies.length > 0) {
            BufferedWriter bufferFile = null;

            filename = this.outputTFFile;
            if (!isEndOfSimulation(time)) {
                filename = filename.replaceAll("TF_species", "TF_species_" + time + "s");
            }
            //Construct the BufferedWriter object
            if (this.outputPath.isEmpty()) {
                bufferFile = new BufferedWriter(new FileWriter(filename));
            } else {
                bufferFile = new BufferedWriter(new FileWriter(new File(this.outputPath, filename)));
            }
            bufferFile.write(TFspecies[0].headerToString(reduced));
            bufferFile.newLine();
            for (TFSpecies tFspecy : TFspecies) {
                bufferFile.write(tFspecy.toString(this, reduced));
                bufferFile.newLine();
            }
            bufferFile.flush();
            bufferFile.close();

        }
    }

    /**
     * records last sliding lengths
     */
    public void recordLastSlidingLengths() {
        if (this.ip.OUTPUT_SLIDING_LENGTHS.value) {
            for (DBP value : this.dbp) {
                if (value.getPosition() != Constants.NONE) {
                    TFspecies[value.speciesID].slidingLength.add(value.getSlidingLength());
                    TFspecies[value.speciesID].slidingEvents.add(value.getSlidingEvents());
                    TFspecies[value.speciesID].observedSlidingLength.add(value.getObservedSlidingLength());
                }
            }
        }
    }

    /**
     * returns the number of TF species
     */
    public int getNoOfDBPspecies() {
        return this.TFspecies.length;
    }


    /**
     * prints target site information
     *
     * @throws FileNotFoundException
     */
    public void printTargetSitesInformation() {

        if (this.ip.ENSEMBLE_SIZE.value == 1 && this.ip.FOLLOW_TS.value) {
            BufferedWriter bufferFile = null;
            try {
                //Construct the BufferedWriter object
                if (this.outputPath.isEmpty()) {
                    bufferFile = new BufferedWriter(new FileWriter(this.outputTargetSiteFile));
                } else {
                    bufferFile = new BufferedWriter(new FileWriter(new File(this.outputPath, this.outputTargetSiteFile)));
                }
                //header
                String str = "\"targetSite\", \"firstReached\", \"timesReached\", \"timeOccupied\"";
                bufferFile.write(str);
                bufferFile.newLine();
                bufferFile.write(tsg.toString());
                bufferFile.flush();
                bufferFile.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

    }

    /**
     * prints final informations
     */
    private void printFinalInfo(double elapsedTimeSec) {
        if (!this.isInDebugMode()) {
            if (this.isPartialSimulation) {
                String str = this.ensemble + ", " + this.cellTime + ", ";
                str += (this.totalElapsedTime + elapsedTimeSec) + ", ";
                printDebugInfo(str);
            }
        }
    }

    /**
     * the final debug info to be printed
     */
    private void printFinalDebugInfo(double elapsedTimeSec) {
        printDebugInfo("=====================================================");
        printDebugInfo("Elapsed time: " + Utils.roundTwoDecimals(this.totalElapsedTime + elapsedTimeSec) + " sec");
        printDebugInfo("=====================================================");
        printDebugInfo("Collisions: " + this.dna.collisionsCountTotal + " ");
        printDebugInfo("=====================================================");

        this.printDebugInfo("The simulator performed " + this.eventsCount + " events in "
                + Utils.roundTwoDecimals(this.totalElapsedTime + elapsedTimeSec) + " seconds, which means it performs "
                + Utils.roundTwoDecimals(this.eventsCount / (this.totalElapsedTime + elapsedTimeSec)) + " events/s ");

        this.printSteadyStates(this.cellTime);
        this.printOutputFiles();

    }


    /**
     * returns the last free molecule of a specific TF species
     *
     * @param TFspeciesID the TF specie
     * @return -1 if nothing is found
     */
    public int getFreeTFmolecule(int TFspeciesID) {
        int result = Constants.NONE;
        if (!this.freeTFmolecules.get(TFspeciesID).isEmpty()) {
            result = this.freeTFmolecules.get(TFspeciesID).get(this.freeTFmolecules.get(TFspeciesID).size() - 1);
        }
        return result;
    }


    /**
     * checkes whether any TF molecule that is free can bind
     */
    public boolean canTFMoleculeBind() {
        boolean result = false;
        for (int i = 0; i < this.freeTFmolecules.size() && !result; i++) {
            if (this.freeTFmolecules.get(i).size() > 0 && dna.effectiveTFavailabilitySum[i] > 0) {
                result = true;
            }
        }
        return result;
    }

    //updateDNA bound time
    public void updateDNABoundTime() {
        double timeBound;

        //update the bound time for all molecules that where bound
        for (DBP value : dbp) {
            if (value.getPosition() != Constants.NONE && value.getTimeOfLastPositionChange() < this.cellTime) {
                timeBound = this.cellTime - value.getTimeOfLastPositionChange();
                value.updateBoundTime(this, timeBound, value.getDirection(), value.getPosition());

                //update target sites statistics
                int tsID = dna.isTargetSite[value.speciesID][value.getPosition()][value.getDirection()];
                if (tsID != Constants.NONE) {
                    tsg.updateTargetSiteStatistics(tsID, this.cellTime, true);
                }
            }
        }
    }

    /**
     * method that computes the average bound time
     */
    public void computeDNABoundTime() {
        //TF species bound time
        for (int j = 0; j < TFspecies.length; j++) {
            TFspecies[j].timeBoundAvg = 0;
            for (int i = 0; i < this.dna.strand.length; i++) {
                for (int dir = 0; dir < this.TFreadingDirection; dir++) {
                    TFspecies[j].timeBoundAvg += dna.effectiveTFOccupancy[j][i][dir];
                }
            }
            TFspecies[j].timeBoundAvg /= (TFspecies[j].copyNumber * this.totalSimulatedTime);
            TFspecies[j].slidingEventsPerBinding =
                    ((double) (TFspecies[j].countTFSlideLeftEvents + TFspecies[j].countTFSlideRightEvents)
                            / (TFspecies[j].countTFHoppingEvents + TFspecies[j].countTFBindingEvents));
            TFspecies[j].slidingLengthPerBinding = Math.sqrt(2 * TFspecies[j].slidingEventsPerBinding);
            TFspecies[j].observedSlidingLengthPerBinding =
                    Math.sqrt(2 * (double) (TFspecies[j].countTFSlideLeftEvents + TFspecies[j].countTFSlideRightEvents
                            + TFspecies[j].countTFHoppingEvents) / (TFspecies[j].countTFBindingEvents));
            TFspecies[j].residenceTimePerBinding =
                    ((TFspecies[j].timeBoundAvg * this.totalSimulatedTime)
                            / ((double) TFspecies[j].countTFBindingEvents / TFspecies[j].copyNumber));
        }

    }


    /**
     * Always treat de-serialization as a full-blown constructor, by
     * validating the final state of the de-serialized object.
     */
    private void readObject(ObjectInputStream in) throws ClassNotFoundException, IOException {
        //always perform the default de-serialization first
        in.defaultReadObject();
    }

    /**
     * This is the default implementation of writeObject.
     * Customise if necessary.
     */
    private void writeObject(ObjectOutputStream out) throws IOException {
        //perform the default serialization for all non-transient, non-static fields
        out.defaultWriteObject();
    }


    /**
     * stops the simulation and prints an error message
     *
     * @param str the error message to print
     */
    public void stopSimulation(String str) {
        this.printDebugInfo(str);
        System.err.println(str);
        try {
            System.exit(0);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * updates statistics about a target site
     */
    public void updateTargetSiteStatistics(int tsID, double time, boolean bound) {

        this.tsg.updateTargetSiteStatistics(tsID, time, bound);

        //add a line in the target sites follow file with the new occupancy
        if (this.ip.FOLLOW_TS.value) {
            String str = time + ", ";
            str += this.tsg.getTargetSiteGroupsOccupancyString();
            this.TargetSiteFollowLines.add(str);
            if (this.TargetSiteFollowLines.size() >= Constants.MAX_STRING_LENGTH) {
                this.printTargetSiteToFollowInfo();
            }
        }
    }

}
