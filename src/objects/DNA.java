package objects;

import environment.Cell;
import utils.*;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
 * a class that stores information on DNA strand and affinities of TFs for DNA
 *
 * @author n.r.zabet@gen.cam.ac.uk (original contributor)
 * @author fedor.garbuzov@mail.ioffe.ru (modified the code)
 */
public class DNA implements Serializable {
    /**
     *
     */
    private static final long serialVersionUID = 4312661038245046220L;
    public byte[] strand;
    public int[] occupied;
    public byte[] closed;

    public boolean isRandom;

    public int DNAsectorsCount;
    public int DNAsectorSize;
    public int[] sectorID;

    public int[] TFsize;
    public double[][][] TFavgMoveRate; //TFSpecies; position; direction
    public boolean[][] effectiveTFavailability; // speciesID, position
    public int[][] effectiveTFsectorsAvailabilitySum; //species id, sectors
    public int[] effectiveTFavailabilitySum;
    public int[] effectiveTFavailabilityMaxSum;
    public double[][][] effectiveTFOccupancy; //TFSpecies; position; direction
    public int[][][] finalTFOccupancy; //TFSpecies; position; direction
    public double[][][] firstReached; //TFSpecies; position; direction

    public double[][][] TFSlideLeftNo; //TFSpecies; position; direction
    public double[][][] TFSlideRightNo; //TFSpecies; position; direction


    public double[] bpFreq;

    public int[][][] isTargetSite;
    public boolean areTargetSites;
    public int TFdirections;
    public ArrayList<TargetSite> ts;
    public String description;
    public DNAregion region;
    public boolean useSubSequence;
    public DNAregion subsequence;
    public int copyNumber;
    public boolean isReflexive;
    public boolean isAbsorbing;
    public boolean isPeriodic;
    public HashMap<Integer, Integer> TFidPos;
    public HashMap<Integer, Integer> TFposID;
    public HashMap<Integer, String> TFposName;
    public int[] collisionsCount;
    public long collisionsCountTotal;

    public Integer currentRepressedLength;
    public ArrayList<RepressionData> repressionData; // FG: time, length, repressed repressor score, repressed act. len
    public double repressedRepScore; // FG: sum of the repressed repressor sites scores
    public double repressedActScore; // FG: sum of the repressed repressor sites scores

    public DNA() {
        strand = new byte[0];
        region = new DNAregion("", 0, 0);
        occupied = new int[0];
        closed = new byte[0];
        isRandom = true;
    }

    /**
     * class constructor - generates everything random
     */
    public DNA(Random generator, int length, double proportionOfA, double proportionOfT, double proportionOfC,
               double proportionOfG, String boundaryCondition)
    {
        strand = CellUtils.generateRandomDNASequence(generator, length, proportionOfA, proportionOfT, proportionOfC,
                proportionOfG);
        computeBPfreq();
        parseDescription("randomly generated", length, boundaryCondition);

        TFidPos = new HashMap<Integer, Integer>();
        TFposID = new HashMap<Integer, Integer>();
        TFposName = new HashMap<Integer, String>();

        ts = new ArrayList<TargetSite>();
        setOccupancyAndClosenessVectorsFree();
        isRandom = true;

        currentRepressedLength = 0;
        repressionData = new ArrayList<>();
        repressedRepScore = 0.0;
        repressedActScore = 0.0;
    }


    /**
     * class constructor
     *
     * @param dna the DNA sequence used for initialization
     */
    public DNA(DNA dna) {
        parseDescription(dna.description, dna.strand.length, "");
        int start = (int) (this.subsequence.start - this.region.start);
        int end = (int) Math.min(Math.min(this.subsequence.end - this.region.start,
                this.region.end - this.region.start), dna.strand.length + start);

        this.strand = new byte[end - start];
        if (end - start >= 0) System.arraycopy(dna.strand, start, this.strand, 0, end - start);
        computeBPfreq();

        isRandom = false;

        TFidPos = new HashMap<Integer, Integer>();
        TFposID = new HashMap<Integer, Integer>();
        TFposName = new HashMap<Integer, String>();

        ts = new ArrayList<TargetSite>();
        setOccupancyAndClosenessVectorsFree();

        currentRepressedLength = 0;
        repressionData = new ArrayList<>();
        repressedRepScore = 0.0;
        repressedActScore = 0.0;
    }

    /**
     * this method parses the description field in a fasta file
     * The description of the fasta file starts with ">" character.
     * This its followed by the name of the chromosome or strain and the position where this sequence was extracted,
     * e.g. MG1655:0..4639675
     * optionally the user can specify the following parameters delimited by ";" character
     * - the subsequence form this fasta file to be considered, e.g. subsequence = "0..100000";
     * - this DNA copy number, e.g. copy = 1;
     * - what boundary condition is used (periodic, reflexive, absorbing), e.g. boundary = "absorbing"
     */
    public void parseDescription(String description, int length, String boundaryCondition) {
        description = description.trim();
        this.description = description;
        String[] buffer = new String[1];
        String[] params;
        DNAregion region = new DNAregion("", 0, length);
        DNAregion subsequence;
        int copyNumber = 1;
        this.isReflexive = true;
        this.isAbsorbing = false;
        this.isPeriodic = false;

        if (this.description.contains(Constants.DNA_FASTA_DELIMITER)) {
            buffer = this.description.split(Constants.DNA_FASTA_DELIMITER);
        } else if (!description.isEmpty()) {
            buffer[0] = description;
        }

        if (buffer.length > 0 && !buffer[0].isEmpty()) {
            region = new DNAregion(buffer[0], "", 0, length, false, false);
        }
        subsequence = new DNAregion("", region.start, region.end);

        if (buffer.length > 1) {
            for (int i = 1; i < buffer.length; i++) {
                params = Utils.extractParameterFromCommandLine(buffer[i] + ";", Constants.PARAMS_FILE_ASSIGNMENT_CHAR);
                if (params.length == 2) {
                    if (params[0].equalsIgnoreCase(Constants.DNA_FASTA_SUBSEQUENCE) && !params[1].isEmpty()) {
                        subsequence = new DNAregion(params[1], "", 0, length, false, false);
                    } else if (params[0].equalsIgnoreCase(Constants.DNA_FASTA_COPY_NUMBER) && !params[1].isEmpty()) {
                        copyNumber = Utils.parseInteger(params[1], copyNumber);
                    } else if (params[0].equalsIgnoreCase(Constants.DNA_FASTA_BOUNDARY) && !params[1].isEmpty()) {
                        setBoundaryCondition(params[1]);
                    }
                }
            }
        }

        if (!boundaryCondition.isEmpty()) {
            setBoundaryCondition(boundaryCondition);
        }

        this.region = region;
        this.subsequence = subsequence;
        this.useSubSequence = subsequence.isSubSequenceOf(region);
        this.copyNumber = copyNumber;

    }

    /**
     * sets the boundary condition based on the string parsed as the argument (absorbing/reflexive/periodic)
     */
    private void setBoundaryCondition(String boundaryCondition) {
        if (boundaryCondition.equalsIgnoreCase(Constants.DNA_FASTA_BOUNDARY_ABSORBING)) {
            isReflexive = false;
            isAbsorbing = true;
            isPeriodic = false;
        } else if (boundaryCondition.equalsIgnoreCase(Constants.DNA_FASTA_BOUNDARY_PERIODIC)) {
            isReflexive = false;
            isAbsorbing = false;
            isPeriodic = true;
        } else if (boundaryCondition.equalsIgnoreCase(Constants.DNA_FASTA_BOUNDARY_REFLEXIVE)) {
            isReflexive = true;
            isAbsorbing = false;
            isPeriodic = false;
        }

    }

    /**
     * generates a description string for the DNA
     */
    public String generateDescriptionString() {
        StringBuilder strBuff = new StringBuilder();
        if (useSubSequence) {
            strBuff.append(this.subsequence.toString());
        } else {
            strBuff.append(this.region.toString());
        }
        if (this.copyNumber != 1) {
            strBuff.append(Constants.DNA_FASTA_DELIMITER);
            strBuff.append(Constants.DNA_FASTA_COPY_NUMBER);
            strBuff.append(Constants.PARAMS_FILE_ASSIGNMENT_CHAR);
            strBuff.append(copyNumber);
        }
        if (this.isReflexive || this.isAbsorbing || this.isPeriodic) {
            strBuff.append(Constants.DNA_FASTA_DELIMITER);
            strBuff.append(Constants.DNA_FASTA_BOUNDARY);
            strBuff.append(Constants.PARAMS_FILE_ASSIGNMENT_CHAR);
            if (this.isReflexive) {
                strBuff.append(Constants.DNA_FASTA_BOUNDARY_REFLEXIVE);
            } else if (this.isAbsorbing) {
                strBuff.append(Constants.DNA_FASTA_BOUNDARY_ABSORBING);
            } else {
                strBuff.append(Constants.DNA_FASTA_BOUNDARY_PERIODIC);
            }
        }
        return strBuff.toString();
    }

    /**
     * initiates the occupancy vector with -1 which stands for free region
     */
    private void setOccupancyAndClosenessVectorsFree() {
        this.occupied = new int[strand.length];
        this.closed = new byte[strand.length];
        freeDNA(Constants.FIRST, strand.length);
        openDNA(Constants.FIRST, strand.length);
    }

    /**
     * computes bp frequency
     */
    private void computeBPfreq() {
        long[] countBP = new long[CellUtils.bps.numberOfBP]; // filled with 0 by default

        bpFreq = new double[CellUtils.bps.numberOfBP];
        for (int i = 0; i < countBP.length; i++) {
            bpFreq[i] = 0;
        }

        if (strand.length > 0) {
            //count
            for (int i = 0; i < strand.length; i++) {
                if (strand[i] > 3) {
                    System.out.println("Error: unknown bp " + strand[i] + " " + i);
                }
                countBP[strand[i]]++;
            }
            //compute frequency
            for (int i = 0; i < countBP.length; i++) {
                bpFreq[i] = (double) countBP[i] / strand.length;
            }
        }
    }


    /**
     * sets the sequence provided as a parameter as the current one
     */
    public void loadSequence(ArrayList<Byte> seq) {
        this.strand = new byte[seq.size()];
        for (int i = 0; i < seq.size(); i++) {
            this.strand[i] = seq.get(i);
        }
    }


    /**
     * FG
     */
    private void recomputeTFAffinityLandscapeForClosedRegions(int startPos, int endPos, int speciesID) {
        boolean stop = false;
        int start = startPos, end;
        while (!stop) {
            while (start < endPos && closed[start] == Constants.BP_IS_OPEN) {
                start++;
            }
            if (start < endPos) {
                end = start + 1;
                while (end < endPos && closed[end] != Constants.BP_IS_OPEN) {
                    end++;
                }
                closeRegionInAffinityLandscape(start, end, speciesID);
                start = end;
            } else {
                stop = true;
            }
        }
    }

    /**
     * FG
     */
    private void recomputeTFAffinityLandscapeForRepressedRegions(int startPos, int endPos, int speciesID) {
        boolean stop = false;
        int start = startPos, end;
        while (!stop) {
            while (start < endPos && closed[start] != Constants.BP_IS_REPRESSED) {
                start++;
            }
            if (start < endPos) {
                end = start + 1;
                while (end < endPos && closed[end] == Constants.BP_IS_REPRESSED) {
                    end++;
                }
                closeRegionInAffinityLandscape(start, end, speciesID);
                start = end;
            } else {
                stop = true;
            }
        }
    }

    /**
     * compute the affinity landscape for a list of TF species
     */
    public void computeTFaffinityLandscape(Cell n, Random generator, TFSpecies[] TFspecies, //double TFSpecificWaitingTime,
                                           int TFdirections, int DNAsectorSize, boolean printFinalOccupancy)
    {
        if (TFspecies != null && TFspecies.length > 0) {
            TFsize = new int[TFspecies.length];
            double[][] TFaffinitiesLR = new double[TFspecies.length][strand.length];
            double[][] TFaffinitiesRL = new double[TFspecies.length][strand.length];
            effectiveTFavailability = new boolean[TFspecies.length][strand.length];
            effectiveTFavailabilitySum = new int[TFspecies.length];
            effectiveTFavailabilityMaxSum = new int[TFspecies.length];
            TFavgMoveRate = new double[TFspecies.length][strand.length][TFdirections];
            effectiveTFOccupancy = new double[TFspecies.length][strand.length][TFdirections];
            this.TFSlideLeftNo = new double[TFspecies.length][strand.length][TFdirections];
            this.TFSlideRightNo = new double[TFspecies.length][strand.length][TFdirections];


            if (printFinalOccupancy) {
                this.finalTFOccupancy = new int[TFspecies.length][strand.length][TFdirections];
                for (int i = 0; i < TFspecies.length; i++) {
                    for (int j = 0; j < strand.length; j++) {
                        for (int k = 0; k < TFdirections; k++) {
                            this.finalTFOccupancy[i][j][k] = 0;
                        }
                    }
                }

            }

            this.TFdirections = TFdirections;
            isTargetSite = new int[TFspecies.length][strand.length][TFdirections];
            //this.TFSpecificWaitingTime = TFSpecificWaitingTime;
            firstReached = new double[TFspecies.length][strand.length][TFdirections];

            this.collisionsCount = new int[strand.length];
            for (int i = 0; i < strand.length; i++) {
                this.collisionsCount[i] = 0;
            }

            //sectors
            if (DNAsectorSize == 0) {
                DNAsectorSize = (int) Math.floor(Math.sqrt(strand.length));
            } else if (DNAsectorSize > strand.length || DNAsectorSize < 0) {
                DNAsectorSize = strand.length;
            }
            this.DNAsectorSize = DNAsectorSize;
            this.DNAsectorsCount = (int) Math.ceil((double) strand.length / DNAsectorSize);
            this.effectiveTFsectorsAvailabilitySum = new int[TFspecies.length][this.DNAsectorsCount];
            sectorID = new int[strand.length];
            for (int i = 0; i < strand.length; i++) {
                sectorID[i] = i / this.DNAsectorSize;
            }

            for (int i = 0; i < TFspecies.length; i++) {
                TFsize[i] = TFspecies[i].sizeTotal;
                //read the affinities LR
                if (TFspecies[i].landscapeFile != null && !TFspecies[i].landscapeFile.isEmpty()) {
                    //read from a file
                    try {
                        ArrayList<String> lines = Utils.readLinesFromFile(TFspecies[i].landscapeFile);
                        //if not enough lines in the affinity file then stop the simulation
                        if ((lines.size() - TFspecies[i].landscapeEscapeLines) < strand.length) {
                            n.stopSimulation("error while parsing the affinity landscape file "
                                    + TFspecies[i].landscapeFile + ": insufficient entries in the affinity file");
                        }

                        long currentPos;
                        int actualPos = Constants.NONE;
                        double affinityLR, affinityRL;
                        double[] bufferValues;

                        int maxCol = Math.max(TFspecies[i].landscapeAffinityColLR, TFspecies[i].landscapeAffinityColRL);
                        maxCol = Math.max(maxCol, TFspecies[i].landscapePosCol);

                        for (int k = TFspecies[i].landscapeEscapeLines; k < lines.size() && actualPos < this.strand.length; k++) {
                            bufferValues = Utils.parseCSVline(lines.get(k), Constants.AFFINITY_CSV_FILE_DELIMITER, Constants.NONE);
                            //Utils.containsValue(bufferValues, Constants.NONE) ||
                            if (bufferValues == null || bufferValues.length <= maxCol) {
                                n.stopSimulation("error while parsing the affinity landscape file "
                                        + TFspecies[i].landscapeFile + ": could not parse line " + i + ": " + lines.get(k));
                            }
                            assert bufferValues != null;
                            currentPos = (int) Math.round(bufferValues[TFspecies[i].landscapePosCol]);
                            actualPos = (int) (currentPos - this.subsequence.start);
                            if (actualPos >= 0 && actualPos < this.strand.length) {
                                affinityLR = bufferValues[TFspecies[i].landscapeAffinityColLR];
                                affinityRL = bufferValues[TFspecies[i].landscapeAffinityColRL];
                                TFaffinitiesLR[i][actualPos] = affinityLR;
                                TFaffinitiesRL[i][actualPos] = affinityRL;
                            }
                        }
                    } catch (IOException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                } else if (TFspecies[i].seqsFile != null && !TFspecies[i].seqsFile.isEmpty()) {
                    // read affinities from sequence file.
                    try {
                        HashMap<DNAsequence, Double> seqsAffinities = new HashMap<DNAsequence, Double>();
                        ArrayList<String> lines = Utils.readLinesFromFile(TFspecies[i].seqsFile);
                        String[] cells;
                        double affinityValue;
                        boolean correct = true;
                        String currentLine = "";
                        int sizeInBP = TFspecies[i].sizeTotal - TFspecies[i].sizeLeft - TFspecies[i].sizeRight;
                        for (int k = TFspecies[i].seqsEscapeLines; k < lines.size() && correct; k++) {
                            currentLine = lines.get(k);
                            cells = currentLine.split(Constants.AFFINITY_CSV_FILE_DELIMITER);

                            if (cells.length == 2) {
                                cells[0] = cells[0].trim();
                                if (cells[0].length() != sizeInBP) {
                                    correct = false;
                                } else {
                                    affinityValue = Utils.parseDouble(cells[1], Constants.NONE);
                                    if (affinityValue == Constants.NONE) {
                                        correct = false;
                                    } else {
                                        seqsAffinities.put(new DNAsequence(cells[0]), affinityValue);
                                    }
                                }
                            } else {
                                correct = false;
                            }
                        }
                        if (!correct) {
                            n.stopSimulation("error while parsing sequences file " + TFspecies[i].seqsFile + " at " +
                                    "line: " + currentLine);
                        }
                        TFaffinitiesLR[i] = CellUtils.computeTFAffinities(strand, seqsAffinities,
                                TFspecies[i].seqsDefaultValue, TFspecies[i].sizeLeft,
                                TFspecies[i].sizeTotal - TFspecies[i].sizeLeft - TFspecies[i].sizeRight,
                                TFspecies[i].sizeTotal, 0);
                        TFaffinitiesRL[i] = CellUtils.computeTFAffinities(strand, seqsAffinities,
                                TFspecies[i].seqsDefaultValue, TFspecies[i].sizeLeft,
                                TFspecies[i].sizeTotal - TFspecies[i].sizeLeft - TFspecies[i].sizeRight,
                                TFspecies[i].sizeTotal, 1);


                    } catch (IOException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }

                } else if (TFspecies[i].pfm != null && TFspecies[i].pfm.motifSize > 0) {
                    //generate from PWM
                    TFaffinitiesLR[i] = CellUtils.computeTFAffinities(generator, strand, TFspecies[i].pfm,
                            TFspecies[i].sizeLeft, TFspecies[i].sizeTotal, TFspecies[i].es, 0,
                            TFspecies[i].affinityLandscapeRoughness);
                } else {
                    //generate from sequence
                    TFaffinitiesLR[i] = CellUtils.computeTFAffinities(generator, strand, TFspecies[i].dbd,
                            TFspecies[i].sizeLeft, TFspecies[i].sizeTotal, TFspecies[i].es, 0,
                            TFspecies[i].affinityLandscapeRoughness);
                }

                effectiveTFavailabilitySum[i] = 0;
                effectiveTFavailabilityMaxSum[i] = 0;
                for (int j = 0; j < this.DNAsectorsCount; j++) {
                    this.effectiveTFsectorsAvailabilitySum[i][j] = 0;
                }

                for (int j = 0; j < strand.length; j++) {
                    if (TFaffinitiesLR[i][j] != Constants.NONE) {
                        effectiveTFavailability[i][j] = true;
                    } else {
                        TFaffinitiesLR[i][j] = 0;
                        effectiveTFavailability[i][j] = false;
                    }

                    TFavgMoveRate[i][j][0] =
                            CellUtils.computeAvgMoveRate(TFspecies[i].specificWaitingTime, -TFaffinitiesLR[i][j]);
                    this.isTargetSite[i][j][0] = Constants.NONE;
                    firstReached[i][j][0] = Constants.NONE;
                    effectiveTFOccupancy[i][j][0] = 0;
                }

                // FG: compute TF availability sum before closing inaccessible regions
                for (int j = 0; j < strand.length; j++) {
                    if (effectiveTFavailability[i][j]) {
                        effectiveTFavailabilitySum[i]++;
                        // increase the number of available spots for sectors position/ sectorSize
                        effectiveTFsectorsAvailabilitySum[i][this.sectorID[j]]++;
                    }
                }

                // FG: recompute affinity landscape due to the closed regions (here, from btrack file)
                recomputeTFAffinityLandscapeForClosedRegions(0, strand.length, i);

                // FG: compute maximal TF availability sum
                for (int j = 0; j < strand.length; j++) {
                    if (effectiveTFavailability[i][j]) {
                        effectiveTFavailabilityMaxSum[i]++;
                    }
                }

                // read the affinities RL
                if (TFdirections == 2) {
                    //read the affinities LR
                    if (TFspecies[i].landscapeFile != null && !TFspecies[i].landscapeFile.isEmpty()) {
                        // do nothing it was already processed
                    } else if (TFspecies[i].pfm != null && TFspecies[i].pfm.motifSize > 0) {
                        TFaffinitiesRL[i] = CellUtils.computeTFAffinities(generator, strand, TFspecies[i].pfm,
                                TFspecies[i].sizeLeft, TFspecies[i].sizeTotal, TFspecies[i].es, 1,
                                TFspecies[i].affinityLandscapeRoughness);
                    } else {
                        TFaffinitiesRL[i] = CellUtils.computeTFAffinities(generator, strand, TFspecies[i].dbd,
                                TFspecies[i].sizeLeft, TFspecies[i].sizeTotal, TFspecies[i].es, 1,
                                TFspecies[i].affinityLandscapeRoughness);
                    }

                    for (int j = 0; j < strand.length; j++) {
                        TFavgMoveRate[i][j][1] =
                                CellUtils.computeAvgMoveRate(TFspecies[i].specificWaitingTime, -TFaffinitiesRL[i][j]);
                        effectiveTFOccupancy[i][j][1] = 0;
                        this.isTargetSite[i][j][1] = Constants.NONE;
                        firstReached[i][j][1] = Constants.NONE;
                    }
                }

                TFidPos.put(TFspecies[i].id, i);
                TFposID.put(i, TFspecies[i].id);
                TFposName.put(i, TFspecies[i].name);

                //slide left and right probabilities
                double intervalLength = TFspecies[i].slideRightNo - TFspecies[i].hopNo;
                double affinityRightLeftRatio;
                for (int j = 0; j < strand.length; j++) {
                    for (int dir = 0; dir < TFdirections; dir++) {
                        //initialise
                        this.TFSlideLeftNo[i][j][dir] = TFspecies[i].slideLeftNo;
                        this.TFSlideRightNo[i][j][dir] = TFspecies[i].slideRightNo;
                        if (j > 0 && j < strand.length - 1 && TFspecies[i].isBiasedRandomWalk) {
                            affinityRightLeftRatio = TFavgMoveRate[i][j - 1][dir] / TFavgMoveRate[i][j + 1][dir];
                            this.TFSlideLeftNo[i][j][dir] = intervalLength / (1 + affinityRightLeftRatio);
                            this.TFSlideRightNo[i][j][dir] =
                                    (affinityRightLeftRatio * intervalLength) / (1 + affinityRightLeftRatio);
                        }
                    }
                }
            }
        }

        this.areTargetSites = false;
        int startRel, endRel, direction;
        if (n.tsg != null && n.tsg.ts != null && !n.tsg.ts.isEmpty()) {
            for (int i = 0; i < n.tsg.ts.size(); i++) {
                startRel = n.tsg.ts.get(i).relStart;
                endRel = n.tsg.ts.get(i).relEnd;
                direction = n.tsg.ts.get(i).region.direction;
                if (direction >= 0 && direction < TFdirections) {
                    for (int k = startRel; k < endRel; k++) {
                        this.isTargetSite[n.tsg.ts.get(i).TFid][k][direction] = n.tsg.ts.get(i).targetSiteID;
                        this.areTargetSites = true;
                    }
                } else {
                    //if there is no direction specified then mark all directions
                    for (int k = startRel; k < endRel; k++) {
                        for (int dir = 0; dir < TFdirections; dir++) {
                            this.isTargetSite[n.tsg.ts.get(i).TFid][k][dir] = n.tsg.ts.get(i).targetSiteID;
                            this.areTargetSites = true;
                        }
                    }
                }
            }

        }
    }


    //to string methods

    /**
     * prints a string with the DNA sequence in letters
     */
    public void printDNAstrand(String path, String filename) {
        CellUtils.printSequence(path, filename, generateDescriptionString(), strand);
    }

    /**
     * prints a string with the DNA  TF affinities
     */
    public void printAffinities(String path, String filename, int start, int end, TFSpecies[] tfs,
                                boolean fullOccupancy, int wigStepSize, double wigThreshold, boolean printBindingEnergy)
    {
        double[][][] bufferTFaffinity = new double[TFavgMoveRate.length][strand.length][this.TFdirections];
        double[][] avg = new double[effectiveTFOccupancy.length][this.TFdirections];

        double[][] cutoff = new double[TFavgMoveRate.length][this.TFdirections];
        int max;

        //normalise the TF occupancy
        for (int dir = 0; dir < TFdirections; dir++) {
            for (int i = 0; i < TFavgMoveRate.length; i++) {
                //init occupancy
                for (int j = 0; j < strand.length; j++) {
                    bufferTFaffinity[i][j][dir] = 0;
                }

                // add occupancy on the entire length of the TFs
                for (int j = 0; j < strand.length; j++) {
                    max = Math.min(strand.length, j + 1);
                    if (fullOccupancy) {
                        max = Math.min(strand.length, j + TFsize[i]);
                    }
                    for (int k = j; k < max; k++) {
                        if (k < strand.length - TFsize[i] + 1) {
                            double tfWaitingTime = 1.0 / TFavgMoveRate[i][j][dir];
                            if (printBindingEnergy) {
                                bufferTFaffinity[i][k][dir] +=
                                        CellUtils.computeBindingEnergy(tfs[i].specificWaitingTime, tfWaitingTime);
                            } else {
                                bufferTFaffinity[i][k][dir] += tfWaitingTime;
                            }
                        }

                    }
                }

                //normalise
                if (!printBindingEnergy) {
                    cutoff[i][dir] = 0;
                    for (int j = 0; j < strand.length; j++) {
                        if (cutoff[i][dir] < bufferTFaffinity[i][j][dir]) {
                            cutoff[i][dir] = bufferTFaffinity[i][j][dir];
                        }
                    }
                    cutoff[i][dir] = cutoff[i][dir] * wigThreshold;
                    avg[i][dir] = 0;
                    for (int j = 0; j < strand.length; j++) {
                        avg[i][dir] += bufferTFaffinity[i][j][dir];
                    }
                    avg[i][dir] /= strand.length;
                }
            }
        }


        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path + filename));

            //check start end output
            if (start < 0 || start >= strand.length) {
                start = 0;
            }

            end = Math.min(end, strand.length);
            if (end < start) {
                end = strand.length;
            }

            double affinity;
            int steps;

            StringBuilder strBuf = new StringBuilder();
            strBuf.append("fixedStep  chrom=");
            strBuf.append(this.region.chromosome);
            strBuf.append("  start=");
            strBuf.append(start);
            strBuf.append("  step=");
            strBuf.append(wigStepSize);
            strBuf.append("  span=");
            strBuf.append(wigStepSize);
            strBuf.append("\n");
            strBuf.append("\"position\"");
            if (this.TFdirections == 1) {
                for (int j = 0; j < TFsize.length; j++) {
                    strBuf.append(", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("\"");
                }
            } else {
                for (int j = 0; j < this.TFsize.length; j++) {
                    strBuf.append(", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("5'3'\", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("3'5'\"");
                }
            }
            out.write(strBuf.toString());
            out.newLine();

            //info
            for (int i = start; i < end; i = i + wigStepSize) {
                strBuf.delete(0, strBuf.length());
                strBuf.append((i + this.subsequence.start));
                for (int j = 0; j < this.TFavgMoveRate.length; j++) {
                    for (int dir = 0; dir < this.TFdirections; dir++) {
                        affinity = 0;
                        steps = 0;
                        for (int k = i; k < Math.min(i + wigStepSize, end); k++) {
                            if ((wigThreshold >= 0 && bufferTFaffinity[j][k][dir] > cutoff[j][dir] && !printBindingEnergy)
                                    || (wigThreshold >= 0 && printBindingEnergy)
                                    || (wigThreshold <= 0 && bufferTFaffinity[j][k][dir] > avg[j][dir])) {
                                affinity += bufferTFaffinity[j][k][dir];
                            }

                            steps++;
                        }
                        strBuf.append(", ");
                        strBuf.append((affinity / steps));
                    }
                }
                out.write(strBuf.toString());
                out.newLine();
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * updates final occupancy
     *
     * @param n the cell
     */
    public void updateFinalPosition(Cell n) {
        for (int i = 0; i < n.dbp.length; i++) {
            if (n.dbp[i].getPosition() != Constants.NONE) {
                this.finalTFOccupancy[n.dbp[i].speciesID][n.dbp[i].getPosition()][n.dbp[i].getDirection()]++;
            }
        }

    }


    /**
     * prints the where TFs are bound to the DNA at the end of the simulation. the molecules are represented by the
     * id of the species and only the first bp is marked as occupied
     */
    public void printFinalPosition(String path, String filename, int start, int end, boolean fullOccupancy,
                                   int wigStepSize, double wigThreshold)
    {
        int[][][] bufferTFoccupancy = new int[finalTFOccupancy.length][strand.length][this.TFdirections];

        double[][] cutoff = new double[finalTFOccupancy.length][this.TFdirections];
        double[][] avg = new double[finalTFOccupancy.length][this.TFdirections];

        int max;
        //normalise the TF occupancy
        for (int dir = 0; dir < TFdirections; dir++) {
            for (int i = 0; i < finalTFOccupancy.length; i++) {
                //init occupancy
                for (int j = 0; j < strand.length; j++) {
                    bufferTFoccupancy[i][j][dir] = 0;
                }

                // add occupancy on the entire length of the TFs
                for (int j = 0; j < strand.length; j++) {
                    max = Math.min(strand.length, j + 1);
                    if (fullOccupancy) {
                        max = Math.min(strand.length, j + TFsize[i]);
                    }
                    for (int k = j; k < max; k++) {
                        bufferTFoccupancy[i][k][dir] += finalTFOccupancy[i][j][dir];

                    }
                }

                //normalise
                cutoff[i][dir] = 0;
                for (int j = 0; j < strand.length; j++) {
                    if (cutoff[i][dir] < bufferTFoccupancy[i][j][dir]) {
                        cutoff[i][dir] = bufferTFoccupancy[i][j][dir];
                    }
                }
                cutoff[i][dir] = cutoff[i][dir] * wigThreshold;


                //computes the threshold of the wig if this is set to autoselect
                avg[i][dir] = 0;
                for (int j = 0; j < strand.length; j++) {
                    avg[i][dir] += bufferTFoccupancy[i][j][dir];
                }
                avg[i][dir] /= strand.length;
            }
        }

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path + filename));

            //check start end output
            if (start < 0 || start >= strand.length) {
                start = 0;
            }

            end = Math.min(end, strand.length);
            if (end < start) {
                end = strand.length;
            }

            //header
            double occupancy;
            int steps;
            StringBuilder strBuf = new StringBuilder();
            strBuf.append("fixedStep  chrom=");
            strBuf.append(this.region.chromosome);
            strBuf.append("  start=");
            strBuf.append(start);
            strBuf.append("  step=");
            strBuf.append(wigStepSize);
            strBuf.append("  span=");
            strBuf.append(wigStepSize);
            strBuf.append("\n");
            strBuf.append("\"position\"");
            if (this.TFdirections == 1) {
                for (int j = 0; j < TFsize.length; j++) {
                    strBuf.append(", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("\"");
                }
            } else {
                for (int j = 0; j < this.TFsize.length; j++) {
                    strBuf.append(", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("5'3'\", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("3'5'\"");
                }
            }
            out.write(strBuf.toString());
            out.newLine();

            //info
            for (int i = start; i < end; i = i + wigStepSize) {
                strBuf.delete(0, strBuf.length());
                strBuf.append((i + this.subsequence.start));

                //TF occupancy
                for (int j = 0; j < bufferTFoccupancy.length; j++) {
                    for (int dir = 0; dir < this.TFdirections; dir++) {
                        occupancy = 0;
                        steps = 0;
                        for (int k = i; k < Math.min(i + wigStepSize, end); k++) {
                            if ((wigThreshold >= 0 && bufferTFoccupancy[j][k][dir] > cutoff[j][dir]) || (wigThreshold <= 0 && bufferTFoccupancy[j][k][dir] > avg[j][dir])) {
                                occupancy += bufferTFoccupancy[j][k][dir];
                            }
                            steps++;
                        }
                        strBuf.append(", ");
                        strBuf.append((occupancy / steps));
                    }
                }
                out.write(strBuf.toString());
                out.newLine();
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * binds a protein to the dna
     *
     * @param proteinID      the ID of the protein to be bound
     * @param position       the position on the dna to bind to
     * @param size           the size of the protein
     * @param checkOccupancy this is true whether this method first checks the occupancy state of the dna
     * @return -1 attempt to bind outside the DNA, proteinID if bound on the DNA, or the ID of the first protein
     * blocking the binding on the DNA
     */
    public int bindMolecule(int proteinID, int position, int size, boolean checkOccupancy) {
        int canBind = canBind(proteinID, position, size, checkOccupancy);
        //bind the protein
        if (canBind == proteinID) {
            occupyDNA(proteinID, position, size);
            recomputeTFAffinityLandscapeOnBinding(position, size);
        }
        return canBind;
    }

    /**
     * checks whether a molecule can bind
     */
    private int canBind(int proteinID, int position, int size, boolean checkOccupancy) {
        int canBind = proteinID;

        //attempt to bind outside the strand interval
        if (position < 0 || position > this.strand.length - size) {
            canBind = Constants.NONE;
        }

        //the DNA is already occupied
        if (canBind == proteinID && checkOccupancy) {
            canBind = getBoundProtein(position, size);
            int bpStatus = checkAvailability(position, size);
            if (canBind == Constants.NONE && bpStatus == Constants.BP_IS_OPEN) {
                canBind = proteinID;
            }
        }
        return canBind;
    }


    /**
     * unbinds a TF molecule  from the dna
     *
     * @param proteinID the ID of the protein to be bound
     * @param position  the position on the dna to bind to
     * @param size      the size of the protein
     * @return -1 attempt to bind outside the DNA, proteinID if bound on the DNA, or the ID of the first protein
     * blocking the binding on the DNA
     */
    public boolean unbindMolecule(Cell n, int proteinID, int position, int size) {
        boolean canUnbind = false;

        //attempt to unbind from outside the strand interval
        if (position >= 0 || position <= this.strand.length - size) {
            canUnbind = true;
            freeDNA(position, size);
            if (!n.dbp[proteinID].isRepressingDNA()) {
                recomputeTFAffinityLandscapeOnUnbinding(position, size);
            }
        }

        return canUnbind;
    }


    /** FG
     * Make region inaccessible by free TFs
     * @param left      inclusive
     * @param right     exclusive
     * @param speciesID species ID (TF type)
     */
    private void closeRegionInAffinityLandscape(int left, int right, int speciesID) {
        int start = Math.max(0, left - TFsize[speciesID] + 1);
        int end = Math.min(right, strand.length);
        for (int j = start; j < end; j++) {
            if (effectiveTFavailability[speciesID][j]) {
                this.effectiveTFsectorsAvailabilitySum[speciesID][this.sectorID[j]]--;
                this.effectiveTFavailabilitySum[speciesID]--;
                effectiveTFavailability[speciesID][j] = false;
            }
        }
    }

    /** FG
     * recomputes the TF affinity landscape when a molecule binds from the DNA
     */
    private void recomputeTFAffinityLandscapeOnBinding(int position, int size) {
        recomputeTFAffinityLandscapeOnRepression(position, position + size - 1);
    }

    /** FG
     * recomputes the TF affinity landscape when a molecule represses the DNA region
     */
    private void recomputeTFAffinityLandscapeOnRepression(int boundaryLeft, int boundaryRight) {
        //recompute affinity landscape for each TF species
        for (int speciesID = 0; speciesID < TFsize.length; speciesID++) {
            closeRegionInAffinityLandscape(boundaryLeft, boundaryRight + 1, speciesID);
        }
    }

    /** FG
     * recomputes the TF affinity landscape when the DNA region is derepressed
     */
    private void recomputeTFAffinityLandscapeOnDerepression(Cell n, int boundaryLeft, int boundaryRight, int proteinID,
                                                            double time)
    {
        int start, end, boundMoleculeID;

        //recompute affinity landscape for each TF species
        for (int speciesID = 0; speciesID < TFsize.length; speciesID++) {
            start = Math.max(0, boundaryLeft - TFsize[speciesID] + 1);
            end = Math.min(this.strand.length - TFsize[speciesID] + 1, boundaryRight + 1);
            for (int bpIdx = start; bpIdx < end; bpIdx++) {
                if (!effectiveTFavailability[speciesID][bpIdx] && closed[bpIdx] == Constants.BP_IS_OPEN) {
                    this.effectiveTFsectorsAvailabilitySum[speciesID][this.sectorID[bpIdx]]++;
                    this.effectiveTFavailabilitySum[speciesID]++;
                    effectiveTFavailability[speciesID][bpIdx] = true;
                }
            }
            recomputeTFAffinityLandscapeForClosedRegions(start, end + TFsize[speciesID] - 1, speciesID);
        }

        // find all repressors in the repressed region and close chromatin in their repression regions
        start = updateLeftBoundary(boundaryLeft - n.maxRepressionLeftSize - n.maxTFSize);
        end = updateRightBoundary(boundaryRight + n.maxRepressionRightSize + n.maxTFSize);
        for (int bpIdx = start; bpIdx <= end; bpIdx++) {
            boundMoleculeID = getBoundMolecule(bpIdx);
            if (boundMoleculeID != Constants.NONE) {
                int size = n.dbp[boundMoleculeID].size;
                for (int speciesID = 0; speciesID < TFsize.length; speciesID++) {
                    closeRegionInAffinityLandscape(bpIdx, bpIdx + size, speciesID);
                }
                if (n.dbp[boundMoleculeID].isRepressingDNA() && boundMoleculeID != proteinID) {
                    int speciesID = n.dbp[boundMoleculeID].speciesID;
                    int boundaryLeft_ = bpIdx - n.TFspecies[speciesID].repressionLeftSize;
                    int boundaryRight_ = bpIdx + size + n.TFspecies[speciesID].repressionRightSize;
                    boundaryLeft_ = updateLeftBoundary(boundaryLeft_);
                    boundaryRight_ = updateRightBoundary(boundaryRight_);
                    for (int speciesID_ = 0; speciesID_ < TFsize.length; speciesID_++) {
                        closeRegionInAffinityLandscape(boundaryLeft_, boundaryRight_, speciesID_);
                    }
                    this.repressDNA(boundaryLeft_, boundaryRight_ - 1, n);
                    if (n.isInDebugMode()) {
                        n.printDebugInfo(time + ": TF " + boundMoleculeID + " of type " + n.TFspecies[speciesID].name +
                                " continued repression from position " + boundaryLeft_ + " to " + (boundaryRight_ - 1));
                    }
                }
                bpIdx += size - 1;
            }
        }
    }

    /**
     * recomputes the TF affinity landscape when a molecule unbinds from the DNA
     */
    private void recomputeTFAffinityLandscapeOnUnbinding(int position, int size) {
        int start, end, pos;
        int[] buffer;

        //recompute affinity landscape for each TF species
        for (int i = 0; i < TFsize.length; i++) {

            start = Math.max(0, position - TFsize[i] + 1);
            if (!this.effectiveTFavailability[i][Math.max(0, start - 1)]) {
                buffer = this.findLastBoundMolecule(start, position - 1);
                if (buffer[0] != Constants.NONE) {
                    start = buffer[1] + 2;
                }
                pos = this.findLastClosedBP(start, position - 1);
                if (pos != Constants.NONE) {
                    start = Math.max(pos + 1, start);
                }
            }

            end = Math.min(position + size, this.strand.length - 1);
            if (!this.effectiveTFavailability[i][end]) {
                buffer = this.findFirstBoundMolecule(end, Math.min(end + TFsize[i] - 1, strand.length - 1));
                if (buffer[0] != Constants.NONE) {
                    end = buffer[1] - TFsize[i];
                }
                pos = this.findFirstClosedBP(position + size, Math.min(end + TFsize[i] - 1, strand.length - 1));
                if (pos != Constants.NONE) {
                    end = Math.min(pos - TFsize[i] + 1, end);
                }
            }

            end = Math.min(strand.length - TFsize[i] + 1, end);
            for (int j = start; j < end; j++) {
                if (!effectiveTFavailability[i][j] && closed[j] == Constants.BP_IS_OPEN) {
                    effectiveTFavailability[i][j] = true;
                    this.effectiveTFsectorsAvailabilitySum[i][this.sectorID[j]]++;
                    this.effectiveTFavailabilitySum[i]++;
                }
            }
            recomputeTFAffinityLandscapeForRepressedRegions(start, end + TFsize[i] - 1, i);
        }
    }


    /** FG
     * Mark open bps between given boundaries as repressed
     * @param boundaryLeft  inclusive
     * @param boundaryRight inclusive
     * @param n             cell
     */
    public void repressDNA(int boundaryLeft, int boundaryRight, Cell n) {
        for (int pos = boundaryLeft; pos <= boundaryRight; pos++) {
            if (this.closed[pos] == Constants.BP_IS_OPEN) {
                this.closed[pos] = Constants.BP_IS_REPRESSED;
                if (n.ip.OUTPUT_REPRESSED_LENGTHS.value) {
                    this.currentRepressedLength++;
                    modifyRepressionScore(n, pos, true);
                }
            }
        }
    }

    public void modifyRepressionScore(Cell n, int pos, boolean increase) {
        for (int speciesID = 0; speciesID < n.TFspecies.length; speciesID++) {
            for (int dir = 0; dir < this.TFdirections; dir++) {
                double moveRate = n.TFspecies[speciesID].calcMoveRate(n.dna.TFavgMoveRate[speciesID][pos][dir], false);
                double scoreIncrement = increase ? (1.0 / moveRate) : (-1.0 / moveRate);
                if (n.TFspecies[speciesID].isRepressor()) {
                    this.repressedRepScore += scoreIncrement;
                } else {
                    this.repressedActScore += scoreIncrement;
                }
            }
        }
    }

    /** FG
     * Mark open bps between given boundaries as repressed and recompute affinity landscape
     * @param boundaryLeft  inclusive
     * @param boundaryRight inclusive
     * @param n             cell
     */
    public void repress(int boundaryLeft, int boundaryRight, Cell n) {
        repressDNA(boundaryLeft, boundaryRight, n);
        this.recomputeTFAffinityLandscapeOnRepression(boundaryLeft, boundaryRight);
    }


    /** FG
     * Update bp index so that it does not go outside the strand
     * @param left bp strand index
     * @return max of 0 (leftmost bp strand index) and left
     */
    public int updateLeftBoundary(int left) {
        return Math.max(0, left);
    }

    /** FG
     * Update bp index so that it does not go outside the strand
     * @param right bp strand index
     * @return min of strand length (rightmost bp strand index) and right
     */
    public int updateRightBoundary(int right) {
        return Math.min(this.strand.length - 1, right);
    }

    /** FG
     * Mark the repressed bps between given boundaries as open
     * @param boundaryLeft  inclusive
     * @param boundaryRight inclusive
     * @param n             cell
     */
    public void derepressDNA(int boundaryLeft, int boundaryRight, Cell n) {
        for (int pos = boundaryLeft; pos <= boundaryRight; pos++) {
            if (this.closed[pos] == Constants.BP_IS_REPRESSED) {
                this.closed[pos] = Constants.BP_IS_OPEN;
                if (n.ip.OUTPUT_REPRESSED_LENGTHS.value) {
                    this.currentRepressedLength--;
                    modifyRepressionScore(n, pos, false);
                }
            }
        }
        assert (this.currentRepressedLength >= 0
                && this.repressedRepScore >= -Constants.DOUBLE_ZERO * Constants.UPDATE_REPRESSION_SCORES_EVERY
                && this.repressedActScore >= -Constants.DOUBLE_ZERO * Constants.UPDATE_REPRESSION_SCORES_EVERY);
    }

    /** FG
     * Mark the repressed bps between given boundaries as open and recompute affinity landscape
     * @param boundaryLeft  inclusive
     * @param boundaryRight inclusive
     * @param n             cell
     */
    public void derepress(Cell n, int boundaryLeft, int boundaryRight, int proteinID, double time) {
        derepressDNA(boundaryLeft, boundaryRight, n);
        this.recomputeTFAffinityLandscapeOnDerepression(n, boundaryLeft, boundaryRight, proteinID, time);
    }

    /**
     * recomputes the TF affinity landscape when a TF molecule slides left on the DNA
     */
    private void recomputeTFAffinityLandscapeOnTFSlideLeft(int position, int moleculeSize, int stepSize) {
        int start, end, pos;
        int[] buffer;

        //recompute affinity landscape for each TF species
        for (int i = 0; i < TFsize.length; i++) {
            //Remove affinities from the left side
            start = Math.max(0, position - stepSize - TFsize[i] + 1);
            end = position - TFsize[i] + 1;
            for (int j = start; j < end; j++) {
                if (effectiveTFavailability[i][j]) {
                    this.effectiveTFsectorsAvailabilitySum[i][this.sectorID[j]]--;
                    this.effectiveTFavailabilitySum[i]--;
                    effectiveTFavailability[i][j] = false;
                }
            }
            //Reallocate affinities to the right side
            start = Math.max(0, position + moleculeSize - stepSize);
            end = position + moleculeSize;
            if (end < strand.length && !this.effectiveTFavailability[i][end]) {
                buffer = this.findFirstBoundMolecule(end, Math.min(end + TFsize[i], strand.length));
                if (buffer[0] != Constants.NONE) {
                    end = buffer[1] - TFsize[i];
                }
                pos = this.findFirstClosedBP(position + moleculeSize, Math.min(end + TFsize[i], strand.length));
                if (pos != Constants.NONE) {
                    end = Math.min(pos - TFsize[i] + 1, end);
                }
            }
            int maxPos = strand.length - TFsize[i] + 1;
            for (int j = start; j < end; j++) {
                if (j < maxPos && !effectiveTFavailability[i][j] && closed[j] == Constants.BP_IS_OPEN) {
                    effectiveTFavailability[i][j] = true;
                    this.effectiveTFsectorsAvailabilitySum[i][this.sectorID[j]]++;
                    this.effectiveTFavailabilitySum[i]++;
                }
            }
        }

    }


    /**
     * recomputes the TF affinity landscape when a TF molecule slides right on the DNA
     */
    private void recomputeTFAffinityLandscapeOnTFSlideRight(int position, int moleculeSize, int stepSize) {
        int start, end, pos;
        int[] buffer;

        //recompute affinity landscape for each TF species
        for (int i = 0; i < TFsize.length; i++) {
            //Remove affinities from the right side
            start = Math.min(position + moleculeSize, strand.length);
            end = Math.min(start + stepSize, strand.length);
            for (int j = start; j < end; j++) {
                if (effectiveTFavailability[i][j]) {
                    this.effectiveTFsectorsAvailabilitySum[i][this.sectorID[j]]--;
                    this.effectiveTFavailabilitySum[i]--;
                    effectiveTFavailability[i][j] = false;
                }
            }
            //Reallocate affinities to the left side
            start = Math.max(0, position - TFsize[i] + 1);
            end = position + stepSize - TFsize[i] + 1;
            if (!this.effectiveTFavailability[i][Math.max(0, start - 1)]) {
                buffer = this.findLastBoundMolecule(start, position - 1);
                if (buffer[0] != Constants.NONE) {
                    start = buffer[1] + 2; // +2 is because of the peculiar findLastBoundMolecule function
                    // buffer[1]+2 gives the first free position
                }
                pos = this.findLastClosedBP(start, position - 1);
                if (pos != Constants.NONE) {
                    start = Math.max(pos + 1, start);
                }
            }
            int maxPos = strand.length - TFsize[i] + 1;
            for (int j = start; j < end; j++) {
                if (j < maxPos && !effectiveTFavailability[i][j] && closed[j] == Constants.BP_IS_OPEN) {
                    effectiveTFavailability[i][j] = true;
                    this.effectiveTFsectorsAvailabilitySum[i][this.sectorID[j]]++;
                    this.effectiveTFavailabilitySum[i]++;
                }
            }
        }
    }


    /**
     * searches the last bound protein on the DNA within a specified interval
     *
     * @param start inclusive
     * @param end   inclusive
     * @return an array which contains on first position the found molecule and on the second one the position
     */
    public int[] findLastBoundMolecule(int start, int end) {
        start = Math.max(0, start);
        end = Math.min(end, strand.length - 1);

        int[] result = new int[2];
        result[0] = Constants.NONE;
        result[1] = end;

        while (result[1] >= start && result[0] == Constants.NONE) {
            result[0] = this.occupied[result[1]];
            result[1]--;
        }

        return result;
    }

    /**
     * searches the last closed protein on the DNA within a specified interval
     *
     * @param start inclusive
     * @param end   inclusive
     * @return position of the last closed bp
     */
    public int findLastClosedBP(int start, int end) {
        start = Math.max(0, start);
        end = Math.min(end, strand.length - 1);
        int pos;
        for (pos = end; pos >= start; pos--) {
            if (closed[pos] != Constants.BP_IS_OPEN) {
                return pos;
            }
        }
        return Constants.NONE;
    }

    /**
     * searches the first closed protein on the DNA within a specified interval
     *
     * @param start inclusive
     * @param end   inclusive
     * @return position of the first closed bp
     */
    public int findFirstClosedBP(int start, int end) {
        start = Math.max(0, start);
        end = Math.min(end, strand.length - 1);
        int pos;
        for (pos = start; pos <= end; pos++) {
            if (closed[pos] != Constants.BP_IS_OPEN) {
                return pos;
            }
        }
        return Constants.NONE;
    }

    /**
     * searches the first bound protein on the DNA within a specified interval
     *
     * @param start inclusive
     * @param end   inclusive
     * @return an array which contains on first position the found molecule and on the second one the position
     */
    public int[] findFirstBoundMolecule(int start, int end) {
        start = Math.max(0, start);
        end = Math.min(end, strand.length - 1);

        int[] result = new int[2];
        result[0] = Constants.NONE;
        result[1] = start;
        while (result[1] <= end && result[0] == Constants.NONE) {
            result[0] = this.occupied[result[1]];
            result[1]++;
        }

        return result;
    }

    /**
     * slides to right a protein
     *
     * @param proteinID      the bound protein
     * @param position       the position where it is bound
     * @param proteinSize    the size of the protein
     * @param stepSize       the size of the step
     * @param checkOccupancy whether to check occupancy or not
     * @return -1 attempt to bind outside the DNA, proteinID if bound on the DNA, or the ID of the first protein
     * blocking the binding on the DNA
     */
    public int slideRight(int proteinID, int position, int proteinSize, int stepSize, boolean checkOccupancy) {
        int canSlide = proteinID;

        int newPosition = position + stepSize;

        //attempt to bind outside the strand interval, if reflexive don't allow it
        if (this.isReflexive && (newPosition < 0 || newPosition >= this.strand.length - proteinSize)) {
            canSlide = Constants.NONE;
        }

        //the DNA is already occupied
        if (canSlide == proteinID && checkOccupancy) {
            canSlide = getBoundProtein(position + proteinSize, stepSize);
            int bpStatus = checkAvailability(position + proteinSize, stepSize);
            if (canSlide == Constants.NONE && bpStatus == Constants.BP_IS_OPEN) {
                canSlide = proteinID;
            }
        }

        //bind the protein
        if (canSlide == proteinID) {
            occupyDNA(proteinID, position + proteinSize, stepSize);
            freeDNA(position, stepSize);
            this.recomputeTFAffinityLandscapeOnTFSlideRight(position, proteinSize, stepSize);
        }

        return canSlide;
    }


    /**
     * slides to left a protein
     *
     * @param proteinID      the bound protein
     * @param position       the position where it is bound
     * @param proteinSize    the size of the protein
     * @param stepSize       the size of the step
     * @param checkOccupancy whether to check occupancy or not
     * @return -1 if attempted to slide outside the DNA or to the closed DNA regions, proteinID if sliding is
     * successful, or the ID of the first blocking protein
     */
    public int slideLeft(int proteinID, int position, int proteinSize, int stepSize, boolean checkOccupancy) {
        int canSlide = proteinID;

        int newPosition = position - stepSize;

        //attempt to bind outside the strand interval, if reflexive don't allow it
        if (this.isReflexive && newPosition < 0 || newPosition >= this.strand.length - proteinSize) {
            canSlide = Constants.NONE;
        }

        //the DNA is already occupied
        if (canSlide == proteinID && checkOccupancy) {
            canSlide = getBoundProtein(position - stepSize, stepSize);
            int bpStatus = checkAvailability(position - stepSize, stepSize);
            if (canSlide == Constants.NONE && bpStatus == Constants.BP_IS_OPEN) {
                canSlide = proteinID;
            }
        }

        //bind the protein
        if (canSlide == proteinID) {
            occupyDNA(proteinID, position - stepSize, stepSize);
            freeDNA(position + proteinSize - stepSize, stepSize);
            this.recomputeTFAffinityLandscapeOnTFSlideLeft(position, proteinSize, stepSize);
        }

        return canSlide;
    }

    /**
     * occupies a DNA region with a protein
     *
     * @param proteinID the ID of the protein to be bound
     * @param position  the position on the dna to bind to
     * @param size      the size of the protein
     */
    public void occupyDNA(int proteinID, int position, int size) {
        for (int i = 0; i < size; i++) {
            this.occupied[i + position] = proteinID;
        }
    }

    /**
     * frees a DNA region
     *
     * @param position the position on the dna to bind to
     * @param size     the size of the protein
     */
    public void freeDNA(int position, int size) {
        int end = Math.min(position + size, strand.length);
        int start = Math.max(0, position);

        for (int i = start; i < end; i++) {
            this.occupied[i] = Constants.NONE;
        }
    }

    /**
     * FG
     *
     * @param position starting position
     * @param size     length of the sequence to open
     */
    public void openDNA(int position, int size) {
        int end = Math.min(position + size, strand.length);
        int start = Math.max(0, position);

        for (int i = start; i < end; i++) {
            this.closed[i] = Constants.BP_IS_OPEN;
        }
    }

    /**
     * returns the left neighbour position of a TF
     *
     * @param position the position
     */
    public int getLeftNeighbour(int position) {
        return position > 0 ? this.occupied[position - 1] : Constants.NONE;
    }


    /**
     * returns the right neighbour position of a TF
     *
     * @param position the position
     */
    public int getRightNeighbour(int position) {
        return position < strand.length - 1 ? this.occupied[position + 1] : Constants.NONE;
    }

    /**
     * checks whether a protein can bind or not to a position, it returns the id of the first encountered protein
     *
     * @param position the start position on the DNA
     * @param size     the size in bp to check
     * @return true if the specified cluster is free or false otherwise
     */
    public int getBoundProtein(int position, int size) {
        int proteinID = Constants.NONE;

        for (int i = 0; i < size && proteinID == Constants.NONE; i++) {
            proteinID = this.occupied[i + position];
        }
        return proteinID;
    }

    /**
     * FG
     *
     * @param position left boundary of the sequence
     * @param size     length of the sequence
     * @return status (all bps in the sequence are open or not)
     */
    public int checkAvailability(int position, int size) {
        int bpSate = Constants.BP_IS_OPEN;
        for (int i = 0; i < size && bpSate == Constants.BP_IS_OPEN; i++) {
            bpSate = this.closed[position + i];
        }
        return bpSate;
    }

    /**
     * prints a string with the DNA TF affinities
     */
    public void printDNAoccupancy(String path, String filename, int start, int end, boolean fullOccupancy,
                                  int wigStepSize, double wigThreshold)
    {
        double[][][] bufferTFoccupancy = new double[effectiveTFOccupancy.length][strand.length][this.TFdirections];

        double[][] cutoff = new double[effectiveTFOccupancy.length][this.TFdirections];
        double[][] avg = new double[effectiveTFOccupancy.length][this.TFdirections];

        int max;
        //normalise the TF occupancy
        for (int dir = 0; dir < TFdirections; dir++) {
            for (int i = 0; i < effectiveTFOccupancy.length; i++) {
                //init occupancy
                for (int j = 0; j < strand.length; j++) {
                    bufferTFoccupancy[i][j][dir] = 0;
                }

                // add occupancy on the entire length of the TFs
                for (int j = 0; j < strand.length; j++) {
                    max = Math.min(strand.length, j + 1);
                    if (fullOccupancy) {
                        max = Math.min(strand.length, j + TFsize[i]);
                    }
                    for (int k = j; k < max; k++) {
                        bufferTFoccupancy[i][k][dir] += effectiveTFOccupancy[i][j][dir];

                    }
                }

                //normalise
                cutoff[i][dir] = 0;
                for (int j = 0; j < strand.length; j++) {
                    if (cutoff[i][dir] < bufferTFoccupancy[i][j][dir]) {
                        cutoff[i][dir] = bufferTFoccupancy[i][j][dir];
                    }
                }
                cutoff[i][dir] = cutoff[i][dir] * wigThreshold;


                //computes the threshold of the wig if this is set to autoselect
                avg[i][dir] = 0;
                for (int j = 0; j < strand.length; j++) {
                    avg[i][dir] += bufferTFoccupancy[i][j][dir];
                }
                avg[i][dir] /= strand.length;
            }
        }

        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path + filename));

            //check start end output
            if (start < 0 || start >= strand.length) {
                start = 0;
            }

            end = Math.min(end, strand.length);
            if (end < start) {
                end = strand.length;
            }

            //header
            int collisionCount;
            double occupancy;
            int steps;
            StringBuilder strBuf = new StringBuilder();
            strBuf.append("fixedStep  chrom=");
            strBuf.append(this.region.chromosome);
            strBuf.append("  start=");
            strBuf.append(start);
            strBuf.append("  step=");
            strBuf.append(wigStepSize);
            strBuf.append("  span=");
            strBuf.append(wigStepSize);
            strBuf.append("\n");
            strBuf.append("\"position\", \"collisionsCount\"");
            if (this.TFdirections == 1) {
                for (int j = 0; j < TFsize.length; j++) {
                    strBuf.append(", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("\"");
                }
            } else {
                for (int j = 0; j < this.TFsize.length; j++) {
                    strBuf.append(", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("5'3'\", \"");
                    strBuf.append(TFposName.get(j));
                    strBuf.append("3'5'\"");
                }
            }
            out.write(strBuf.toString());
            out.newLine();

            //info
            for (int i = start; i < end; i = i + wigStepSize) {
                //collision count
                collisionCount = 0;
                steps = 0;
                for (int k = i; k < Math.min(i + wigStepSize, end); k++) {
                    collisionCount += this.collisionsCount[k];
                    steps++;
                }
                strBuf.delete(0, strBuf.length());
                strBuf.append((i + this.subsequence.start));
                strBuf.append(", ");
                strBuf.append(((double) collisionCount / (steps)));

                //TF occupancy
                for (int j = 0; j < bufferTFoccupancy.length; j++) {
                    for (int dir = 0; dir < this.TFdirections; dir++) {
                        occupancy = 0;
                        steps = 0;
                        for (int k = i; k < Math.min(i + wigStepSize, end); k++) {
                            if ((wigThreshold >= 0 && bufferTFoccupancy[j][k][dir] > cutoff[j][dir])
                                    || (wigThreshold <= 0 && bufferTFoccupancy[j][k][dir] > avg[j][dir])) {
                                occupancy += bufferTFoccupancy[j][k][dir];
                            }
                            steps++;
                        }
                        strBuf.append(", ");
                        strBuf.append((occupancy / steps));
                    }
                }
                out.write(strBuf.toString());
                out.newLine();
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * returns a molecule which is bound at starting position specified as a parameter
     *
     * @param position the starting position
     */
    public int getBoundMolecule(int position) {
        return this.occupied[position] != Constants.NONE && (position == 0 || occupied[position - 1] != occupied[position]) ? occupied[position] : Constants.NONE;
    }
}
