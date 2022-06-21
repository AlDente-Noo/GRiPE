package objects;

import environment.Cell;
import utils.CellUtils;
import utils.Constants;
import utils.Utils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * class that contains the description of TF species
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
public class TFSpecies implements Serializable {

    private static final long serialVersionUID = -8664545127198799155L;
    private final HashMap<Integer, Integer> directCooperativeSpecies;
    public int id; // the id of the TF species
    public String name; // the id of the TF species
    public byte[] dbd; // the sequence regonise at DNA binding domain
    public PFM pfm;
    public String landscapeFile;
    public int landscapePosCol;
    public int landscapeAffinityColLR;
    public int landscapeAffinityColRL;
    public int landscapeEscapeLines;
    public String seqsFile;
    public double seqsDefaultValue;
    public int seqsEscapeLines;
    public String dbdFile;
    public int copyNumber; // number of molecules of this species
    public double es; // specific binding energy per nucleotide
    public int sizeLeft;
    public int sizeRight;
    public int sizeTotal;
    public double assocRate;
    public DNAregion initialDrop;
    public DNAregion relInitialDrop;
    public double maxAffinity;
    public double timeBoundAvg;
    public double timeBoundVar;
    public boolean isCognate;
    public double unBindingProbability;
    public double slideLeftProbability;
    public double slideRightProbability;
    public double jumpingProbability;
    public double hopSTDdisplacement;
    public double specificWaitingTime;
    public int stepLeftSize;
    public int stepRightSize;
    public int uncorrelatedDisplacementSize;
    public boolean stallsIfBlocked;
    public double collisionUnbindingProbability;
    public double affinityLandscapeRoughness;
    public double repressionRate;
    public double repressionAttenuationFactor;
    public int repressionLeftSize;
    public int repressionRightSize;
    public double jumpNo;
    public double hopNo;
    public double slideLeftNo;
    public double slideRightNo;
    public double preboundProportion;
    public boolean preboundToHighestAffinity;
    public long countTFBindingEvents;
    public long countTFUnbindingEvents;
    public long countTFSlideLeftEvents;
    public long countTFSlideRightEvents;
    public long countTFHoppingEvents;
    public long countTFRepressionEvents;
    public long countTFDerepressionEvents;
    public long countTFforcedJumpsEvents;
    public long countTFHopsOutside;
    public double slidingEventsPerBinding;
    public double slidingLengthPerBinding;
    public double observedSlidingLengthPerBinding;
    public double residenceTimePerBinding;
    public boolean isCooperative;
    public boolean hasDNAbasedCooperativity;
    public boolean hasDirectCooperativity;
    public ArrayList<TFcooperativity> TFcoop;
    public int[][] isCooperativeSite;
    public ArrayList<Integer> slidingLength;
    public ArrayList<Integer> observedSlidingLength;
    public ArrayList<Integer> slidingEvents;
    public boolean isBiasedRandomWalk;
    public boolean isTwoStateRandomWalk;
    public boolean isImmobile;

    //public double moveRateThreshold;

    /**
     * class constructor for the random tfs
     */
    public TFSpecies(int id, byte[] DNAstrand, int pos, int dbdLength, int copyNumber, double es, int sizeLeft,
                     int sizeRight, double assocRate, DNAregion initialDrop, boolean isCognate,
                     double unBindingProbability, double slideLeftProbability, double slideRightProbability,
                     double jumpingProbability, double hopSTDdisplacement, double specificWaitingTime,
                     int stepLeftSize, int stepRightSize, int uncorrelatedDisplacementSize, boolean stallsIfBlocked,
                     double collisionUnbindingProbability, double affinityLandscapeRoughness,
                     double preboundProportion, boolean preboundToHighestAffinity, boolean isImmobile,
                     boolean isBiasedRandomWalk, boolean isTwoStateRandomWalk, double repressionRate,
                     double repressionAttenuationFactor, int repressionLenLeft, int repressionLenRight) {
        this.id = id;
        name = "TF" + id;
        this.es = es;
        this.copyNumber = copyNumber;
        pfm = null;
        this.landscapeFile = "";
        this.seqsFile = "";
        dbd = new byte[dbdLength];
        maxAffinity = 0;
        for (int i = pos; i < pos + dbdLength && i < DNAstrand.length; i++) {
            dbd[i - pos] = DNAstrand[i];
            if (dbd[i - pos] != CellUtils.bpANYID) {
                maxAffinity += es;
            }
        }
        this.sizeLeft = sizeLeft;
        this.sizeRight = sizeRight;
        this.sizeTotal = sizeLeft + sizeRight + dbdLength;
        this.assocRate = assocRate;
        this.timeBoundAvg = 0;
        this.timeBoundVar = 0;

        this.initialDrop = initialDrop;
        this.isCognate = isCognate;

        this.unBindingProbability = unBindingProbability;
        this.slideLeftProbability = slideLeftProbability;
        this.slideRightProbability = slideRightProbability;
        this.jumpingProbability = jumpingProbability;
        this.repressionRate = repressionRate;
        this.repressionAttenuationFactor = repressionAttenuationFactor;
        this.repressionLeftSize = repressionLenLeft;
        this.repressionRightSize = repressionLenRight;
        this.hopSTDdisplacement = hopSTDdisplacement;
        this.specificWaitingTime = specificWaitingTime;
        this.stepLeftSize = stepLeftSize;
        this.stepRightSize = stepRightSize;
        this.uncorrelatedDisplacementSize = uncorrelatedDisplacementSize;
        this.stallsIfBlocked = stallsIfBlocked;
        this.collisionUnbindingProbability = collisionUnbindingProbability;
        this.affinityLandscapeRoughness = affinityLandscapeRoughness;
        this.preboundProportion = preboundProportion;
        this.preboundToHighestAffinity = preboundToHighestAffinity;
        jumpNo = this.unBindingProbability * this.jumpingProbability;
        hopNo = this.unBindingProbability;
        slideLeftNo = hopNo + this.slideLeftProbability;
        slideRightNo = slideLeftNo + this.slideRightProbability;

        this.countTFBindingEvents = 0;
        this.countTFUnbindingEvents = 0;
        this.countTFSlideLeftEvents = 0;
        this.countTFSlideRightEvents = 0;
        this.countTFHoppingEvents = 0;
        this.countTFforcedJumpsEvents = 0;
        this.countTFHopsOutside = 0;
        this.countTFRepressionEvents = 0;

        this.isCooperative = false;
        hasDNAbasedCooperativity = false;
        hasDirectCooperativity = false;
        TFcoop = new ArrayList<TFcooperativity>();
        directCooperativeSpecies = new HashMap<Integer, Integer>();
        slidingLength = new ArrayList<Integer>();
        slidingEvents = new ArrayList<Integer>();
        observedSlidingLength = new ArrayList<Integer>();

        this.isImmobile = isImmobile;

        this.isBiasedRandomWalk = isBiasedRandomWalk;
        this.isTwoStateRandomWalk = isTwoStateRandomWalk;
        this.dbdFile = "";
    }


    /**
     * class constructor
     */
    public TFSpecies(DNAregion dnaRegion, int id, String name, String dbd, int copyNumber, double es, int sizeLeft,
                     int sizeRight, double assocRate, DNAregion initialDrop, boolean isCognate,
                     double unBindingProbability, double slideLeftProbability, double slideRightProbability,
                     double jumpingProbability, double hopSTDdisplacement, double specificWaitingTime,
                     int stepLeftSize, int stepRightSize, int uncorrelatedDisplacementSize, boolean stallsIfBlocked,
                     double collisionUnbindingProbability, double affinityLandscapeRoughness,
                     double preboundProportion, boolean preboundToHighestAffinity, boolean isImmobile,
                     boolean isBiasedRandomWalk, boolean isTwoStateRandomWalk, double repressionRate,
                     double repressionAttenuationFactor, int repressionLenLeft, int repressionLenRight, Cell n) {
        this.id = id;
        this.name = name;
        this.isCognate = false;

        this.dbd = new byte[0];
        int sizeInBP = 0;
        pfm = null;
        this.landscapeFile = "";
        this.seqsFile = "";
        this.dbdFile = "";

        this.sizeLeft = sizeLeft;
        this.sizeRight = sizeRight;

        sizeInBP = parseDBD(n, dbd);

        this.copyNumber = copyNumber;
        this.es = es;

        this.sizeTotal = sizeLeft + sizeRight + sizeInBP;
        this.assocRate = assocRate;
        this.timeBoundAvg = 0;
        this.timeBoundVar = 0;

        //initial drop
        relInitialDrop = new DNAregion("", initialDrop.start - dnaRegion.start, initialDrop.end - dnaRegion.start);
        if (relInitialDrop.start > 0 && relInitialDrop.end - 1 < dnaRegion.size()) {
            this.initialDrop = initialDrop;
        } else {
            this.initialDrop = new DNAregion("", 0, dnaRegion.size());
            this.relInitialDrop = new DNAregion("", 0, dnaRegion.size());
        }

        this.unBindingProbability = unBindingProbability;
        this.slideLeftProbability = slideLeftProbability;
        this.slideRightProbability = slideRightProbability;
        this.jumpingProbability = jumpingProbability;
        this.repressionRate = repressionRate;
        this.repressionAttenuationFactor = repressionAttenuationFactor;
        this.repressionLeftSize = repressionLenLeft;
        this.repressionRightSize = repressionLenRight;
        this.hopSTDdisplacement = hopSTDdisplacement;
        this.specificWaitingTime = specificWaitingTime;
        this.stepLeftSize = stepLeftSize;
        this.stepRightSize = stepRightSize;
        this.uncorrelatedDisplacementSize = uncorrelatedDisplacementSize;
        this.stallsIfBlocked = stallsIfBlocked;
        this.collisionUnbindingProbability = collisionUnbindingProbability;
        this.affinityLandscapeRoughness = affinityLandscapeRoughness;
        this.preboundProportion = preboundProportion;
        this.preboundToHighestAffinity = preboundToHighestAffinity;
        jumpNo = this.unBindingProbability * this.jumpingProbability;
        hopNo = this.unBindingProbability;
        slideLeftNo = hopNo + this.slideLeftProbability;
        slideRightNo = slideLeftNo + this.slideRightProbability;

        this.countTFBindingEvents = 0;
        this.countTFUnbindingEvents = 0;
        this.countTFSlideLeftEvents = 0;
        this.countTFSlideRightEvents = 0;
        this.countTFHoppingEvents = 0;
        this.countTFforcedJumpsEvents = 0;
        this.countTFHopsOutside = 0;
        this.countTFRepressionEvents = 0;

        this.isCooperative = true;
        TFcoop = new ArrayList<TFcooperativity>();
        directCooperativeSpecies = new HashMap<Integer, Integer>();

        this.isImmobile = isImmobile;

        slidingLength = new ArrayList<Integer>();
        observedSlidingLength = new ArrayList<Integer>();
        slidingEvents = new ArrayList<Integer>();
        this.isBiasedRandomWalk = isBiasedRandomWalk;
        this.isTwoStateRandomWalk = isTwoStateRandomWalk;
    }

    public boolean isRepressor() {
        return this.repressionRate > 0.0;
    }

    /**
     * parse a String into DBD
     *
     * @param n   the cell
     * @param dbd the DNA binding domain text
     */
    public int parseDBD(Cell n, String dbd) {
        int sizeInBP = 0;
        if (dbd.startsWith(Constants.DBD_TYPE_SEQ)) {
            ArrayList<Byte> bufferDBD;
            dbd = dbd.replaceAll(Constants.DBD_TYPE_SEQ, "");
            bufferDBD = CellUtils.getSequenceIDs(dbd);

            this.dbd = new byte[bufferDBD.size()];
            this.maxAffinity = 0;
            for (int i = 0; i < bufferDBD.size(); i++) {
                this.dbd[i] = bufferDBD.get(i);
                if (this.dbd[i] != CellUtils.bpANYID) {
                    maxAffinity += es;
                }
            }
            sizeInBP = this.dbd.length;
            if (this.dbd.length > 0) {
                this.isCognate = true;
            }
        } else if (dbd.startsWith(Constants.DBD_TYPE_PFM) || dbd.startsWith(Constants.DBD_TYPE_PWM)) {
            pfm = new PFM(dbd, n);
            if (this.pfm.isCorrect) {
                this.isCognate = true;
            }
            sizeInBP = this.pfm.motifSize;
        } else if (dbd.startsWith(Constants.DBD_TYPE_LANDSCAPE)) {
            String buffer = dbd.replaceAll(Constants.DBD_TYPE_LANDSCAPE, "");
            this.dbdFile = dbd;
            if (buffer != null && !buffer.isEmpty() && buffer.contains(Constants.DBD_TYPE_SEPARATOR)) {

                String[] bufferStr = buffer.split(Constants.DBD_TYPE_SEPARATOR);

                if (bufferStr.length >= 6) {
                    //is cognate
                    this.isCognate = Utils.parseBoolean(bufferStr[1], true);

                    //DBD size in bp
                    sizeInBP = Utils.parseInteger(bufferStr[0], Constants.NONE);
                    if (sizeInBP < 0 || (sizeInBP == 0 && sizeLeft == 0 && sizeRight == 0)) {
                        n.stopSimulation("Could not load the affinity landscape " + dbd + ": no motif size");
                    }

                    if (bufferStr.length == 7) {
                        this.landscapePosCol = Utils.parseInteger(bufferStr[2], Constants.NONE);
                        this.landscapeAffinityColLR = Utils.parseInteger(bufferStr[3], Constants.NONE);
                        this.landscapeAffinityColRL = Utils.parseInteger(bufferStr[4], Constants.NONE);
                        this.landscapeFile = bufferStr[5];
                        this.landscapeEscapeLines = Utils.parseInteger(bufferStr[6], 0);
                        this.landscapePosCol--;
                        if (this.landscapePosCol < 0) {
                            n.stopSimulation("Could not load the affinity landscape " + dbd + ": no column for DNA " +
                                    "position");

                        }
                    } else {
                        this.landscapePosCol = Constants.NONE;
                        this.landscapeAffinityColLR = Utils.parseInteger(bufferStr[2], Constants.NONE);
                        this.landscapeAffinityColRL = Utils.parseInteger(bufferStr[3], Constants.NONE);
                        this.landscapeFile = bufferStr[4];
                        this.landscapeEscapeLines = Utils.parseInteger(bufferStr[5], 0);
                    }
                    this.landscapeAffinityColLR--;
                    this.landscapeAffinityColRL--;

                    if (this.landscapeAffinityColLR < 0) {
                        n.stopSimulation("Could not load the affinity landscape " + dbd + ": no column for affinity " +
                                "in 3'-5' direction");

                    }

                    if (this.landscapeAffinityColRL < 0) {
                        n.stopSimulation("Could not load the affinity landscape " + dbd + ": no column for affinity " +
                                "in 5'-3' direction");
                    }

                    if (this.landscapeFile == null || this.landscapeFile.isEmpty()) {
                        n.stopSimulation("Could not load the affinity landscape " + dbd + ": no file containing the " +
                                "landscape");

                    } else {
                        File f = new File(this.landscapeFile);
                        if (!f.exists() || f.isDirectory() || !f.canRead() || f.length() == 0) {
                            n.stopSimulation("Could not load the affinity landscape " + dbd + ": cannot find the " +
                                    "landscape file: " + landscapeFile);
                        }
                    }
                }
            } else {
                n.stopSimulation("Could not load the affinity landscape " + dbd + ": could not split the text by " + Constants.DBD_TYPE_SEPARATOR);
            }

            //if(!isLandscapeSpecificationCorrect){
            //	n.stopSimulation("Could not load the affinity landscape: "+dbd);
            //}

        } else if (dbd.startsWith(Constants.DBD_TYPE_SEQS)) {
            String buffer = dbd.replaceAll(Constants.DBD_TYPE_SEQS, "");
            this.dbdFile = dbd;

            if (buffer != null && !buffer.isEmpty() && buffer.contains(Constants.DBD_TYPE_SEPARATOR)) {

                String[] bufferStr = buffer.split(Constants.DBD_TYPE_SEPARATOR);

                if (bufferStr.length == 5) {
                    //is cognate
                    this.isCognate = Utils.parseBoolean(bufferStr[1], true);

                    //DBD size in bp
                    sizeInBP = Utils.parseInteger(bufferStr[0], Constants.NONE);
                    if (sizeInBP < 0 || (sizeInBP == 0 && sizeLeft == 0 && sizeRight == 0)) {
                        n.stopSimulation("Could not load the affinity landscape " + dbd + ": no motif size");
                    }

                    this.seqsDefaultValue = Utils.parseDouble(bufferStr[2], Constants.NONE);
                    this.seqsFile = bufferStr[3];
                    this.seqsEscapeLines = Utils.parseInteger(bufferStr[4], 0);

                    if (this.seqsDefaultValue < 0) {
                        n.stopSimulation("Could not load the sequence file " + dbd + ": no default value");
                    }

                    if (this.seqsFile == null || this.seqsFile.isEmpty()) {
                        n.stopSimulation("Could not load the sequence file " + dbd + ": no file containing the " +
                                "sequences");
                    } else {
                        File f = new File(this.seqsFile);
                        if (!f.exists() || f.isDirectory() || !f.canRead() || f.length() == 0) {
                            n.stopSimulation("Could not load the sequence file " + dbd + ": cannot find the sequence " +
                                    "file: " + seqsFile);
                        }
                    }
                }

            } else {
                n.stopSimulation("Could not load the sequence file" + dbd + ": could not split the text by " + Constants.DBD_TYPE_SEPARATOR);
            }
        }
        return sizeInBP;

    }


    /**
     * @return returns the string of the current TF species
     */
    public String toString(Cell n, boolean reduced) {

        //"\"name\", \"DBD\",
        StringBuilder str = new StringBuilder("\"" + name + "\", \"");
        if (pfm != null && pfm.isCorrect && pfm.motifSize > 0) {
            str.append(pfm.toString(n.dna.bpFreq));
        } else if (dbd.length > 0) {
            str.append("SEQ:");
            for (byte b : dbd) {
                str.append(CellUtils.bps.bps[b]);
            }
        } else {
            str.append(this.dbdFile);
        }

        //\"ES\", \"COPYNUMBER\", \"SIZELEFT\", \"SIZERIGHT\", \"ASSOCRATE\", \"INITIALDROP\",
        // \"UNBINDINGPROBABILITY\", \"SLIDELEFTPROBABILITY\", \"SLIDERIGHTPROBABILITY\", \"JUMPINGPROBABILITY\",
        // \"HOPSTDDISPLACEMENT\", \"SPECIFICWAITINGTIME\", \"STEPLEFTSIZE\", \"STEPRIGHTSIZE\",
        // \"UNCORRELATEDDISPLACEMENTSIZE\", \"STALLSIFBLOCKED\", \"COLLISIONUNBINDPROBABILITY\",
        // \"AFFINITYLANDSCAPEROUGHNESS\", \"PREBOUNDPROPORTION\", \"PREBOUNDTOHIGHESTAFFINITY\"
        str.append("\",").append(es).append(", ").append(copyNumber).append(", ").append(this.sizeLeft).append(", ").append(this.sizeRight).append(", ").append(this.assocRate).append(", \"").append(this.initialDrop).append("\"").append(", ").append(unBindingProbability).append(", ").append(slideLeftProbability).append(", ").append(slideRightProbability).append(", ").append(jumpingProbability).append(", ").append(hopSTDdisplacement).append(", ").append(specificWaitingTime).append(", ").append(stepLeftSize).append(", ").append(stepRightSize).append(", ").append(uncorrelatedDisplacementSize).append(", ").append(stallsIfBlocked).append(", ").append(collisionUnbindingProbability).append(", ").append(affinityLandscapeRoughness).append(", ").append(preboundProportion).append(", ").append(preboundToHighestAffinity).append(", ").append(isImmobile);

        //"ISBIASEDRANDOMWALK\", \"ISTWOSTATERANDOMWALK\",
        str.append(", ").append(this.isBiasedRandomWalk).append(", ").append(this.isTwoStateRandomWalk);

        //"REPRESSIONRATE", "DEREPRESSIONATTENUATIONFACTOR", "REPRLENLEFT", "REPRLENRIGHT"
        str.append(", ").append(this.repressionRate).append(", ").append(this.repressionAttenuationFactor).append(", "
        ).append(this.repressionLeftSize).append(", ").append(this.repressionRightSize);

        if (!reduced) {
            //"eventsBindingTotal\", \"eventsUnbindingTotal\", \"eventsSlideLeftTotal\", \"eventsSlideRightTotal\",
            // \"eventsSlideTotal\", \"eventsHoppingTotal\", \"eventsForcedJumps\", \"eventsHopOutsideDNA\",
            // \"collisionsCount\",
            str.append(", ").append(this.countTFBindingEvents).append(", ").append(this.countTFUnbindingEvents).append(", ").append(this.countTFSlideLeftEvents).append(", ").append(this.countTFSlideRightEvents).append(", ").append(this.countTFSlideLeftEvents + this.countTFSlideRightEvents).append(", ").append(this.countTFHoppingEvents).append(", ").append(this.countTFforcedJumpsEvents).append(", ").append(this.countTFHopsOutside).append(", ").append(n.dna.collisionsCountTotal);

            //     \"sizeTotal\"" , \"isCognate\"
            str.append(", ").append(this.sizeTotal).append(", ").append(this.isCognate);

            //, timeBoundAvg \"residenceTimePerBinding\", \"slidingEventsPerBinding\", \"slidingLengthPerBinding\",
            // \"observedSlidingLengthPerBinding\"";
            if (this.slidingLength != null && !slidingLength.isEmpty()) {
                str.append(", ").append(this.timeBoundAvg).append(", ").append(residenceTimePerBinding).append(", ").append(Utils.computeMean(slidingEvents)).append(", ").append(Utils.computeMean(this.slidingLength)).append(", ").append(Utils.computeMean(this.observedSlidingLength));
            } else {
                str.append(", ").append(this.timeBoundAvg).append(", ").append(residenceTimePerBinding).append(", ").append(slidingEventsPerBinding).append(", ").append(this.slidingLengthPerBinding).append(", ").append(this.observedSlidingLengthPerBinding);
            }
            str.append(", \"").append(this.getCoopString()).append("\"");

            str.append(", ").append(this.countTFRepressionEvents).append(", ").append(this.countTFDerepressionEvents);
        }
        return str.toString();
    }

    /**
     * @return returns the string header of the info of the current TF species
     */
    public String headerToString(boolean reduced) {
        String str = "\"name\", \"DBD\", \"ES\", \"COPYNUMBER\", \"SIZELEFT\", \"SIZERIGHT\", \"ASSOCRATE\", " +
                "\"INITIALDROP\", \"UNBINDINGPROBABILITY\", \"SLIDELEFTPROBABILITY\", \"SLIDERIGHTPROBABILITY\", " +
                "\"JUMPINGPROBABILITY\", \"HOPSTDDISPLACEMENT\", \"SPECIFICWAITINGTIME\", \"STEPLEFTSIZE\", " +
                "\"STEPRIGHTSIZE\", \"UNCORRELATEDDISPLACEMENTSIZE\", \"STALLSIFBLOCKED\", " +
                "\"COLLISIONUNBINDPROBABILITY\", \"AFFINITYLANDSCAPEROUGHNESS\", \"PREBOUNDPROPORTION\", " +
                "\"PREBOUNDTOHIGHESTAFFINITY\", \"TFISIMMOBILE\", \"ISBIASEDRANDOMWALK\", \"ISTWOSTATERANDOMWALK\", " +
                "\"REPRESSIONRATE\", \"DEREPRESSIONATTENUATIONFACTOR\", \"REPRLENLEFT\", \"REPRLENRIGHT\"";
        if (!reduced) {
            str += ", \"eventsBindingTotal\", \"eventsUnbindingTotal\", \"eventsSlideLeftTotal\", " +
                    "\"eventsSlideRightTotal\", \"eventsSlideTotal\", \"eventsHoppingTotal\", \"eventsForcedJumps\", " +
                    "\"eventsHopOutsideDNA\", \"collisionsCount\", \"sizeTotal\", \"isCognate\", \"timeBoundAvg\", " +
                    "\"residenceTimePerBinding\", \"slidingEventsPerBinding\", \"slidingLengthPerBinding\", " +
                    "\"observedSlidingLengthPerBinding\", \"cooperativity\"";
            str += ", \"eventsRepression\", \"eventsDerepression\"";
        }
        return str;
    }


    /**
     * adds TF coopearativity
     */
    public void addCooperativity(TFcooperativity coop, int DNAlength, int TFdirections) {
        this.TFcoop.add(coop);

        this.isCooperative = true;
        if (!this.hasDNAbasedCooperativity && coop.type == Constants.TF_COOPERATIVITY_TYPE_DNA) {
            this.hasDNAbasedCooperativity = true;
        }
        if (!this.hasDirectCooperativity && coop.type == Constants.TF_COOPERATIVITY_TYPE_DIRECT) {
            this.hasDirectCooperativity = true;
        }

        //if there is a cooperative site then mark this on the array
        if (coop.type == Constants.TF_COOPERATIVITY_TYPE_DNA) {
            if (this.isCooperativeSite == null) {
                intiateCooperativeSite(DNAlength, TFdirections);
            }
            setSiteAsCooperative(coop.region0, this.TFcoop.size() - 1);
        }

        //if there is a cooperative site then mark this on the array
        if (coop.type == Constants.TF_COOPERATIVITY_TYPE_DIRECT) {
            directCooperativeSpecies.put(coop.species1ID, this.TFcoop.size() - 1);
        }

    }


    /**
     * initialise the cooperative site
     */
    private void intiateCooperativeSite(int DNAlength, int TFdirections) {
        isCooperativeSite = new int[DNAlength][TFdirections];
        for (int i = 0; i < DNAlength; i++) {
            for (int j = 0; j < TFdirections; j++) {
                this.isCooperativeSite[i][j] = Constants.NONE;
            }
        }
    }


    /**
     * sets a site as cooperative
     */
    private void setSiteAsCooperative(DNAregion region, int coopID) {
        int start = (int) Math.max(0, region.start);
        int end = (int) Math.min(1, region.end);
        end = Math.max(start + 1, end);

        int startDir = 0, endDir = this.isCooperativeSite[0].length;
        if (this.TFcoop.get(coopID).direction0 != Constants.NONE) {
            startDir = this.TFcoop.get(coopID).direction0;
            endDir = this.TFcoop.get(coopID).direction0 + 1;
        }

        for (int i = start; i < end; i++) {
            for (int j = startDir; j < endDir; j++) {
                this.isCooperativeSite[i][j] = coopID;
            }
        }
    }

    /**
     * returns true if this site is marked as being cooperative
     */
    public boolean isCooperativeSite(int position, int dir) {
        return this.isCooperativeSite[position][dir] != Constants.NONE;
    }

    /**
     * returns the species that is affected by binding this molecule to region0
     */
    public TFcooperativity getCooperativity(int position, int direction) {
        return this.TFcoop.get(this.isCooperativeSite[position][direction]);
    }

    /**
     * prints the sliding lengths to a file
     */
    public void printSlidingLengths(String path, String filename) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path + filename));
            out.write("\"no\", \"slidingLength\", \"slidingEvents\" \n");
            StringBuilder strBuf = new StringBuilder();

            for (int i = 0; i < this.slidingLength.size(); i++) {
                strBuf.delete(0, strBuf.length());
                strBuf.append(i);
                strBuf.append(", ");
                strBuf.append(this.slidingLength.get(i));
                strBuf.append(", ");
                strBuf.append(this.slidingEvents.get(i));
                out.write(strBuf.toString());
                out.newLine();
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * prints the sliding lengths to a file
     */
    public void printObservedSlidingLengths(String path, String filename) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(path + filename));
            StringBuilder strBuf = new StringBuilder();

            out.write("\"no\", \"observedSlidingLength\" \n");

            for (int i = 0; i < this.observedSlidingLength.size(); i++) {
                strBuf.delete(0, strBuf.length());
                strBuf.append(i);
                strBuf.append(", ");
                strBuf.append(this.observedSlidingLength.get(i));
                out.write(strBuf.toString());
                out.newLine();
            }

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    /**
     * returns the id of the cooperativity for a right neighbour
     *
     * @param speciesID1
     * @param direction0
     * @param direction1
     * @return
     */
    public TFcooperativity getDirectCooperativityRight(int speciesID1, int direction0, int direction1) {
        TFcooperativity result = null;

        if (directCooperativeSpecies.containsKey(speciesID1)) {
            int coopID = directCooperativeSpecies.get(speciesID1);
            if ((TFcoop.get(coopID).direction0 == Constants.NONE || TFcoop.get(coopID).direction0 == direction0) && (TFcoop.get(coopID).direction1 == Constants.NONE || TFcoop.get(coopID).direction1 == direction1)) {
                result = TFcoop.get(coopID);
            }
        }

        return result;
    }


    /**
     * returns the id of the cooperativity for a left neighbour
     *
     * @param speciesID1
     * @param direction0
     * @param direction1
     * @return
     */
    public TFcooperativity getDirectCooperativityLeft(int speciesID1, int direction0, int direction1) {
        TFcooperativity result = null;

        if (directCooperativeSpecies.containsKey(speciesID1)) {
            int coopID = directCooperativeSpecies.get(speciesID1);
            if ((TFcoop.get(coopID).direction0 == Constants.NONE || TFcoop.get(coopID).direction0 == direction1) && (TFcoop.get(coopID).direction1 == Constants.NONE || TFcoop.get(coopID).direction1 == direction0)) {
                result = TFcoop.get(coopID);
            }
        }

        return result;
    }

    /**
     * generates the string description of the cooperativity
     */
    private String getCoopString() {
        StringBuilder str = new StringBuilder();
        if (this.TFcoop == null || this.TFcoop.isEmpty()) {
            str = new StringBuilder("none");
        } else {
            for (TFcooperativity tFcooperativity : this.TFcoop) {
                str.append("; ").append(tFcooperativity.toString());
            }
            str = new StringBuilder(str.substring(2));
        }
        return str.toString();
    }

}

