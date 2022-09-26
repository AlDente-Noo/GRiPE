package utils;

import environment.Cell;
import objects.DNAregion;
import objects.TFSpecies;

import java.util.ArrayList;


/**
 * class that parses a TF description file
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
public class TFfileParser {
    public ArrayList<TFSpecies> data;
    public boolean loaded;
    public boolean hasVariousRowSizes;
    public boolean hasAllHeaders;
    public boolean parsed;


    public TFfileParser(Cell n, String filename, int defaultCopyNumber, double defaultEs, int defaultSizeLeft,
                        int defaultSizeRight, double defaultAssocRate, int DNAsize,
                        double defaultUnBindingProbability, double defaultSlideLeftProbability,
                        double defaultSlideRightProbability, double defaultJumpingProbability,
                        double defaultHopSTDdisplacement, double defaultSpecificWaitingTime, int defaultStepLeftSize,
                        int defaultStepRightSize, int defaultUncorrelatedDisplacementSize,
                        boolean defaultStallsIfBlocked, double defaultCollisionUnbindingProbability,
                        double defaultAffinityLandscapeRoughness, double defaultPreboundProportion,
                        boolean defaultPreboundToHighestAffinity, boolean defaultIsImmobile, DNAregion dnaRegion,
                        boolean deafultIsBiasedRandomWalk, boolean defaultIsTwoStateRandomWalk, double defaultSpecificEnergyThreshold,
                        double defaultRepressionRate, double defaultDerepressionAttenuationFactor, int defaultReprLenLeft,
                        int defaultReprLenRight) {
        parsed = false;
        data = new ArrayList<TFSpecies>();

        CSVparser csv = new CSVparser(filename, true, true);

        hasAllHeaders = false;
        // manage to parse csv
        this.loaded = csv.loaded;
        this.hasVariousRowSizes = csv.hasVariousRowSizes;

        if (csv.loaded && !csv.hasVariousRowSizes) {
            hasAllHeaders = true;
            // check header
            for (String key : Constants.PARSER_TF_CSV_FILE_HEADER) {
                if (!csv.header.containsKey(key)) {
                    hasAllHeaders = false;
                    n.printDebugInfo("Error when parsing the TF csv file " + filename + ": missing column " + key);
                }
            }

            if (!csv.data.isEmpty()) {
                parsed = true;
            }

            int id; // the id of the TF species
            String name; // the id of the TF species
            byte[] dbd; // the sequence recognise at DNA binding domain
            int copyNumber; // number of molecules of this species
            double es; // specific binding energy per nucleotide
            String bufferDBD;
            id = 0;
            String cellContent;

            int sizeLeft, sizeRight;
            double assocRate;
            DNAregion bufferDNAregion;
            boolean isCognate;

            double unBindingProbability;
            double slideLeftProbability;
            double slideRightProbability;
            double jumpingProbability;
            double hopSTDdisplacement;
            double specificWaitingTime;
            int stepLeftSize;
            int stepRightSize;
            int uncorrelatedDisplacementSize;
            boolean stallsIfBlocked;
            double collisionUnbindingProbability;
            double affinityLandscapeRoughness;
            double preboundProportion;
            boolean preboundToHighestAffinity;
            boolean isImmobile;
            boolean isBiasedRandomWalk;
            boolean isTwoStateRandomWalk;
            double specificEnergyThreshold;

            double repressionRate;
            double derepressionAttenuationFactor;
            int reprLenLeft;
            int reprLenRight;

            for (ArrayList<String> buffer : csv.data) {
                //init params
                name = "";
                bufferDBD = "";
                dbd = new byte[0];
                copyNumber = defaultCopyNumber;
                es = defaultEs;
                sizeLeft = defaultSizeLeft;
                sizeRight = defaultSizeRight;
                assocRate = defaultAssocRate;
                bufferDNAregion = new DNAregion("", "", 0, DNAsize, true, false);
                unBindingProbability = defaultUnBindingProbability;
                slideLeftProbability = defaultSlideLeftProbability;
                slideRightProbability = defaultSlideRightProbability;
                jumpingProbability = defaultJumpingProbability;
                hopSTDdisplacement = defaultHopSTDdisplacement;
                specificWaitingTime = defaultSpecificWaitingTime;
                stepLeftSize = defaultStepLeftSize;
                stepRightSize = defaultStepRightSize;
                uncorrelatedDisplacementSize = defaultUncorrelatedDisplacementSize;
                stallsIfBlocked = defaultStallsIfBlocked;
                collisionUnbindingProbability = defaultCollisionUnbindingProbability;
                affinityLandscapeRoughness = defaultAffinityLandscapeRoughness;
                preboundProportion = defaultPreboundProportion;
                preboundToHighestAffinity = defaultPreboundToHighestAffinity;
                isImmobile = defaultIsImmobile;
                isBiasedRandomWalk = deafultIsBiasedRandomWalk;
                isTwoStateRandomWalk = defaultIsTwoStateRandomWalk;
                specificEnergyThreshold = defaultSpecificEnergyThreshold;
                repressionRate = defaultRepressionRate;
                derepressionAttenuationFactor = defaultDerepressionAttenuationFactor;
                reprLenLeft = defaultReprLenLeft;
                reprLenRight = defaultReprLenRight;

                //NAME
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[0])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[0]));
                    name = cellContent;
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF name at line: " + buffer);
                }

                //DBD
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[1])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[1]));
                    bufferDBD = cellContent;
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF DNA binding domain at line: " + buffer);
                }

                //ES
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[2])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[2]));
                    es = Utils.parseDouble(cellContent, es);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF ES (energy penalty for a nucleotide " +
                            "missmatch) at line: " + buffer);
                }

                //COPYNUMBER
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[3])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[3]));
                    copyNumber = Utils.parseInteger(cellContent, copyNumber);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF copy number at line: " + buffer);
                }

                //TF_SIZE_LEFT
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[4])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[4]));
                    sizeLeft = Utils.parseInteger(cellContent, sizeLeft);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF size left at line: " + buffer);
                }

                //TF_SIZE_RIGHT
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[5])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[5]));
                    sizeRight = Utils.parseInteger(cellContent, sizeRight);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF size right at line: " + buffer);
                }

                //TF_ASSOC_RATE
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[6])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[6]));
                    assocRate = Utils.parseDouble(cellContent, assocRate);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF association rate at line: " + buffer);
                }

                //INITIAL DROP
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[7])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[7]));
                    bufferDNAregion = new DNAregion(cellContent, "", 0, DNAsize, true, false);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF initial drop interval at line: " + buffer);
                }

                //IS_COGNATE
                isCognate = dbd != null && dbd.length > 0;

                //unBindingProbability
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[8])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[8]));
                    unBindingProbability = Utils.parseDouble(cellContent, unBindingProbability);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF unbinding probability at line: " + buffer);
                }
                //slideLeftProbability
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[9])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[9]));
                    slideLeftProbability = Utils.parseDouble(cellContent, slideLeftProbability);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF slide left probability at line: " + buffer);
                }
                //slideRightProbability
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[10])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[10]));
                    slideRightProbability = Utils.parseDouble(cellContent, slideRightProbability);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF slide right probability at line: " + buffer);
                }
                //jumpingProbability
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[11])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[11]));
                    jumpingProbability = Utils.parseDouble(cellContent, jumpingProbability);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF jumping probability at line: " + buffer);
                }

                //hopSTDdisplacement
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[12])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[12]));
                    hopSTDdisplacement = Utils.parseDouble(cellContent, hopSTDdisplacement);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF hop standard displacement at line: " + buffer);
                }

                //specificWaitingTime
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[13])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[13]));
                    specificWaitingTime = Utils.parseDouble(cellContent, specificWaitingTime);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF specific waiting time at line: " + buffer);
                }

                //stepLeftSize
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[14])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[14]));
                    stepLeftSize = Utils.parseInteger(cellContent, stepLeftSize);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF step left size at line: " + buffer);
                }

                //stepRightSize
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[15])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[15]));
                    stepRightSize = Utils.parseInteger(cellContent, stepRightSize);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF step right size at line: " + buffer);
                }

                //uncorrelatedDisplacementSize
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[16])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[16]));
                    uncorrelatedDisplacementSize = Utils.parseInteger(cellContent, uncorrelatedDisplacementSize);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF uncorrelated displacement size at line: " + buffer);
                }

                //stallsIfBlocked
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[17])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[17]));
                    stallsIfBlocked = Utils.parseBoolean(cellContent, stallsIfBlocked);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF stalls if blocked at line: " + buffer);
                }

                //collision unbinding probability
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[18])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[18]));
                    collisionUnbindingProbability = Utils.parseDouble(cellContent, collisionUnbindingProbability);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF collision unbinding probability at line: " + buffer);
                }

                //affinityLandscapeRoughness
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[19])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[19]));
                    affinityLandscapeRoughness = Utils.parseDouble(cellContent, affinityLandscapeRoughness);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF affinity landscape roughness at line: " + buffer);
                }

                //preboundProportion
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[20])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[20]));
                    preboundProportion = Utils.parseDouble(cellContent, preboundProportion);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF prebound proportion at line: " + buffer);
                }

                //preboundToHighestAffinity
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[21])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[21]));
                    preboundToHighestAffinity = Utils.parseBoolean(cellContent, preboundToHighestAffinity);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF preboundToHighestAffinity at line: " + buffer);
                }

                //isImmobile
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[22])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[22]));
                    isImmobile = Utils.parseBoolean(cellContent, isImmobile);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF isImmobile at line: " + buffer);
                }

                //is biased random walk
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[23])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[23]));
                    isBiasedRandomWalk = Utils.parseBoolean(cellContent, isBiasedRandomWalk);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF is biased random walk at line: " + buffer);
                }

                //is two state random walk
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[24])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[24]));
                    isTwoStateRandomWalk = Utils.parseBoolean(cellContent, isTwoStateRandomWalk);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF is two state random walk at line: " + buffer);
                }

                // FG: repression rates and lengths
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[25])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[25]));
                    repressionRate = Utils.parseDouble(cellContent, repressionRate);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF repression rate at line: " + buffer);
                }
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[26])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[26]));
                    derepressionAttenuationFactor = Utils.parseDouble(cellContent, derepressionAttenuationFactor);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF derepression rate at line: " + buffer);
                }
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[27])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[27]));
                    reprLenLeft = Utils.parseInteger(cellContent, reprLenLeft);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF repression len left at line: " + buffer);
                }
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[28])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[28]));
                    reprLenRight = Utils.parseInteger(cellContent, reprLenRight);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses TF repression len right at line: " + buffer);
                }

                //FG
                if (csv.header.containsKey(Constants.PARSER_TF_CSV_FILE_HEADER[29])) {
                    cellContent = buffer.get(csv.header.get(Constants.PARSER_TF_CSV_FILE_HEADER[29]));
                    specificEnergyThreshold = Utils.parseDouble(cellContent, specificEnergyThreshold);
                } else {
                    n.printDebugInfo("TF file " + filename + " misses specific PWM threshold at line: " + buffer);
                }

                if (copyNumber > 0) {
                    data.add(new TFSpecies(dnaRegion, id, name, bufferDBD, copyNumber, es, sizeLeft, sizeRight,
                            assocRate, bufferDNAregion, isCognate, unBindingProbability, slideLeftProbability,
                            slideRightProbability, jumpingProbability, hopSTDdisplacement, specificWaitingTime,
                            stepLeftSize, stepRightSize, uncorrelatedDisplacementSize, stallsIfBlocked,
                            collisionUnbindingProbability, affinityLandscapeRoughness, preboundProportion,
                            preboundToHighestAffinity, isImmobile, isBiasedRandomWalk, isTwoStateRandomWalk,
                            specificEnergyThreshold,
                            repressionRate, derepressionAttenuationFactor, reprLenLeft, reprLenRight, n));
                    id++;
                }

            }

        } else {
            n.printDebugInfo("Could not load the TF csv file: '" + filename + "'");
        }

    }


}
