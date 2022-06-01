package agents;

import environment.Cell;
import event.RepressionEvent;
import utils.Constants;

import java.io.Serializable;

/**
 * This class specifies Remodeller which closes or opens chromatin (DNA).
 * Closing and opening of DNA (repression or 'unrepression') is invoked by repressor TFs.
 *
 * @author fedor.garbuzov@mail.ioffe.ru
 */
public class Remodeller implements Serializable {

    private static final long serialVersionUID = 7085739298818754958L;
    protected double unrepressionRate; // the rate at which the DNA is opened once repressor is unbound

    /**
     * constructor
     *
     * @param unrepressionRate the rate at which the DNA is opened once repressor is unbound
     */
    public Remodeller(double unrepressionRate) {
        this.unrepressionRate = unrepressionRate;
    }

    /**
     * performs DNA closing (repression) or opening ('unrepression')
     *
     * @param n  pointer to the environment
     * @param re repression or 'unrepression' event
     */
    public void act(Cell n, RepressionEvent re) {
        assert re.proteinID != Constants.NONE;
        int speciesID = n.dbp[re.proteinID].speciesID;
        if (re.nextAction == Constants.EVENT_TF_REPRESSION) {
            if (!n.dbp[re.proteinID].isRepressed(n)) {
                // close chromatin
                n.dna.repress(n, re.boundaryLeft, re.boundaryRight);
                // update repression state
                n.dbp[re.proteinID].repressesDNA = true;
                n.dbp[re.proteinID].updateRepressionRate(n);
                n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(Constants.NONE, n);
                n.TFspecies[speciesID].countTFRepressionEvents++;
                if (n.isInDebugMode()) {
                    n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name +
                            " started repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
                }
            } else if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name + " at " +
                        "position " + n.dbp[re.proteinID].getPosition() + " attempted to start repression but is " +
                        "repressed.");
            }
        } else {
            assert (re.nextAction == Constants.EVENT_TF_UNREPRESSION);
            // change repression status if only this unrepression event is for the nearby region
            // this prevents changing repression status if this unrepression event was scheduled at the last unbinding
            if ((re.boundaryLeft + n.TFspecies[speciesID].repressionLeftSize) == n.dbp[re.proteinID].position) {
                n.dbp[re.proteinID].repressesDNA = false;
            }
            // FG: update repression rate
            n.dbp[re.proteinID].updateRepressionRate(n);

            n.dna.unrepress(n, re.boundaryLeft, re.boundaryRight);
            n.TFspecies[speciesID].countTFUnrepressionEvents++;

            // find all repressors in the repressed region and close chromatin in their repression regions
            int boundMoleculeID, boundSpeciesID;
            int size, boundaryLeft, boundaryRight;
            int bpStart = n.dna.updateLeftBoundary(re.boundaryLeft - n.maxRepressionLeftSize - n.maxTFSize);
            int bpEnd = n.dna.updateRightBoundary(re.boundaryRight + n.maxRepressionRightOrTFSize);
            for (int bpIdx = bpStart; bpIdx <= bpEnd; bpIdx++) {
                boundMoleculeID = n.dna.getBoundMolecule(bpIdx);
                if (boundMoleculeID != Constants.NONE) {
                    if (n.dbp[boundMoleculeID].isRepressingDNA() && boundMoleculeID != re.proteinID) {
                        boundSpeciesID = n.dbp[boundMoleculeID].speciesID;
                        size = n.TFspecies[boundSpeciesID].sizeTotal;
                        boundaryLeft = n.dna.updateLeftBoundary(bpIdx - n.TFspecies[boundSpeciesID].repressionLeftSize);
                        boundaryRight =
                                n.dna.updateRightBoundary(bpIdx + size - 1 + n.TFspecies[boundSpeciesID].repressionRightSize);
                        n.dna.repressDNA(n, boundaryLeft, boundaryRight);
                    }
                }
            }

            n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(Constants.NONE, n);

            if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name + " " +
                        "ended repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
            }
        }
    }
}
