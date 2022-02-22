package agents;

import environment.Cell;
import event.RepressionEvent;
import utils.Constants;

import java.io.Serializable;

public class Remodeller implements Serializable {

    private static final long serialVersionUID = 7085739298818754958L;
    protected double unrepressionRate;

    public Remodeller(double unrepressionRate) {
        this.unrepressionRate = unrepressionRate;
    }

    public double getRepressionRate() {
        return unrepressionRate;
    }

    public void act(Cell n, RepressionEvent re) {
        assert re.proteinID != Constants.NONE;
        int speciesID = n.dbp[re.proteinID].speciesID;
        if (re.nextAction == Constants.EVENT_TF_REPRESSION){
            if (!n.dbp[re.proteinID].isRepressed(n)) {
                n.dna.repress(n, re.boundaryLeft, re.boundaryRight);
                // FG: update repression state
                n.dbp[re.proteinID].repressesDNA = true;
                n.dbp[re.proteinID].updateRepressionRate(n);
                n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(Constants.NONE, n);
                n.TFspecies[speciesID].countTFRepressionEvents++;
                if (n.isInDebugMode()) {
                    n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name + " started repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
                }
            } else if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name + " at position " + n.dbp[re.proteinID].getPosition() + " attempted to start repression but is repressed.");
            }
        } else {
            assert (re.nextAction == Constants.EVENT_TF_UNREPRESSION);
            // FG: update repression rate
            n.dbp[re.proteinID].repressesDNA = false;
            n.dbp[re.proteinID].updateRepressionRate(n);

            n.dna.unrepress(n, re.boundaryLeft, re.boundaryRight);
            n.TFspecies[speciesID].countTFUnrepressionEvents++;

            // find all repressors in the repressed region and apply repression in their repression regions
            int boundMoleculeID, boundSpeciesID;
            int size, boundaryLeft, boundaryRight;
            int bpStart = n.dna.updateLeftBoundary(re.boundaryLeft - n.maxRepressionLeftSize);
            int bpEnd = n.dna.updateRightBoundary(re.boundaryRight + n.maxRepressionRightOrTFSize);
            for (int bpIdx = bpStart; bpIdx <= bpEnd; bpIdx++) {
                boundMoleculeID = n.dna.getBoundMolecule(bpIdx);
                if (boundMoleculeID != Constants.NONE) {
                    if (n.dbp[boundMoleculeID].isRepressingDNA() && boundMoleculeID != re.proteinID) {
                        boundSpeciesID = n.dbp[boundMoleculeID].speciesID;
                        size = n.TFspecies[boundSpeciesID].sizeTotal;
                        boundaryLeft = n.dna.updateLeftBoundary(bpIdx - n.TFspecies[boundSpeciesID].repressionLeftSize);
                        boundaryRight = n.dna.updateRightBoundary(bpIdx + size - 1 + n.TFspecies[boundSpeciesID].repressionRightSize);
                        n.dna.repressDNA(n, boundaryLeft, boundaryRight);
                    }
                }
            }

            n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(Constants.NONE, n);

            if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name + " ended repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
            }
        }
    }
}
