package agents;

import environment.Cell;
import event.RepressionEvent;
import utils.Constants;
import utils.Pair;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * This class specifies Remodeller which closes or opens chromatin (DNA).
 * Closing and opening of DNA (repression or derepression) is invoked by repressors (TFs).
 *
 * @author fedor.garbuzov@mail.ioffe.ru
 */
public class Remodeller implements Serializable {

    private static final long serialVersionUID = 7085739298818754958L;
    protected double derepressionRate; // the rate at which the DNA is opened once repressor is unbound

    /**
     * constructor
     *
     * @param derepressionRate the rate at which the DNA is opened once repressor is unbound
     */
    public Remodeller(double derepressionRate) {
        this.derepressionRate = derepressionRate;
    }

    public boolean derepressionIsInTheRegionOfMolecule(Cell n, RepressionEvent re) {
        int speciesID = n.dbp[re.proteinID].speciesID;
        return (re.boundaryLeft + n.TFspecies[speciesID].repressionLeftSize) == n.dbp[re.proteinID].position;
    }

    /**
     * performs DNA closing (repression) or opening (derepression)
     *
     * @param n  the environment
     * @param re repression or derepression event
     */
    public void act(Cell n, RepressionEvent re) {
        assert re.proteinID != Constants.NONE;
        int speciesID = n.dbp[re.proteinID].speciesID;
        if (re.nextAction == Constants.EVENT_TF_REPRESSION) {
            assert !n.dbp[re.proteinID].isRepressingDNA();
            if (!n.dbp[re.proteinID].isRepressed(n)) {
                // close chromatin
                n.dna.repress(re.boundaryLeft, re.boundaryRight);
                // update repression state
                n.dbp[re.proteinID].repressesDNA = true;
                n.dbp[re.proteinID].updateRepressionRate(n);
                if (n.isInDebugMode()) {
                    n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name +
                            " started repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
                }
                n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(n);
                n.TFspecies[speciesID].countTFRepressionEvents++;
                n.dna.repressedLength.add(new Pair<>(re.time, n.dna.currentRepressedLength));
            } else if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name + " at " +
                        "position " + n.dbp[re.proteinID].getPosition() + " attempted to start repression but is " +
                        "repressed.");
            }
        } else {
            assert (re.nextAction == Constants.EVENT_TF_DEREPRESSION);
            // FG: this prevents changing repression status if this derepression event was scheduled at the last unbinding
            if (re.scheduleNextEvent) {
                n.dbp[re.proteinID].repressesDNA = false;
            }
            // FG: update repression rate
            n.dbp[re.proteinID].updateRepressionRate(n);

            n.dna.derepress(n, re.boundaryLeft, re.boundaryRight, re.proteinID);
            n.TFspecies[speciesID].countTFDerepressionEvents++;
            n.dna.repressedLength.add(new Pair<>(re.time, n.dna.currentRepressedLength));

            if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name + " " +
                        "ended repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
            }

            n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(n);
        }
    }
}
