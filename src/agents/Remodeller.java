package agents;

import environment.Cell;
import event.RepressionEvent;
import utils.Constants;
import utils.Pair;
import utils.RepressionData;

import java.io.Serializable;

/**
 * Remodeller closes and opens chromatin (DNA).
 * Closing and opening of DNA (repression or derepression) is invoked by repressors (TFs).
 *
 * @author fedor.garbuzov@mail.ioffe.ru
 */
public class Remodeller implements Serializable {

    private static final long serialVersionUID = 7085739298818754958L;
    protected final double derepressionRate; // the rate at which the DNA is opened once repressor is unbound
    private double totalRepressionAndDerepressionEvents;

    /**
     * constructor
     *
     * @param derepressionRate the rate at which the DNA is opened once repressor is unbound
     */
    public Remodeller(double derepressionRate) {
        this.derepressionRate = derepressionRate;
        this.totalRepressionAndDerepressionEvents = 0;
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
                n.dna.repress(re.boundaryLeft, re.boundaryRight, n);
                // update repression state
                n.dbp[re.proteinID].repressesDNA = true;
                n.dbp[re.proteinID].updateRepressionRate(n);
                if (n.isInDebugMode()) {
                    n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name
                            + " started repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
                }
                n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(n);
                n.TFspecies[speciesID].countTFRepressionEvents++;
                totalRepressionAndDerepressionEvents++;
            } else if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name
                        + " at position " + n.dbp[re.proteinID].getPosition()
                        + " attempted to start repression but is repressed.");
            }
        } else {
            assert (re.nextAction == Constants.EVENT_TF_DEREPRESSION);
            // this prevents changing repression status if this derepression event was scheduled on the unbinding
            if (re.scheduleNextEvent) {
                n.dbp[re.proteinID].repressesDNA = false;
            }
            // update repression rate
            n.dbp[re.proteinID].updateRepressionRate(n);

            // the status message is shown before actual derepressing (DNA.derepress function) because
            // during derepressing new messages can be shown which should be after this one
            if (n.isInDebugMode()) {
                n.printDebugInfo(re.time + ": TF " + re.proteinID + " of type " + n.TFspecies[speciesID].name
                        + " ended repression from position " + re.boundaryLeft + " to " + re.boundaryRight);
            }

            n.dna.derepress(n, re.boundaryLeft, re.boundaryRight, re.proteinID, re.time);

            n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(n);
            n.TFspecies[speciesID].countTFDerepressionEvents++;
            totalRepressionAndDerepressionEvents++;
        }

        if (n.ip.OUTPUT_REPRESSED_LENGTHS.value) {
            // completely recalculate repression data to prevent accumulation of the rounding error
            if (totalRepressionAndDerepressionEvents % Constants.UPDATE_REPRESSION_SCORES_EVERY == 0) {
                n.dna.repressedRepScore = 0.0;
                n.dna.repressedActScore = 0.0;
                for (int pos = 0; pos < n.dna.strand.length; pos++) {
                    if (n.dna.closed[pos] == Constants.BP_IS_REPRESSED) {
                        n.dna.modifyRepressionScore(n, pos, true);
                    }
                }
            }
            // save repression data
            n.dna.repressionData.add(new RepressionData(re.time, n.dna.currentRepressedLength,
                    n.dna.repressedRepScore, n.dna.repressedActScore));
        }
    }

}
