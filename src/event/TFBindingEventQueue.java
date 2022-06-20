package event;

import environment.Cell;

import java.io.Serializable;

/**
 * class that contains the list of binding events (3D diffusion)
 *
 * @author n.r.zabet@gen.cam.ac.uk
 * @author fedor.garbuzov@mail.ioffe.ru
 */
public class TFBindingEventQueue implements Serializable {

    private static final long serialVersionUID = 2270248942901610094L;
    protected double[] proteinBindingPropensity;
    protected double proteinBindingPropensitySum;
    private ProteinEvent bindingEvent;

    public TFBindingEventQueue(Cell n) {
        this.bindingEvent = null;
        this.proteinBindingPropensitySum = 0;
        this.proteinBindingPropensity = new double[n.getNoOfDBPspecies()];
        for (int i = 0; i < proteinBindingPropensity.length; i++) {
            this.proteinBindingPropensity[i] = computePropensity(i, n);
            this.proteinBindingPropensitySum += this.proteinBindingPropensity[i];
        }
        this.updateProteinBindingPropensities(n);
    }

    /**
     * returns the next binding event
     */
    public ProteinEvent peek() {
        return bindingEvent;
    }

    /**
     * returns the next binding event and removes it from the list
     */
    public ProteinEvent pop() {
        ProteinEvent pe = bindingEvent;
        bindingEvent = null;
        return pe;
    }

    /**
     * replaces the current binding event with a new one
     *
     * @param pe the new event
     */
    public void add(ProteinEvent pe) {
        bindingEvent = pe;
    }


    /**
     * returns true if the binding event is null
     */
    public boolean isEmpty() {
        return bindingEvent == null;
    }

    /**
     * computes the sum of all propensities as in Gillespie algorithm
     */
    private double propensitySum() {
        double result = 0;
        for (double d : proteinBindingPropensity) {
            result += d;
        }
        return result;
    }

    /**
     * returns the propensity sum
     */
    public double getPropensitySum() {
        return this.propensitySum();
    }

    /**
     * returns the array of all propensities
     */
    public double[] getPropensities() {
        return this.proteinBindingPropensity;
    }

    /**
     * deletes current protein binding event
     */
    public void clear() {
        this.bindingEvent = null;
    }

    /**
     * updates propensities once a TF abundance has changed
     */
    public void updateProteinBindingPropensities(Cell n) {
        if (n.isInDebugMode()) {
            n.printDebugInfo("Full update of propensities for TF binding");
        }
        this.proteinBindingPropensitySum = 0;
        for (int i = 0; i < proteinBindingPropensity.length; i++) {
            this.proteinBindingPropensity[i] = computePropensity(i, n);
            this.proteinBindingPropensitySum += this.proteinBindingPropensity[i];
        }
        //when the propensities are updated the binding is deleted in order to force its regeneration in the simulator
        this.clear();
    }


    /**
     * computes the propensity of binding a TF to the DNA
     *
     * @param TFspeciesID the id of the specie for which the propensity is computed
     * @param n           the cell
     */
    public double computePropensity(int TFspeciesID, Cell n) {
        return (n.TFspecies[TFspeciesID].assocRate * n.freeTFmolecules.get(TFspeciesID).size()
                * n.dna.effectiveTFavailabilitySum[TFspeciesID] / n.dna.effectiveTFavailabilityMaxSum[TFspeciesID]);
    }

}
