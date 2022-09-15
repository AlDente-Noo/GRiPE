package event;

import utils.Constants;

/**
 * Class that describes a protein event. It is an instantiation of Event
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
public class ProteinEvent extends Event {

    private static final long serialVersionUID = 6542151547203073970L;
    public int proteinID;
    public int position;
    public boolean isTF;
    public boolean isHoppingEvent;
    public double propensity;

    public ProteinEvent() {
        super(Double.MAX_VALUE, Constants.NONE);
        this.proteinID = Constants.NONE;
        this.position = Constants.NONE;
        this.isTF = false;
        this.isHoppingEvent = false;
        propensity = 0.0;
    }

    /**
     * class constructor
     *
     * @param time      the time of the event
     * @param proteinID the id of the protein which is affected
     * @param position  the position on the DNA
     * @param isTF      true if this is a TF event or false otherwise
     */
    public ProteinEvent(double time, int proteinID, int position, boolean isTF, int nextAction,
                        boolean isHoppingEvent, double propensity) {
        super(time, nextAction);
        this.proteinID = proteinID;
        this.position = position;
        this.isTF = isTF;
        this.isHoppingEvent = isHoppingEvent;
        this.propensity = propensity;
    }

    /**
     * returns true if the current event is a binding one
     */
    public boolean isBindingEvent() {
        return position >= 0;
    }

    /**
     * generates the description string of current event
     */
    public String toString() {
        String stateStr = "" + time + ": ";
        stateStr += "TF";
        stateStr += " " + proteinID + " to position " + this.position;
        stateStr += " through an event of type " + nextAction;
        return stateStr;
    }

    /**
     * compares whether this event equals the one supplied as an argument
     */
    public boolean isEqualTo(ProteinEvent pe) {
        return this.proteinID == pe.proteinID && this.time == pe.time;
    }

    public boolean isEmpty() {
        return time == Double.MAX_VALUE && nextAction == Constants.NONE
                && proteinID == Constants.NONE && propensity == 0.0;
    }

}
