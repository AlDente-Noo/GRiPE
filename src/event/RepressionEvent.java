package event;

import utils.Constants;

public class RepressionEvent extends Event {
    public int proteinID;
    public int boundaryLeft, boundaryRight;
    public double propensity;
    public boolean scheduleNextEvent;

    public RepressionEvent() {
        super(Double.MAX_VALUE, Constants.NONE);
        proteinID = Constants.NONE;
        boundaryRight = Constants.NONE;
        boundaryLeft = Constants.NONE;
        propensity = 0.0;
        scheduleNextEvent = false;
    }

    /**
     * class constructor
     */
    public RepressionEvent(double time, int nextAction, int proteinID, int boundaryLeft, int boundaryRight,
                           double propensity, boolean scheduleNextEvent) {
        super(time, nextAction);
        this.proteinID = proteinID;
        this.boundaryLeft = boundaryLeft;
        this.boundaryRight = boundaryRight;
        this.propensity = propensity;
        this.scheduleNextEvent = scheduleNextEvent;
    }

    public String toString() {
        String stateStr = "" + time + ": ";
        String actionStr;
        assert (nextAction == Constants.EVENT_TF_REPRESSION || nextAction == Constants.EVENT_TF_DEREPRESSION);
        if (nextAction == Constants.EVENT_TF_REPRESSION) {
            actionStr = "Repression";
        } else {
            actionStr = "Derepression";
        }
        stateStr += actionStr + " by TF " + proteinID + " from position " + boundaryLeft + " to " + boundaryRight;
        return stateStr;
    }

    public boolean isEmpty() {
        return time == Double.MAX_VALUE && nextAction == Constants.NONE && proteinID == Constants.NONE && propensity == 0.0;
    }
}
