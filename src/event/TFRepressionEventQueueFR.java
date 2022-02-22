package event;

import java.util.PriorityQueue;

import environment.Cell;
import utils.Constants;
import utils.Gillespie;
import utils.Utils;

public class TFRepressionEventQueueFR extends TFEventQueue {

    private PriorityQueue<RepressionEvent> events;

    public TFRepressionEventQueueFR(Cell n) {
        events = new PriorityQueue<RepressionEvent>();
    }

    public void add(Event e) {
        events.add((RepressionEvent) e);
    }

    public RepressionEvent peek() {
        return events.peek();
    }

    public RepressionEvent pop() {
        return events.poll();
    }


    /**
     * returns true if the list of events is empty or false otherwise
     */
    public boolean isEmpty(){
        return events.isEmpty();
    }

    /**
     * returns the number of events in the list
     */
    public int size(){
        return events.size();
    }

    public RepressionEvent createNextEvent(Cell n, int moleculeID, double time) {
        if (n.dbp[moleculeID].getRepressionRate() < Constants.DOUBLE_ZERO) {
            return new RepressionEvent();
        }
        int speciesID = n.dbp[moleculeID].speciesID;
        int position = n.dbp[moleculeID].getPosition();
        assert position != Constants.NONE;
        int boundaryLeft  = n.dna.updateLeftBoundary(position - n.TFspecies[speciesID].repressionLeftSize);
        int boundaryRight = n.dna.updateRightBoundary(position + n.TFspecies[speciesID].sizeTotal - 1 + n.TFspecies[speciesID].repressionRightSize);
        double propensity = n.dbp[moleculeID].getRepressionRate();
        double nextTime = Gillespie.computeNextReactionTime(propensity, n.randomGenerator);
        int nextAction;
        if (n.dbp[moleculeID].isRepressingDNA()) {
            nextAction = Constants.EVENT_TF_UNREPRESSION;
        }
        else {
            nextAction = Constants.EVENT_TF_REPRESSION;
        }
        return new RepressionEvent(time + nextTime, nextAction, moleculeID, boundaryLeft, boundaryRight, propensity);
    }

    public void scheduleNextEvent(Cell n, int moleculeID, Event e) {
        RepressionEvent re = (RepressionEvent) e;
        if (!re.isEmpty()) {
            n.dbp[moleculeID].pe = null;
            n.dbp[moleculeID].re = re;
            this.add(re);
        }
    }


    public void updateNextEvent(Cell n, int moleculeID, double time) {
        boolean removed = this.events.remove(n.dbp[moleculeID].re);
        if(removed){
            RepressionEvent re = this.createNextEvent(n, moleculeID, time);
            this.scheduleNextEvent(n, moleculeID, re);
        }
    }

    //public double getNextEventTime(int moleculeID);
    
}
