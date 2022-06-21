package objects;

import environment.Cell;
import utils.Constants;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * contains a list of target sites and a list of the target sites groups
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
public class TargetSitesAndGroups implements Serializable {

    private static final long serialVersionUID = -8051357555423486023L;
    public ArrayList<TargetSite> ts;
    public ArrayList<TargetSitesGroup> tsg;
    public boolean[] occupancy;

    /**
     * class constructor
     */
    public TargetSitesAndGroups() {
        ts = new ArrayList<TargetSite>();
        tsg = new ArrayList<TargetSitesGroup>();
    }

    /**
     * this is true if no TS is created
     */
    public boolean isEmpty() {
        return ts.isEmpty();
    }


    /**
     * adds a target site i
     * if the target site does not exists is added to the list and its ID is returned otherwise the ID of the
     * existing target site is returned;
     *
     * @param ts the target site
     */
    public int addTargetSite(TargetSite ts) {
        int result = Constants.NONE;

        for (int i = 0; i < this.ts.size() && result == Constants.NONE; i++) {
            if (this.ts.get(i).equals(ts)) {
                result = i;
            }
        }

        if (result == Constants.NONE) {
            ts.targetSiteID = this.ts.size();
            this.ts.add(ts);
            occupancy = new boolean[this.ts.size()];
            Arrays.fill(occupancy, false);
            result = this.ts.size() - 1;
        }

        return result;

    }


    /**
     * generate the text string
     */

    public String toString() {
        StringBuilder text = new StringBuilder();
        for (TargetSitesGroup targetSitesGroup : tsg) {
            text.append("\"").append(targetSitesGroup.toString()).append("\", ").append(targetSitesGroup.firstTimeReached).append(", ").append(targetSitesGroup.timesReached).append(", ").append(targetSitesGroup.timeOccupied).append("\n");
        }
        return text.toString();
    }


    /**
     * adds a new group
     *
     * @param n    the cell
     * @param pos  the position on the DNA of the target site
     * @param TFid the ID of the TF
     */
    public void addGroup(Cell n, int pos, int TFid) {
        int newID = this.addTargetSite(new TargetSite(n, this.ts.size(), "", pos, pos + 1, n.dna.region, TFid,
                n.TFspecies[TFid].sizeTotal, n.dna.strand.length));
        TargetSitesGroup bufferTSG = new TargetSitesGroup(tsg.size(), ts.get(newID).toString());
        bufferTSG.addTargetSite(newID);
        bufferTSG.generateRPN(newID + "");
        bufferTSG.text = this.ts.get(newID).toString();
        tsg.add(bufferTSG);
    }

    /**
     * updates statistics about a target site
     *
     * @param tsID the ID of the target site
     * @param time the current simulation time
     */
    public void updateTargetSiteStatistics(int tsID, double time, boolean bound) {

        this.occupancy[tsID] = bound;

        //check for each group the current TF belongs to the state
        boolean evaluateGroup;
        for (int i : ts.get(tsID).group) {
            evaluateGroup = tsg.get(i).evaluateRPNTree(occupancy);
            if (evaluateGroup && !tsg.get(i).isOccupied) {
                tsg.get(i).updateTimesReachedStatistics(time);
            }

            if (tsg.get(i).isOccupied) {
                tsg.get(i).updateOccupancyStatistics(time);
            }

            tsg.get(i).isOccupied = evaluateGroup;
            tsg.get(i).lastTimeUpdate = time;
        }

    }


    /**
     * verifies if there are target sites to be reached
     *
     * @return true if there are unreached site and false otherwise
     */
    public boolean areTargetSitesToBeReached() {
        boolean result = false;

        for (int tsID = 0; tsID < ts.size() && !result; tsID++) {
            for (int i = 0; i < ts.get(tsID).group.size() && !result; i++) {
                if (tsg.get(i).timesReached == 0) {
                    result = true;
                }
            }
        }

        return result;

    }

    /**
     * returns a string with the names of all target site groups
     */
    public String getTargetSiteGroupsString() {
        StringBuilder text = new StringBuilder();
        for (int i = 0; i < tsg.size(); i++) {
            if (i > 0) {
                text.append(", ");
            }
            text.append("\"").append(tsg.get(i).text).append("\"");
        }
        return text.toString();
    }

    /**
     * returns a string with the occupancy of all target site groups
     */
    public String getTargetSiteGroupsOccupancyString() {
        StringBuilder text = new StringBuilder();
        for (int i = 0; i < this.tsg.size(); i++) {
            if (i > 0) {
                text.append(", ");
            }
            if (this.tsg.get(i).isOccupied) {
                text.append(" 1");
            } else {
                text.append(" 0");
            }
        }
        return text.toString();
    }
}
