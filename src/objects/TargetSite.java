package objects;

import environment.Cell;
import utils.Constants;

import java.io.Serializable;
import java.util.ArrayList;


/**
 * a Target site object in a string of DNA
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
public class TargetSite implements Serializable {

    private static final long serialVersionUID = -4244069998461684773L;
    public int relStart; //position of the site on the DNA
    public int relEnd; //position of the site on the DNA
    public DNAregion region;
    public int TFid; // TF id
    public String TFname;
    public int size;
    public ArrayList<Integer> group;
    public int targetSiteID;


    /**
     * class constructor
     * @param n           the environment which contains all information about the biological system
     * @param targetSiteID the target site ID
     * @param chromosome  the name of the chromosome
     * @param start       the start position of the target site
     * @param end         the end position of the target site
     * @param dnaRegion   dna region which contains the TS
     * @param TFid        the ID of the TF that should bind
     * @param TFsize      the size of the TF
     */
    public TargetSite(Cell n, int targetSiteID, String chromosome, long start, long end, DNAregion dnaRegion,
                      int TFid, int TFsize, int DNAsize) {
        this.region = new DNAregion(chromosome, start, end);
        this.relStart = (int) (region.start - dnaRegion.start);
        this.relEnd = (int) (this.region.end - dnaRegion.start - TFsize);
        rescaleInterval(TFsize, DNAsize);
        this.TFid = TFid;
        this.TFname = n.TFspecies[TFid].name;
        this.size = Math.max(relEnd - relStart, 1);

        group = new ArrayList<Integer>();
        this.targetSiteID = targetSiteID;
    }


    /**
     * class constructor
     * @param targetSiteID the target site ID
     * @param description the text that defines the region
     * @param chromosome  the name of the chromosome
     * @param start       the start position of the target site
     * @param end         the end position of the target site
     * @param dnaRegion   dna region which contains the TS
     * @param TFid        the ID of the TF that should bind
     * @param TFsize      the size of the TF
     */
    public TargetSite(int targetSiteID, String description, String chromosome, long start, long end,
                      DNAregion dnaRegion, int TFid, int TFsize, int DNAsize) {
        this.region = new DNAregion(description, chromosome, start, end, false, true);
        this.relStart = (int) (this.region.start - dnaRegion.start);
        this.relEnd = (int) (this.region.end - dnaRegion.start);
        this.TFid = TFid;
        this.size = Math.max(relEnd - relStart - TFsize, 1);

        rescaleInterval(TFsize, DNAsize);
        group = new ArrayList<Integer>();
        this.targetSiteID = targetSiteID;
    }


    /**
     * class constructor
     *
     * @param n           the environment which contains all information about the biological system
     * @param targetSiteID the target site ID
     * @param description the text that defines the region
     * @param chromosome  the name of the chromosome
     * @param start       the start position of the target site
     * @param end         the end position of the target site
     */
    public TargetSite(Cell n, int targetSiteID, String description, String chromosome, long start, long end) {
        String DNAregionDescription = description, TFstr;
        int delimiterPos, TFsize = 0;

        if (description.contains(Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER)) {
            delimiterPos = description.indexOf(Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER);
            TFstr = description.substring(0, delimiterPos);
            DNAregionDescription = description.substring(delimiterPos + 1);

            this.TFid = n.getTFspeciesID(TFstr);
            if (this.TFid == Constants.NONE) {
                n.stopSimulation("error while parsing target site " + description + "; unkown TF");
            }
            this.TFname = n.TFspecies[TFid].name;
            TFsize = n.TFspecies[this.TFid].sizeTotal;
            this.region = new DNAregion(DNAregionDescription, chromosome, start, end, false, true);
            this.relStart = (int) (this.region.start - n.dna.subsequence.start);
            this.relEnd = (int) (this.region.end - n.dna.subsequence.start - TFsize);
            rescaleInterval(TFsize, n.dna.strand.length);

            this.size = Math.max(relEnd - relStart, 1);
            group = new ArrayList<Integer>();
            this.targetSiteID = targetSiteID;

        } else {
            n.stopSimulation("error while parsing target site " + description + "; no TF species");
        }
    }


    /**
     * returns a string with the location if the provided argument is the id of the TF that needs to bind at current
     * target
     */
    public String toString() {
        return this.TFname + Constants.FASTA_FILE_DESCRIPTION_CHROMOSOME_DELIMITER + this.region.toString();
    }


    private void rescaleInterval(int TFsize, int DNAsize) {
        if (relEnd <= relStart) {
            relEnd = relStart + 1;
        }
        relStart = Math.max(0, relStart);
        relStart = Math.min(relStart, DNAsize - TFsize);
        relEnd = Math.min(DNAsize - TFsize, relEnd);
        relEnd = Math.max(relEnd, 0);
    }

    /**
     * returns true if the two target sites start and end at the same position and have the same direction
     */
    public boolean equals(TargetSite ts) {
        return this.relStart == ts.relStart && this.relEnd == ts.relEnd && this.region.direction == ts.region.direction && this.TFid == ts.TFid;
    }

    /**
     * returns true if this target site is in an AND group
     */
    public boolean isInGroup() {
        return this.group.size() > 0;
    }

}
