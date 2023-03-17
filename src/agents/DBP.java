package agents;

import environment.Cell;
import event.ProteinEvent;
import event.RepressionEvent;
import utils.Constants;
import utils.Utils;

import java.io.Serializable;

/**
 * this class specifies all DNA Binding Proteins
 *
 * @author n.r.zabet@gen.cam.ac.uk (original contributor)
 * @author fedor.garbuzov@mail.ioffe.ru (modified the code)
 */
public abstract class DBP implements Serializable {
    private static final long serialVersionUID = 121399550800359214L;
    public int ID;  // the ID of the TF
    public int size; // = sizeTotal from TFSpecies
    public int speciesID; //the ID of the species to which this TF belongs
    public int TAU;
    public boolean hasDNAbasedCooperativity;
    public boolean hasDirectCooperativity;
    public ProteinEvent pe;
    public RepressionEvent re;
    public double moveRate;
    public double repressionEventRate;
    public int leftNeighbour;
    public int rightNeighbour;
    public int stickToLeft;
    public int stickToRight;
    public int leftMostPosition;
    public int rightMostPosition;
    public int observedLeftMostPosition;
    public int observedRightMostPosition;
    public int slidingEvents;
    public int boundToDNA;
    protected int position; // the position on the DNA  (-1 free)
    protected int lastPosition; // the last position of this TF
    protected double timeOfLastPositionChange; // the last time the position was changed
    protected boolean repressesDNA; // FG: if true, then the TF represses a DNA region
    protected int direction;// the direction on which it is bound 0 for 5' -> 3' or 1 for 3' -> 5'
    protected double timeBound;
    protected boolean wasBound;


    /**
     * empty class constructor
     */
    public DBP() {
        ID = Constants.NONE;
        position = Constants.NONE;
        lastPosition = Constants.NONE;
        timeOfLastPositionChange = 0;
        size = 0;
        timeBound = 0;
        moveRate = Constants.NONE;
        repressionEventRate = Constants.NONE;

        leftNeighbour = Constants.NONE;
        rightNeighbour = Constants.NONE;
        pe = null;
        re = null;
        stickToLeft = Constants.NONE;
        stickToRight = Constants.NONE;

        leftMostPosition = Constants.NONE;
        rightMostPosition = Constants.NONE;
        observedLeftMostPosition = Constants.NONE;
        observedRightMostPosition = Constants.NONE;
        slidingEvents = 0;
        wasBound = false;
        boundToDNA = Constants.NONE;
    }

    /**
     * class constructor
     *
     * @param ID                       the ID of the TF
     * @param position                 the position on the DNA  (-1 free)
     * @param lastPosition             the last position of this TF
     * @param timeOfLastPositionChange the last time the position was changed
     * @param size                     total size of TF (motif size + left and right size)
     */
    public DBP(int ID, int position, int lastPosition, int timeOfLastPositionChange, int size) {
        this.ID = ID;
        this.position = position;
        this.lastPosition = lastPosition;
        this.timeOfLastPositionChange = timeOfLastPositionChange;
        this.size = size;
        timeBound = 0;
        moveRate = Constants.NONE;
        repressionEventRate = Constants.NONE;
        leftNeighbour = Constants.NONE;
        rightNeighbour = Constants.NONE;
        pe = null;
        re = null;
        stickToLeft = Constants.NONE;
        stickToRight = Constants.NONE;
        leftMostPosition = Constants.NONE;
        rightMostPosition = Constants.NONE;
        observedLeftMostPosition = Constants.NONE;
        observedRightMostPosition = Constants.NONE;
        wasBound = false;
        boundToDNA = Constants.NONE;
    }

    public abstract void act(Cell n, ProteinEvent pe);

    public abstract void updateBoundTime(Cell n, double timeBound, int direction, int position);

    /**
     * returns current position
     */
    public int getPosition() {
        return this.position;
    }

    /** FG
     * returns true if TF is currently repressing DNA
     */
    public boolean isRepressingDNA() {
        return repressesDNA;
    }

    /** FG
     * returns true if TF is repressed by another TF
     */
    public boolean isRepressed(Cell n) {
        // free TF cannot be repressed
        if (this.position == Constants.NONE) {
            return false;
        }
        // it is assumed that TF is repressed if at least one of the bps of TF site is under repression
        for (int i = 0; i < this.size; i++) {
            if (n.dna.closed[this.position + i] == Constants.BP_IS_REPRESSED) {
                return true;
            }
        }
        return false;
    }

    public abstract void setPosition(int newPosition, double timeOfLastPositionChange);

    public abstract void clearPosition(double timeOfLastPositionChange);

    /**
     * sets the sliding extremes
     */
    public void initSlidingExtremes() {
        this.leftMostPosition = position;
        this.rightMostPosition = position;
        this.slidingEvents = 0;
    }

    /**
     * sets the sliding extremes when sliding to left
     */
    public void setLeftSlidingExtreme() {
        if (position < this.leftMostPosition) {
            this.leftMostPosition = position;
        }
        slidingEvents++;
    }

    /**
     * sets the sliding extremes when sliding to right
     */
    public void setRightSlidingExtreme() {
        if (position > this.rightMostPosition) {
            this.rightMostPosition = position;
        }
        slidingEvents++;
    }

    /**
     * increments sliding events
     */
    public void incrementSlidingEvents() {
        slidingEvents++;
    }

    /**
     * returns the current slided length
     */
    public int getSlidingLength() {
        return this.rightMostPosition - this.leftMostPosition;
    }

    /**
     * returns the sliding events
     */
    public int getSlidingEvents() {
        return slidingEvents;
    }

    /**
     * sets the sliding extremes
     */
    public void initObservedSlidingExtremes() {
        this.observedLeftMostPosition = position;
        this.observedRightMostPosition = position;
    }

    /**
     * sets the sliding extremes
     */
    public void initObservedSlidingExtremes(int left, int right) {
        this.observedLeftMostPosition = left;
        this.observedRightMostPosition = right;
    }

    /**
     * sets the sliding extremes when sliding to left
     */
    public void setObservedLeftSlidingExtreme() {
        if (position < this.observedLeftMostPosition) {
            this.observedLeftMostPosition = position;
        }
    }

    /**
     * sets the sliding extremes when sliding to right
     */
    public void setObservedRightSlidingExtreme() {
        if (position > this.observedRightMostPosition) {
            this.observedRightMostPosition = position;
        }
    }

    /**
     * returns the observed slided length (with hopping)
     */
    public int getObservedSlidingLength() {
        return this.observedRightMostPosition - this.observedLeftMostPosition;
    }

    /**
     * returns current direction
     */
    public int getDirection() {
        return this.direction;
    }

    /**
     * sets the direction of the molecule
     *
     * @param direction denotes DNA strand (5'->3' (0) or 3'->5' (1))
     */
    public void setDirection(int direction) {
        this.direction = direction;
    }
    public boolean getMoveRateParam(Cell n) {
        boolean isAboutToStay = false;
        if (lastPosition > 0) {
            double lastMR =   n.dna.TFavgMoveRate[speciesID][this.lastPosition][direction];
            double currentMR = n.dna.TFavgMoveRate[speciesID][this.position][direction];
            double maxMR = n.TFspecies[this.speciesID].maxMoveRate;
            if ( (lastMR < maxMR && currentMR >maxMR) || (lastMR > maxMR && currentMR < maxMR)) {
                isAboutToStay = true;
            }
        }
        return isAboutToStay;
    }

    /**
     * returns current move rate (total rate of random walk events)
     */
    public double getMoveRate() {
        return this.moveRate;
    }

    /**
     * sets the move rate of the molecule
     */
    public void setMoveRate(Cell n) {
        boolean flag = getMoveRateParam(n);
        this.moveRate = n.TFspecies[speciesID].calcMoveRate(n.dna.TFavgMoveRate[speciesID][position][direction],flag);
        if(n.isInDebugMode() && flag){
            n.printDebugInfo(pe.time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name
            + " at position " + position + " is switching between search and recognition state and it's moverate is going to be reset to " + this.moveRate);
        }
    }

    /** FG
     * Return repressionRate, which is either rate of repression event if TF is not repressing DNA
     * or rate of 'derepression' event if TF is repressing DNA
     */
    public double getRepressionEventRate() {
        return this.repressionEventRate;
    }

    /** FG
     * Set the repressionRate variable to the proper value
     */
    public void updateRepressionRate(Cell n) {
        if (this.position == Constants.NONE) {
            repressionEventRate = 0.0;
        } else if (this.repressesDNA) {
            repressionEventRate = n.remodeller.derepressionRate / n.TFspecies[speciesID].repressionAttenuationFactor;
        } else {
            repressionEventRate = n.TFspecies[speciesID].repressionRate;
        }
    }


    /**
     * returns the last occupied position
     */
    public int getLastPosition() {
        return lastPosition;
    }

    /**
     * gets the time of last position change
     */
    public double getTimeOfLastPositionChange() {
        return this.timeOfLastPositionChange;
    }


    /**
     * updates the time of last position change
     */
    public void setTimeOfLastPositionChange(double timeOfLastPositionChange) {
        this.timeOfLastPositionChange = timeOfLastPositionChange;
    }

    /**
     * returns the time this molecule was bound
     */
    public double getTimeBound() {
        return this.timeBound;
    }

    public abstract void setCooperativityArea(Cell n, double time);

    public abstract void resetCooperativityArea(Cell n, double time);


    public abstract void setCooperativityLeft(Cell n, double time);

    public abstract void resetCooperativityLeft(Cell n, double time);

    public abstract void setCooperativityRight(Cell n, double time);

    public abstract void resetCooperativityRight(Cell n, double time);


    /**
     * this methods sets the neighbours of the current molecule once it bound to the DNA
     */
    public void setNeighboursOnBinding(Cell n) {
        //set left neighbour
        this.leftNeighbour = n.dna.getLeftNeighbour(this.position);
        if (this.leftNeighbour != Constants.NONE) {
            n.dbp[this.leftNeighbour].rightNeighbour = this.ID;
        }
        //set right neighbour
        this.rightNeighbour = n.dna.getRightNeighbour(this.position + this.size - 1);
        if (this.rightNeighbour != Constants.NONE) {
            n.dbp[this.rightNeighbour].leftNeighbour = this.ID;
        }

    }

    /**
     * resets TF neighbours when current TF unbinds
     */
    public void setNeighboursOnUnbinding(Cell n) {
        //reset left neighbour
        if (this.leftNeighbour != Constants.NONE) {
            n.dbp[this.leftNeighbour].rightNeighbour = Constants.NONE;
            this.leftNeighbour = Constants.NONE;
        }
        //reset right neighbour
        if (this.rightNeighbour != Constants.NONE) {
            n.dbp[this.rightNeighbour].leftNeighbour = Constants.NONE;
            this.rightNeighbour = Constants.NONE;

        }
    }

    /**
     * resets TF neighbours on slide left
     */
    public void setNeighboursOnSlideLeft(Cell n) {
        //set left neighbour
        this.leftNeighbour = n.dna.getLeftNeighbour(this.position);
        if (this.leftNeighbour != Constants.NONE) {
            n.dbp[this.leftNeighbour].rightNeighbour = this.ID;
        }
        //reset right neighbour
        if (this.rightNeighbour != Constants.NONE) {
            n.dbp[this.rightNeighbour].leftNeighbour = Constants.NONE;
            this.rightNeighbour = Constants.NONE;

        }
    }

    /**
     * resets TF neighbours on slide right
     */
    public void setNeighboursOnSlideRight(Cell n) {
        //reset left neighbour
        if (this.leftNeighbour != Constants.NONE) {
            n.dbp[this.leftNeighbour].rightNeighbour = Constants.NONE;
            this.leftNeighbour = Constants.NONE;
        }
        //set right neighbour
        this.rightNeighbour = n.dna.getRightNeighbour(this.position + this.size - 1);
        if (this.rightNeighbour != Constants.NONE) {
            n.dbp[this.rightNeighbour].leftNeighbour = this.ID;
        }
    }

    public abstract boolean unbindMolecule(Cell n, double time);

    public abstract int bindMolecule(Cell n, double time, int newPosition);

    public abstract int bindMolecule(Cell n, double time, int newPosition, int direction);

    public abstract int hopMolecule(Cell n, double time, int newPosition);

    public abstract int slideLeftMolecule(Cell n, double time, int newPosition, boolean isHopEvent, boolean isReflected);

    public abstract int slideRightMolecule(Cell n, double time, int newPosition, boolean isHopEvent, boolean isReflected);

    /**
     * changes the current direction (changes DNA strand)
     */
    public void changeDirection(Cell n) {
        if (n.TFreadingDirection > 1) {
            int newDirection = Utils.generateNextInteger(n.randomGenerator, 0, n.TFreadingDirection);
            if (newDirection != direction) {
                setDirection(newDirection);//draw a random direction
                setMoveRate(n); //get new move rate
            }
        }
    }

    /**
     * returns true if this molecule previously bound to the DNA
     */
    public boolean wasBound() {
        return wasBound;
    }

}
