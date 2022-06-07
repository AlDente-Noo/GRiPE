package agents;

import environment.Cell;
import event.ProteinEvent;
import objects.TFcooperativity;
import utils.Constants;
import utils.Utils;

import java.io.Serializable;

/**
 * class that describe a Transcription Factor which is just an instantiation of DBP
 *
 * @author n.r.zabet@gen.cam.ac.uk
 * @author fedor.garbuzov@mail.ioffe.ru
 */
public class TF extends DBP implements Serializable {

    private static final long serialVersionUID = -7084470998818754958L;
    private TFcooperativity bufferCoop;
    private int boundMolecule;

    /**
     * empty class constructor
     */
    public TF() {
        super();
        speciesID = Constants.NONE;
        direction = Constants.NONE;
        hasDNAbasedCooperativity = false;
        hasDirectCooperativity = false;
        repressesDNA = false;
    }

    /**
     * class constructor
     *
     * @param ID                       the ID of the TF
     * @param position                 the position on the DNA  (-1 free)
     * @param lastPosition             the last position of this TF
     * @param timeOfLastPositionChange the last time the position was changed
     * @param speciesID                the ID of the species to which this TF belongs
     * @param size                     TF size
     * @param direction                DNA strand
     * @param hasDNAbasedCooperativity has DNA based cooperativity
     * @param hasDirectCooperativity   has direct cooperativity
     */
    public TF(int ID, int position, int lastPosition, int timeOfLastPositionChange, int speciesID, int size,
              int direction, boolean hasDNAbasedCooperativity, boolean hasDirectCooperativity) {
        super(ID, position, lastPosition, timeOfLastPositionChange, size);
        this.speciesID = speciesID;
        this.direction = direction;
        this.hasDNAbasedCooperativity = hasDNAbasedCooperativity;
        this.hasDirectCooperativity = hasDirectCooperativity;
        this.repressesDNA = false; // FG: TODO: make it a parameter
    }

    /**
     * performs binding or random walk events
     */
    public void act(Cell n, ProteinEvent pe) {
        double previousTimeOfLastPositionChange = this.timeOfLastPositionChange;
        int oldDirection = this.direction;

        if (pe.nextAction == Constants.EVENT_TF_BINDING) {
            // binding event
            assert this.position == Constants.NONE;
            n.bindTFMolecule(n.dbp[pe.proteinID].speciesID, pe.time);
        } else {
            // random walk event
            if (pe.nextAction == Constants.EVENT_TF_RANDOM_WALK_JUMP) {
                //the TF unbinds
                unbindMolecule(n, pe.time);
            } else if (pe.nextAction == Constants.EVENT_TF_RANDOM_WALK_HOP) {
				assert pe.isHoppingEvent; // FG: seems like if the event is hop then isHoppingEvent is always true
                // FG: prevent hopping of the repressed molecule
                if (this.isRepressingDNA() || this.isRepressed(n)) {
                    if (n.isInDebugMode()) {
                        n.printDebugInfo(pe.time + ": attempted to hop TF " + this.ID + " of type " +
								n.TFspecies[speciesID].name + " from position " + position + " to " + pe.position +
								", but the TF is repressing or repressed");
                    }
                } else {
                    if (pe.position < 0 || pe.position > n.dna.strand.length - this.size) {
						// hop outside of DNA results in unbinding
                        unbindMolecule(n, pe.time);
                        n.TFspecies[speciesID].countTFHopsOutside++;
                    } else if (Math.abs(pe.position - this.position) > n.TFspecies[speciesID].uncorrelatedDisplacementSize) {
                    	// too big hop results in unbinding
                        unbindMolecule(n, pe.time);
                        n.TFspecies[speciesID].countTFforcedJumpsEvents++;
                    } else if (pe.position == this.position) {
						//the TF hops and rebinds at the same position
                        hopMoleculeSamePos(n, pe.time, this.position);
                    } else {
                    	//the TF hops and rebinds at a different position
                        hopMolecule(n, pe.time, pe.position);
                    }
                }
            } else if (pe.nextAction == Constants.EVENT_TF_RANDOM_WALK_SLIDE_LEFT) {
                // FG: prevent sliding of the repressed molecule
                if (this.isRepressingDNA() || this.isRepressed(n)) {
                    if (n.isInDebugMode()) {
                        n.printDebugInfo(pe.time + ": attempted to slide left TF " + this.ID + " of type " +
								n.TFspecies[speciesID].name + " from position " + position + " to " + pe.position +
								", but the TF is repressing or repressed");
                    }
                } else {
                    if (n.dna.isAbsorbing && pe.position < 0) {
                    	// absorbing boundary conditions cause unbinding when DNA end is reached
						if (n.isInDebugMode()) {
							n.printDebugInfo(pe.time + " TF " + this.ID + " reached left limit and, due to absorbing " +
									"boundary condition, it will unbind!");
						}
                        unbindMolecule(n, pe.time);
                    } else if (n.dna.isPeriodic && pe.position < 0) {
                    	// periodic b.c. cause the movement of the molecule to the opposite end of DNA when a DNA end is reached
                        if (n.isInDebugMode()) {
                            n.printDebugInfo(pe.time + " TF " + this.ID + " reached left limit and, due to periodic " +
                                    "boundary condition, it will attempt to bind at the end!");
                        }
                        hopMolecule(n, pe.time, n.dna.strand.length - this.size - 1);
                    } else {
                        if (pe.position >= 0) {
                        	// slide molecule is possible
                            slideLeftMolecule(n, pe.time, pe.position, pe.isHoppingEvent);
                        } else {
                        	// reflective b.c. prevent movement of the molecule outside of DNA
                            if (n.isInDebugMode()) {
                                n.printDebugInfo(pe.time + ": attempted to slide left TF " + this.ID +
										" of type " + n.TFspecies[speciesID].name + " from position " + position +
										", but the left border is reached");
                            }
                            resetPosition(pe.time);
                        }
                        if (pe.isHoppingEvent) {
                        	// sliding event can be invoked by hopping event with small distance (less than slide length)
							// in that case TF unbinds from DNA and can rebind to another strand
                            this.changeDirection(n);
                        }
                    }
                }
            } else if (pe.nextAction == Constants.EVENT_TF_RANDOM_WALK_SLIDE_RIGHT) {
                // similar to the slide left
                assert(this.position != Constants.NONE);
                if (this.isRepressingDNA() || this.isRepressed(n)) {
                    if (n.isInDebugMode()) {
                        n.printDebugInfo(pe.time + ": attempted to slide right TF " + this.ID + " of type " +
								n.TFspecies[speciesID].name + " from position " + position + " to " + pe.position +
								", but the TF is repressing or repressed");
                    }
                } else {
                    if (n.dna.isAbsorbing && pe.position >= (n.dna.strand.length - this.size)) {
                        if (n.isInDebugMode()) {
                            n.printDebugInfo(pe.time + " TF " + this.ID + " reached right limit and, due to absorbing" +
                                    " boundary condition, it will unbind!");
                        }
                        unbindMolecule(n, pe.time);
                    } else if (n.dna.isPeriodic && pe.position >= (n.dna.strand.length - this.size)) {
                        if (n.isInDebugMode()) {
                            n.printDebugInfo(pe.time + " TF " + this.ID + " reached right limit and, due to periodic " +
                                    "boundary condition, it will attempt to bind at the start!");
                        }
                        hopMolecule(n, pe.time, 0);
                    } else {
                        if (pe.position < (n.dna.strand.length - this.size)) {
                            slideRightMolecule(n, pe.time, pe.position, pe.isHoppingEvent);
                        } else {
                            if (n.isInDebugMode()) {
                                n.printDebugInfo(pe.time + ": attempted to slide right TF " + this.ID +
										" of type " + n.TFspecies[speciesID].name + " from position " + position +
										", but the right border is reached");
                            }
                            resetPosition(pe.time);
                        }
                        if (pe.isHoppingEvent) {
                            this.changeDirection(n);
                        }
                    }
                }
            }
        }

        double timeBound = 0;
        // update bound time
        if (this.lastPosition != Constants.NONE) {
            timeBound = this.timeOfLastPositionChange - previousTimeOfLastPositionChange;
            updateBoundTime(n, timeBound, oldDirection, lastPosition);
        }

        // update target sited statistics
		// if random walk or binding or unbinding was performed (which could be prevented by repression)
        if (!this.isRepressed(n) && (!this.isRepressingDNA() || this.position == Constants.NONE)) {
            if (this.position != Constants.NONE && n.dna.isTargetSite[this.speciesID][this.position][this.direction] != Constants.NONE) {
                n.updateTargetSiteStatistics(n.dna.isTargetSite[this.speciesID][this.position][this.direction],
                        this.timeOfLastPositionChange, true, timeBound);
                if (n.runUntilTSReached) {
                    n.areTargetSitesToBeReached = n.tsg.areTargetSitesToBeReached();
                }
            } else if (this.lastPosition != Constants.NONE && n.dna.isTargetSite[this.speciesID][this.lastPosition][oldDirection] != Constants.NONE) {
                n.updateTargetSiteStatistics(n.dna.isTargetSite[this.speciesID][this.lastPosition][oldDirection],
                        this.timeOfLastPositionChange, false, timeBound);
            }
        }
    }

    /**
     * generates a string describing current state
     */
    public String toString() {
        return "TF " + ID + " of type " + speciesID + " at position " + position + " (" + directionToString() + ")" +
				" previously at " + lastPosition + " has the time of last  change " + timeOfLastPositionChange +
				" and size " + size;
    }

    /**
     * returns a string describing the current direction
     */
    public String directionToString() {
        String dir = "none";
        if (this.direction == 0) {
            dir = "5'->3' ";
        }
        if (this.direction == 1) {
            dir = "3'->5' ";
        }

        return dir;
    }

    /**
     * returns a String describing obstacle which prevented movement of TF
     */
    public String moleculeBlockingString(Cell n, int blockingMolecule) {
        String str = "";
        if (blockingMolecule == Constants.NONE) {
            str += "outside of the strand or repressed or closed";
        } else {
            str += "blocked by molecule " + blockingMolecule;
            str += " (TF " + n.TFspecies[n.dbp[blockingMolecule].speciesID].name + ")";
        }

        return str;
    }

    /**
     * updates dna occupancy
     */
    public void updateBoundTime(Cell n, double timeBound, int direction, int position) {
        n.dna.effectiveTFOccupancy[speciesID][position][direction] += timeBound;
    }


    /**
     * sets the position of the TF
     */
    public void setPosition(int newPosition, double timeOfLastPositionChange) {
        this.lastPosition = this.position;
        this.position = newPosition;
        this.timeOfLastPositionChange = timeOfLastPositionChange;
    }

    /**
     * sets the position of the TF as free
     */
    public void clearPosition(double timeOfLastPositionChange) {
        this.lastPosition = this.position;
        this.position = Constants.NONE;
        this.timeOfLastPositionChange = timeOfLastPositionChange;
    }

    /**
     * resets position to current position
     */
    public void resetPosition(double time) {
        this.lastPosition = this.position;
        this.timeOfLastPositionChange = time;
    }

    public void setMoveRate(Cell n) {
        moveRate = n.dna.TFavgMoveRate[speciesID][position][direction];
    }

    /**
     * once a TF moves if it releases a cooperativity site then update the move rate
     *
     * @param n
     */
    public void setCooperativityArea(Cell n, double time) {
        if (n.TFspecies[this.speciesID].isCooperativeSite(this.position, this.direction)) {
            bufferCoop = n.TFspecies[this.speciesID].getCooperativity(this.position, this.direction);
            int startDir = 0, endDir = n.dna.TFdirections;
            if (bufferCoop.direction1 != Constants.NONE) {
                startDir = bufferCoop.direction1;
                endDir = bufferCoop.direction1 + 1;
            }

            for (int i = (int) bufferCoop.region1.start; i < bufferCoop.region1.end; i++) {
                for (int j = startDir; j < endDir; j++) {
                    n.dna.TFavgMoveRate[bufferCoop.species1ID][i][j] /= bufferCoop.affinityIncrease;

                    if (bufferCoop.isReversible) {
                        boundMolecule = n.dna.getBoundMolecule(i);
                        if (boundMolecule != Constants.NONE && n.dbp[boundMolecule].speciesID == bufferCoop.species1ID && n.dbp[boundMolecule].direction == j) {
                            n.dbp[boundMolecule].setMoveRate(n);
                            n.eventQueue.TFRandomWalkEventQueue.updateNextEvent(n, boundMolecule, time);
                        }
                    }
                }
            }
            if (n.isInDebugMode()) {
                n.printDebugInfo("The affinity of the TF " + n.TFspecies[bufferCoop.species1ID].name + " for site " + bufferCoop.region1 + " was increased by " + bufferCoop.affinityIncrease + " due to the binding of TF " + n.TFspecies[bufferCoop.species0ID].name + " to site" + bufferCoop.region0);
            }

        }
    }


    /**
     * once a TF moves if it releases a cooperativity site then update the move rate
     *
     * @param n
     */
    public void resetCooperativityArea(Cell n, double time) {
        if (n.TFspecies[this.speciesID].isCooperativeSite(this.position, this.direction)) {
            bufferCoop = n.TFspecies[this.speciesID].getCooperativity(this.position, this.direction);
            int startDir = 0, endDir = n.dna.TFdirections;
            if (bufferCoop.direction1 != Constants.NONE) {
                startDir = bufferCoop.direction1;
                endDir = bufferCoop.direction1 + 1;
            }

            //double buffer;
            for (int i = (int) bufferCoop.region1.start; i < bufferCoop.region1.end; i++) {
                for (int j = startDir; j < endDir; j++) {
                    n.dna.TFavgMoveRate[bufferCoop.species1ID][i][j] *= bufferCoop.affinityIncrease;

                    if (bufferCoop.isReversible) {
                        boundMolecule = n.dna.getBoundMolecule(i);
                        if (boundMolecule != Constants.NONE && n.dbp[boundMolecule].speciesID == bufferCoop.species1ID && n.dbp[boundMolecule].direction == j) {
                            n.dbp[boundMolecule].setMoveRate(n);
                            n.eventQueue.TFRandomWalkEventQueue.updateNextEvent(n, boundMolecule, time);
                        }
                    }
                }
            }

            if (n.isInDebugMode()) {
                n.printDebugInfo("The affinity of the TF " + n.TFspecies[bufferCoop.species1ID].name + " for site " + bufferCoop.region1 + " was reseted due to the unbinding of TF " + n.TFspecies[bufferCoop.species0ID].name + " from site" + bufferCoop.region0);
            }
        }
    }

    /**
     * once a TF moves if it releases a cooperativity site then update the move rate
     *
     * @param n
     */
    public void setCooperativityRight(Cell n, double time) {
        if (this.rightNeighbour != Constants.NONE && this.stickToRight == Constants.NONE) {
            TFcooperativity coop = n.TFspecies[this.speciesID].getDirectCooperativityRight(this.rightNeighbour,
                    this.direction, n.dbp[this.rightNeighbour].getDirection());

            if (coop != null && n.randomGenerator.nextDouble() < coop.dimerisationProbability) {

                this.stickToRight = this.rightNeighbour;
                n.dbp[this.rightNeighbour].stickToLeft = this.ID;

                //move rate
                double buffer = this.moveRate;
                this.moveRate = coop.computeMoveRate(this.moveRate, n.dbp[this.rightNeighbour].moveRate,
                        coop.affinityIncrease);
                n.dbp[this.rightNeighbour].moveRate = coop.computeMoveRate(n.dbp[this.rightNeighbour].moveRate,
                        buffer, coop.affinityIncrease);
                n.eventQueue.TFRandomWalkEventQueue.updateNextEvent(n, this.rightNeighbour, time);

                if (n.isInDebugMode()) {
                    n.printDebugInfo(time + ": TF molecule " + this.ID + " become cooperative to TF molecule " + this.rightNeighbour);
                }
            }

        }
    }


    /**
     * reset move rates and
     */
    public void resetCooperativityRight(Cell n, double time) {
        if (this.stickToRight != Constants.NONE) {
            this.setMoveRate(n);
            n.dbp[this.stickToRight].setMoveRate(n);
            //reshedule the right neighbour move rate
            n.eventQueue.TFRandomWalkEventQueue.updateNextEvent(n, this.stickToRight, time);

            n.dbp[this.stickToRight].stickToLeft = Constants.NONE;
            this.stickToRight = Constants.NONE;
        }
    }

    /**
     * once a TF moves if it releases a cooperativity site then update the move rate
     *
     * @param n
     */
    public void setCooperativityLeft(Cell n, double time) {
        if (this.leftNeighbour != Constants.NONE && this.stickToLeft == Constants.NONE) {
            TFcooperativity coop = n.TFspecies[this.speciesID].getDirectCooperativityLeft(this.leftNeighbour,
                    this.direction, n.dbp[this.leftNeighbour].getDirection());

            if (coop != null && n.randomGenerator.nextDouble() < coop.dimerisationProbability) {
                this.stickToLeft = this.leftNeighbour;
                n.dbp[this.leftNeighbour].stickToRight = this.ID;

                //move rates
                double buffer = this.moveRate;
                this.moveRate = coop.computeMoveRate(this.moveRate, n.dbp[this.leftNeighbour].moveRate,
                        coop.affinityIncrease);
                n.dbp[this.leftNeighbour].moveRate = coop.computeMoveRate(n.dbp[this.leftNeighbour].moveRate, buffer,
                        coop.affinityIncrease);
                n.eventQueue.TFRandomWalkEventQueue.updateNextEvent(n, this.leftNeighbour, time);

                if (n.isInDebugMode()) {
                    n.printDebugInfo(time + ": TF molecule " + this.ID + " become cooperative to TF molecule " + this.leftNeighbour);
                }
            }
        }
    }

    /**
     * reset move rates and
     */
    public void resetCooperativityLeft(Cell n, double time) {
        if (this.stickToLeft != Constants.NONE) {
            n.dbp[this.stickToLeft].setMoveRate(n);
            //reshedule the right neighbour move rate
            n.eventQueue.TFRandomWalkEventQueue.updateNextEvent(n, this.stickToLeft, time);
            n.dbp[this.stickToLeft].stickToRight = Constants.NONE;
            this.stickToLeft = Constants.NONE;
        }
    }


    /**
     * attempts to unbind the current molecule from the DNA
     *
     * @param time the time
     */
    public boolean unbindMolecule(Cell n, double time) {

        boolean unbound = false;
        int pos = this.position;

        // FG: schedule unrepression event if the unbound molecule is repressing the DNA region
        if (n.dbp[this.ID].repressesDNA) {
            n.dbp[this.ID].repressionRate = n.remodeller.unrepressionRate;
            n.eventQueue.scheduleNextTFRepressionEvent(n, this.ID, time);
            n.dbp[this.ID].repressesDNA = false;
        }

        unbound = n.dna.unbindMolecule(n, this.ID, this.position, this.size);
        assert unbound; // FG: seems like unbinding cannot fail

        //if in debug mode print status
        if (n.isInDebugMode()) {
            if (!unbound) {
                n.printDebugInfo(time + ": attempted to unbind TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " which was bound to the DNA at position " + this.position);
            } else {
                n.printDebugInfo(time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " unbound " +
                        "from the DNA at position " + this.position);
            }
        }

        if (unbound) {
            //count event type
            n.TFspecies[speciesID].countTFUnbindingEvents++;

            //reset cooperativity before deleting position
            if (this.hasDNAbasedCooperativity) {
                this.resetCooperativityArea(n, time);
            }
            if (this.hasDirectCooperativity) {
                this.resetCooperativityLeft(n, time);
                this.resetCooperativityRight(n, time);
            }

            //clear position
            this.clearPosition(time);

            //reset neighbours
            setNeighboursOnUnbinding(n);

            //record sliding length if necessary
            if (n.ip.OUTPUT_SLIDING_LENGTHS.value) {
                n.TFspecies[speciesID].slidingLength.add(getSlidingLength());
                n.TFspecies[speciesID].slidingEvents.add(getSlidingEvents());
                n.TFspecies[speciesID].observedSlidingLength.add(getObservedSlidingLength());
            }

            //update the list of free TFs and the binding event
            n.freeTFmoleculesTotal++;
            n.freeTFmolecules.get(speciesID).add(this.ID);
            n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(speciesID, n);
            // FG: update repression rate
            this.updateRepressionRate(n);
        }
        return unbound;

    }


    /**
     * attempts to slide to right on the DNA a TF molecule
     *
     * @param moleculeID the ID of the TF that slides left
     * @param time       the time
     */
    public int slideRightMolecule(Cell n, double time, int newPosition, boolean isHopEvent) {

        int bound = Constants.NONE;

        //if(n.isInDebugMode()){
        //	n.printDebugInfo("attempted to slide right where there is a molecule "+this.leftNeighbour);
        //}

        //attempt to slide right the molecule on the DNA
        if (this.rightNeighbour != Constants.NONE) {
            bound = this.rightNeighbour;
        } else {
            bound = n.dna.slideRight(n, this.ID, this.position, this.size, newPosition - this.position,
                    n.ip.CHECK_OCCUPANCY_ON_SLIDING.value);
        }

        //print debug info
        if (n.isInDebugMode()) {
            String slideTxt = "slide";
            if (isHopEvent) {
                slideTxt = "hop";
            }
            if (bound != this.ID) {
                n.printDebugInfo(time + ": attempted to " + slideTxt + " right TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " from position " + position + " by " + (newPosition - position) + " bp, but this is " + moleculeBlockingString(n, bound));
            } else {
                n.printDebugInfo(time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " " + slideTxt + " right from position " + position + " by " + (newPosition - position) + " bp");
            }
        }

        if (bound != this.ID) {
            if (bound != Constants.NONE) {
                n.dna.collisionsCount[newPosition]++;
            }
            bound = Constants.NONE;

            //couldn't slide but there is a chance it will release due to the collision
            if (n.TFspecies[speciesID].collisionUnbindingProbability > 0.0 && n.randomGenerator.nextDouble() < n.TFspecies[speciesID].collisionUnbindingProbability) {
                this.unbindMolecule(n, time);
            } else {
                this.resetPosition(time);
                if (isHopEvent) {
                    n.TFspecies[speciesID].countTFHoppingEvents++;
                }
            }

        } else {
            //reset cooperativity if case
            if (hasDNAbasedCooperativity) {
                resetCooperativityArea(n, time);
            }
            //if slide right the left neighbour molecule will lose cooperativity
            if (hasDirectCooperativity) {
                resetCooperativityLeft(n, time);
            }

            setPosition(newPosition, time);
            setMoveRate(n); //get new move rate
            setNeighboursOnSlideRight(n);

            //set cooperativity on slide right
            if (hasDNAbasedCooperativity) {
                setCooperativityArea(n, time);
            }
            //if slide right the right neighbour might gain cooperativity cooperativity
            if (hasDirectCooperativity) {
                setCooperativityRight(n, time);
            }

            //sliding statistics
            if (n.ip.OUTPUT_SLIDING_LENGTHS.value) {
                setObservedRightSlidingExtreme();
                //is a hop event
                if (isHopEvent) {
                    n.TFspecies[speciesID].slidingLength.add(getSlidingLength());
                    n.TFspecies[speciesID].slidingEvents.add(getSlidingEvents());
                    initSlidingExtremes();
                } else {
                    setRightSlidingExtreme();
                }
            }

            if (n.ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value) {
                n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(speciesID, n);
            }

            if (isHopEvent) {
                n.TFspecies[speciesID].countTFHoppingEvents++;
            } else {
                n.TFspecies[speciesID].countTFSlideRightEvents++;
            }


        }
        return bound;

    }


    /**
     * attempts to slide to left on the DNA a TF molecule
     *
     * @param moleculeID the ID of the TF that slides left
     * @param time       the time
     */
    public int slideLeftMolecule(Cell n, double time, int newPosition, boolean isHopEvent) {
        int bound = Constants.NONE;

        //if(n.isInDebugMode()){
        //	n.printDebugInfo("attempted to slide left where there is a molecule "+this.leftNeighbour);
        //}

        if (this.leftNeighbour != Constants.NONE) {
            bound = this.leftNeighbour;
        } else {
            bound = n.dna.slideLeft(n, this.ID, position, size, position - newPosition,
                    n.ip.CHECK_OCCUPANCY_ON_SLIDING.value);
        }

        //print debug info
        if (n.isInDebugMode()) {
            String slideTxt = "slide";
            if (isHopEvent) {
                slideTxt = "hop";
            }
            if (bound != this.ID) {
                n.printDebugInfo(time + ": attempted to " + slideTxt + " left TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " from position " + position + " by " + (position - newPosition) + " bp, but this is " + moleculeBlockingString(n, bound));
            } else {
                n.printDebugInfo(time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " " + slideTxt + " left from position " + position + " by " + (position - newPosition) + " bp");
            }
        }


        if (bound != this.ID) {
            if (bound != Constants.NONE) {
                n.dna.collisionsCount[newPosition]++;
            }
            bound = Constants.NONE;

            //couldn't slide but there is a chance it will release due to the collision
            if (n.TFspecies[speciesID].collisionUnbindingProbability > 0.0 && n.randomGenerator.nextDouble() < n.TFspecies[speciesID].collisionUnbindingProbability) {
                this.unbindMolecule(n, time);
            } else {
                this.resetPosition(time);
                if (isHopEvent) {
                    n.TFspecies[speciesID].countTFHoppingEvents++;
                }
            }

        } else {
            //reset cooperativity if case
            if (hasDNAbasedCooperativity) {
                resetCooperativityArea(n, time);
            }

            //if slide left the right neighbour molecule will lose cooperativity
            if (hasDirectCooperativity) {
                resetCooperativityRight(n, time);
            }

            setPosition(newPosition, time);
            setMoveRate(n); //get new move rate
            setNeighboursOnSlideLeft(n);

            //set cooperativity if the case
            if (hasDNAbasedCooperativity) {
                setCooperativityArea(n, time);
            }
            //if slide right the right neighbour might gain cooperativity
            if (hasDirectCooperativity) {
                setCooperativityLeft(n, time);
            }

            //slidingstatistics
            if (n.ip.OUTPUT_SLIDING_LENGTHS.value) {
                setObservedLeftSlidingExtreme();
                //is a hop event
                if (isHopEvent) {
                    n.TFspecies[speciesID].slidingLength.add(getSlidingLength());
                    n.TFspecies[speciesID].slidingEvents.add(getSlidingEvents());
                    initSlidingExtremes();
                } else {
                    setLeftSlidingExtreme();
                }
            }

            if (n.ip.SLIDING_AND_HOPPING_AFFECTS_TF_ASSOC_RATE.value) {
                n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(speciesID, n);
            }
            if (isHopEvent) {
                n.TFspecies[speciesID].countTFHoppingEvents++;
            } else {
                n.TFspecies[speciesID].countTFSlideLeftEvents++;
            }
        }


        return bound;

    }

    public void hopMoleculeSamePos(Cell n, double time, int position) {
        this.resetPosition(time);
        this.changeDirection(n);
        //hop to the same position
        n.TFspecies[speciesID].countTFHoppingEvents++;
        //set the sliding lengths accordingly
        if (n.ip.OUTPUT_SLIDING_LENGTHS.value) {
            n.TFspecies[speciesID].slidingLength.add(getSlidingLength());
            n.TFspecies[speciesID].slidingEvents.add(getSlidingEvents());
        }
        this.initSlidingExtremes();
        if (n.isInDebugMode()) {
            n.printDebugInfo(time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " hop at the same" +
                    " position " + position);
        }
    }

    /**
     * attempts to rebind a TF to the DNA
     *
     * @param moleculeID the ID of the TF that binds
     * @param time       the time
     */
    public int hopMolecule(Cell n, double time, int newPosition) {
        int bound = Constants.NONE, leftExtreme = Constants.NONE, rightExtreme = Constants.NONE;

        // use slide functions for hops by the distance less than molecule size, otherwise canBind will be false
        if (position < newPosition && newPosition < position + size) {
            return slideRightMolecule(n, time, newPosition, true);
        }
        if (position - size < newPosition && newPosition < position) {
            return slideLeftMolecule(n, time, newPosition, true);
        }

        //check if the TF can rebind
        boolean canBind = n.dna.effectiveTFavailability[this.speciesID][newPosition];

        // if it can rebind or it doesn't matter if it can rebind => unbinds the TF
        if (canBind || !n.TFspecies[speciesID].stallsHoppingIfBlocked) {
            if (n.ip.OUTPUT_SLIDING_LENGTHS.value) {
                leftExtreme = this.observedLeftMostPosition;
                rightExtreme = this.observedRightMostPosition;

            }

            //keep old last position
            int oldLastPosition = this.lastPosition;

            unbindMolecule(n, time);

            //if unbound but can rebind then rebind the TF
            if (canBind) {
                //attempt rebind
                bound = bindMolecule(n, time, newPosition);
                if (n.ip.OUTPUT_SLIDING_LENGTHS.value && leftExtreme != Constants.NONE && rightExtreme != Constants.NONE) {
                    //if this was a hopping then the observed sliding length that was recorded is removed and the
                    // current one is properly initialised
                    this.initObservedSlidingExtremes(leftExtreme, rightExtreme);
                    n.TFspecies[speciesID].observedSlidingLength.remove(n.TFspecies[speciesID].observedSlidingLength.size() - 1);
                }

                //if rebound the last position is reseted to the previous one before unbinding
                this.lastPosition = oldLastPosition;

                n.TFspecies[speciesID].countTFHoppingEvents++;

            }

        }

        if (n.isInDebugMode()) {
            if (bound == ID) {
                n.printDebugInfo(time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " hoped at " +
                        "position " + this.position + " in direction " + this.directionToString());
            } else {
                //int boundProtein = n.dna.getBoundProtein(newPosition,size);
                String str = time + ": failed to hop TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " " +
                        "to position " + newPosition + " which is occupied or repressed, ";
                if (position == Constants.NONE) {
                    str += "molecule released";
                } else {
                    str += "molecule stalled";
                }
                n.printDebugInfo(str);
            }
        }


        return bound;

    }


    /**
     * attempts to bind a TF to the DNA
     *
     * @param moleculeID the ID of the species of the molecule that binds
     * @param time       the time
     */
    public int bindMolecule(Cell n, double time, int newPosition) {
        int bound = Constants.NONE;

        bound = n.dna.bindMolecule(n, ID, newPosition, size, n.ip.CHECK_OCCUPANCY_ON_BINDING.value);
        if (bound != ID) {
            bound = Constants.NONE;
        } else {
            setPosition(newPosition, time);
            setDirection(Utils.generateNextInteger(n.randomGenerator, 0, n.TFreadingDirection));//draw a random
            // direction
            setMoveRate(n); //get new move rate
            this.wasBound = true;

            //set cooperativy after binding
            if (hasDNAbasedCooperativity) {
                setCooperativityArea(n, time);
            }

            if (hasDirectCooperativity) {
                setCooperativityLeft(n, time);
                setCooperativityRight(n, time);
            }

            //init sliding length measurements
            if (n.ip.OUTPUT_SLIDING_LENGTHS.value) {
                initSlidingExtremes();
                initObservedSlidingExtremes();
            }

            //set neighbours
            setNeighboursOnBinding(n);

            //cound binding events
            n.TFspecies[speciesID].countTFBindingEvents++;


            //update the list of free TFs
            n.freeTFmoleculesTotal--;
            n.freeTFmolecules.get(speciesID).remove(n.freeTFmolecules.get(speciesID).size() - 1);

            //update propensities
            n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(speciesID, n);
            // FG: update repression rate
            this.updateRepressionRate(n);
        }


        //if in debug mode print status
        if (n.isInDebugMode()) {
            if (bound != ID) {
                n.printDebugInfo(time + ": failed attempted to bind TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " at position " + this.position + " in direction " + this.directionToString());
            } else {
                n.printDebugInfo(time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " bound at " +
                        "position " + this.position + " in direction " + this.directionToString());
            }
        }

        return bound;

    }


    /**
     * attempts to bind a TF to the DNA
     *
     * @param moleculeID the ID of the species of the molecule that binds
     * @param time       the time
     */
    public int bindMolecule(Cell n, double time, int newPosition, int direction) {
        int bound = Constants.NONE;

        bound = n.dna.bindMolecule(n, ID, newPosition, size, n.ip.CHECK_OCCUPANCY_ON_BINDING.value);
        if (bound != ID) {
            bound = Constants.NONE;
        } else {
            setPosition(newPosition, time);
            setDirection(direction);//draw a random direction
            setMoveRate(n); //get new move rate
            this.wasBound = true;

            //set cooperativy after binding
            if (hasDNAbasedCooperativity) {
                setCooperativityArea(n, time);
            }

            if (hasDirectCooperativity) {
                setCooperativityLeft(n, time);
                setCooperativityRight(n, time);
            }

            //init sliding length measurements
            if (n.ip.OUTPUT_SLIDING_LENGTHS.value) {
                initSlidingExtremes();
                initObservedSlidingExtremes();
            }

            //set neighbours
            setNeighboursOnBinding(n);

            //cound binding events
            n.TFspecies[speciesID].countTFBindingEvents++;


            //update the list of free TFs
            n.freeTFmoleculesTotal--;
            n.freeTFmolecules.get(speciesID).remove(n.freeTFmolecules.get(speciesID).size() - 1);

            //update propensities
            n.eventQueue.TFBindingEventQueue.updateProteinBindingPropensities(speciesID, n);
            // FG: update repression rate and repression propensity
            this.updateRepressionRate(n);
        }


        //if in debug mode print status
        if (n.isInDebugMode()) {
            if (bound != ID) {
                n.printDebugInfo(time + ": failed attempted to bind TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " at position " + this.position + " in direction " + this.directionToString());
            } else {
                n.printDebugInfo(time + ": TF " + this.ID + " of type " + n.TFspecies[speciesID].name + " bound at " +
						"position " + this.position + " in direction " + this.directionToString());
            }
        }

        return bound;

    }

}
