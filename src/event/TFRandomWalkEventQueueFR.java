package event;

import java.util.Iterator;
import java.util.PriorityQueue;

import utils.Constants;
import utils.Gillespie;
import utils.Utils;

import environment.Cell;

/**
 * random walk event class using First Reaction method
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFRandomWalkEventQueueFR extends TFEventQueue {
	/**
	 * 
	 */
	private static final long serialVersionUID = 2144717026947779064L;
	public PriorityQueue<ProteinEvent> randomWalkEvents;
	
	
	/**
	 * class constructor. Initialkises the event list
	 */
	public TFRandomWalkEventQueueFR(Cell n){
		randomWalkEvents = new PriorityQueue<ProteinEvent>();	
	}
	
	/**
	 * adds a new event to the list
	 * @param pe the new protein event
	 */
	public void add(Event pe){
		randomWalkEvents.add( (ProteinEvent) pe);
	}
	
	
	/**
	 * peeks the soonest event
	 * @return
	 */
	public ProteinEvent peek(){
		return randomWalkEvents.peek();
	} 
	
	
	
	
	/**
	 * polls the soonest event
	 * @return
	 */
	public ProteinEvent pop(){
		return randomWalkEvents.poll();
	} 
	
	/**
	 * returns true if the list of events is empty or false otherwise
	 */
	public boolean isEmpty(){
		return randomWalkEvents.isEmpty();
	}
	
	/**
	 * returns the number of events in the list
	 */
	public int size(){
		return randomWalkEvents.size();
	}
	
	public ProteinEvent createNextEvent(Cell n, int moleculeID, double time) {
		if(n.dbp[moleculeID].getPosition()!=Constants.NONE) {

			int position = n.dbp[moleculeID].getPosition();
			int newPosition = position;
			int nextAction = Constants.NONE;
			int speciesID = n.dbp[moleculeID].speciesID;
			int direction = n.dbp[moleculeID].getDirection();
			boolean isHoppingEvent = false;
			//compute the time the TF stays stuck
			double propensity = n.dbp[moleculeID].getMoveRate();
			double nextTime = Gillespie.computeNextReactionTime(propensity, n.randomGenerator);
			double randomNumber = n.randomGenerator.nextDouble() * n.TFspecies[speciesID].slideRightNo;


			/*double slideLeftNo=n.TFspecies[speciesID].slideLeftNo, slideRightNo=n.TFspecies[speciesID].slideRightNo;
			if(position>0 && position<n.dna.strand.length-1 && n.TFspecies[speciesID].isBiasedRandomWalk){
				double intervalLength=n.TFspecies[speciesID].slideRightNo-n.TFspecies[speciesID].hopNo;
				double affinityRightLeftRatio = n.dna.TFavgMoveRate[speciesID][position-1][direction]/ n.dna.TFavgMoveRate[speciesID][position+1][direction];
				slideLeftNo=  intervalLength/(1+affinityRightLeftRatio);
				slideRightNo =  (affinityRightLeftRatio*intervalLength)/(1+affinityRightLeftRatio);
			}*/

			//|| n.dbp[moleculeID].moveRate < n.TFspecies[speciesID].moveRateThreshold

			if (randomNumber < n.TFspecies[speciesID].jumpNo) {
				nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
				newPosition = Constants.NONE;
			} else if (randomNumber < n.TFspecies[speciesID].hopNo) {
				isHoppingEvent = true;
				nextAction = Constants.EVENT_TF_RANDOM_WALK_HOP;
				newPosition = Utils.generateNextNormalDistributedInteger(n.randomGenerator, position, n.TFspecies[speciesID].hopSTDdisplacement);

//				if (newPosition < 0) {
//					//newPosition = 0;
//					nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
//					newPosition = Constants.NONE;
//					n.TFspecies[speciesID].countTFHopsOutside++;
//				} else if (newPosition >= n.dna.strand.length - n.dbp[moleculeID].size) {
//					//newPosition =  n.dna.strand.length-n.dbp[moleculeID].size;
//					nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
//					newPosition = Constants.NONE;
//					n.TFspecies[speciesID].countTFHopsOutside++;
//				} else if (Math.abs(newPosition - position) > n.TFspecies[speciesID].uncorrelatedDisplacementSize) {
//					//the hop size is to big and the hop becomes a jump
//					nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
//					newPosition = Constants.NONE;
//					n.TFspecies[speciesID].countTFforcedJumpsEvents++;
//				} else if (newPosition > position && n.dbp[moleculeID].size > (newPosition - position)) {
//					// the size it too small and the molecule slides to  right
//					nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_RIGHT;
//				} else if (newPosition < position && n.dbp[moleculeID].size > (position - newPosition)) {
//					// the size it too small and the molecule slides to  left
//					nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_LEFT;
//				}
			} else if (randomNumber < n.dna.TFSlideLeftNo[speciesID][position][direction]) {
				nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_LEFT;
				newPosition = position - n.TFspecies[speciesID].stepLeftSize;

			} else if (randomNumber < n.dna.TFSlideRightNo[speciesID][position][direction]) {
				nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_RIGHT;
				newPosition = position + n.TFspecies[speciesID].stepRightSize;
			}

			return new ProteinEvent(time + nextTime, moleculeID, newPosition, true, nextAction, isHoppingEvent, propensity);
		}

		return new ProteinEvent();
	}
	
	/**
	 * schedules the next non-cognate TF random walk event
	 */
	public void scheduleNextEvent(Cell n, int moleculeID, Event e) {
		ProteinEvent pe = (ProteinEvent) e;
		if (!pe.isEmpty()) {
			n.dbp[moleculeID].pe = pe;
			n.dbp[moleculeID].re = null;
			this.add(pe);
		}
	}
	
	/**
	 * updates the event of a bound molecule. Because the DM implementation always checks for the latest move rate and updates it in the list and then redraws a new event we can just schedule the event list
	 */
	public void updateNextEvent(Cell n, int moleculeID, double time){
		boolean removed = this.randomWalkEvents.remove(n.dbp[moleculeID].pe);
		if(removed){
			ProteinEvent pe = this.createNextEvent(n, moleculeID, time);
			this.scheduleNextEvent(n, moleculeID, pe);
		}
	}

	/**
	 * returns the time of the next event for a species 
	 */
	public double getNextEventTime(int moleculeID){
		Iterator<ProteinEvent> itr =randomWalkEvents.iterator();
		ProteinEvent pe=null;
		double result=0;
		
		while(itr.hasNext()) {
		    pe = itr.next();
		    if(pe.proteinID==moleculeID){
		    		result+=pe.time;
		    }
		} 
		
		return result;
	}
	
}