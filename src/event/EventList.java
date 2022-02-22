package event;

import java.io.Serializable;

import utils.Constants;
import utils.Gillespie;
import environment.Cell;

/**
 * class that contains the event list
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class EventList  implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 8623732338279526421L;
	/**
	 * 
	 */
	public TFBindingEventQueue TFBindingEventQueue;
	public TFEventQueue TFRandomWalkEventQueue;
	public TFEventQueue TFRepressionEventQueue;

	
	public EventList(Cell n){
		

		TFBindingEventQueue = new TFBindingEventQueueDNAOccupancy(n);
		
		// TF random walk event list 1D diffusion
		if(n.ip.EVENT_LIST_USES_FR.value){
			//if(n.ip.EVENT_LIST_SUBGROUP_SIZE.value >=0 && n.ip.EVENT_LIST_SUBGROUP_SIZE.value<n.dbp.length){
			//	TFRandomWalkEventQueue =  new TFRandomWalkEventQueueFRopt(n);
			//} else{
				TFRandomWalkEventQueue = new TFRandomWalkEventQueueFR(n);
				TFRepressionEventQueue = new TFRepressionEventQueueFR(n);
			//}
		} else{
			n.printDebugInfo("Warning: Direct Method is not implemented, switching to First Reaction");
			TFRandomWalkEventQueue = new TFRandomWalkEventQueueFR(n);
			TFRepressionEventQueue = new TFRepressionEventQueueFR(n);
			//TFRandomWalkEventQueue =  new TFRandomWalkEventQueueDM(n);
		}
	}
	
	
	
	/**
	 * returns the next TF binding event and removes it from the list
	 * @return
	 */
	public ProteinEvent popNextTFBindingEvent(){
		return TFBindingEventQueue.pop();
	}

	
	
	
	/**
	 * returns the next TF random walk event and removes it from the list
	 * @return
	 */
	public ProteinEvent popNextTFRandomWalkEvent(){
		return (ProteinEvent) TFRandomWalkEventQueue.pop();
	}

	/** FG
	 * returns the next TF repression event and removes it from the list
	 * @return
	 */
	public RepressionEvent popNextTFRepressionEvent(){
		return (RepressionEvent) TFRepressionEventQueue.pop();
	}
	
	/**
	 * returns a number which encodes whether the next event is TF binding ...?
	 * @return
	 */
	public int getNextEventType(){
		int result = Constants.NEXT_EVENT_IS_NONE;
		double nextEventTime =  Double.MAX_VALUE;
		
		if(!TFBindingEventQueue.isEmpty() && nextEventTime > TFBindingEventQueue.peek().time){
			nextEventTime = TFBindingEventQueue.peek().time;
			result = Constants.NEXT_EVENT_IS_TF_BINDING;
		}
		
		if(!TFRandomWalkEventQueue.isEmpty() && nextEventTime > TFRandomWalkEventQueue.peek().time){
			nextEventTime = TFRandomWalkEventQueue.peek().time;
			result = Constants.NEXT_EVENT_IS_TF_RANDOM_WALK;
		}

		if(!TFRepressionEventQueue.isEmpty() && nextEventTime > TFRepressionEventQueue.peek().time){
			nextEventTime = TFRepressionEventQueue.peek().time;
			result = Constants.NEXT_EVENT_IS_TF_REPRESSION;
		}
				
		return result;
	}
	
	/**
	 * returns the soonest event
	 * @return
	 */
	public Event getNextEvent(){
		Event e=null;
		int nextEventType = getNextEventType();

		switch(nextEventType){
			case Constants.NEXT_EVENT_IS_NONE: break; 
			case Constants.NEXT_EVENT_IS_TF_BINDING: e =this.popNextTFBindingEvent(); break;
			case Constants.NEXT_EVENT_IS_TF_RANDOM_WALK: e =this.popNextTFRandomWalkEvent(); break;
			case Constants.NEXT_EVENT_IS_TF_REPRESSION: e =this.popNextTFRepressionEvent(); break;
			default: e=null;
		}
		
		return e;
	}
	
	/**
	 * checks whether there is any event left in the entire list
	 * @return
	 */
	public boolean isEmpty(){
		if((TFBindingEventQueue==null || TFBindingEventQueue.isEmpty())
				&& (TFRandomWalkEventQueue ==null || TFRandomWalkEventQueue.isEmpty())
				&& (TFRepressionEventQueue ==null || TFRepressionEventQueue.isEmpty())){
			return true;
		}
		return false;
	}
	
	
	/**
	 * computes the total number of events in the lists
	 * @return
	 */
	public int size(){
		int result=0;


		if(TFBindingEventQueue!=null && !TFBindingEventQueue.isEmpty()){
			result++;
		}
		
		if(TFRandomWalkEventQueue !=null){
			result+= TFRandomWalkEventQueue.size();
		}
			
		return result;
	}
	
	/**
	 * schedules the next TF binding event
	 */
	public void scheduleNextTFBindingEvent( Cell n, double time){
		// if no TF then halt
		double propensity = this.TFBindingEventQueue.proteinBindingPropensitySum;
		if(n.freeTFmoleculesTotal > 0 && propensity>0){

			//generate next reaction time
			double nextTime = Gillespie.computeNextReactionTime(this.TFBindingEventQueue.proteinBindingPropensitySum, n.randomGenerator);
			
			//find next reaction (FG: randomly choose a TF specie to bind)
			int nextTFspecies = Gillespie.getNextReaction(this.TFBindingEventQueue.proteinBindingPropensitySum*n.randomGenerator.nextDouble(), this.TFBindingEventQueue.proteinBindingPropensity);
			
			if(nextTFspecies>Constants.NONE && nextTFspecies<n.TFspecies.length){
				int TFID = n.getFreeTFmolecule(nextTFspecies);
				if(TFID != Constants.NONE){
					int position =  Constants.NONE;// Gillespie.getNextReaction(n.randomGenerator.nextDouble()*n.dna.effectiveTFAffinitiesSum[nextTFspecies], n.dna.effectiveTFaffinities[nextTFspecies]); 
					this.TFBindingEventQueue.add(new ProteinEvent(time+nextTime, TFID, position, true, Constants.EVENT_TF_BINDING, false, propensity));
				}
			}
		}
	}

	/**
	 * schedules the next TF random walk or repression or unrepression event
	 */
	public void scheduleNextTFOnDNAEvent(Cell n, int moleculeID, double time) {
		double propensitySum, nextTime;
		if(!n.TFspecies[n.dbp[moleculeID].speciesID].isImmobile) {
			ProteinEvent    pe =    (ProteinEvent) TFRandomWalkEventQueue.createNextEvent(n, moleculeID, time);
			RepressionEvent re = (RepressionEvent) TFRepressionEventQueue.createNextEvent(n, moleculeID, time);
			propensitySum = Math.min(pe.propensity + re.propensity, Double.MAX_VALUE);
			if (propensitySum > 0) {
				nextTime = Gillespie.computeNextReactionTime(propensitySum, n.randomGenerator);
				if (n.randomGenerator.nextDouble() * propensitySum < pe.propensity) {
					pe.time = time + nextTime;
					TFRandomWalkEventQueue.scheduleNextEvent(n, moleculeID, pe);
				} else {
					re.time = time + nextTime;
					TFRepressionEventQueue.scheduleNextEvent(n, moleculeID, re);
				}
				//if (time >= 1.54742475190775) {
				//	n.printDebugInfo("reached bad time");
				//}
			}
		}
	}
	/*public void scheduleNextTFOnDNAEvent(Cell n, int moleculeID, double time) {
		if(!n.TFspecies[n.dbp[moleculeID].speciesID].isImmobile) {
			// FG: identify whether the next repression-type event is repression or unrepression
			double moveRate = n.dbp[moleculeID].getMoveRate();
			double reprRate = n.TFspecies[n.dbp[moleculeID].speciesID].repressionRate;
			int nextRepressionEvent = Constants.EVENT_TF_REPRESSION;
			if (((TF) n.dbp[moleculeID]).repression) {
				reprRate = n.TFspecies[n.dbp[moleculeID].speciesID].unrepressionRate;
				nextRepressionEvent = Constants.EVENT_TF_UNREPRESSION;
			}
			// FG: choose between the repression and random walk events
			double r = n.randomGenerator.nextDouble();
			if (r < moveRate / (moveRate + reprRate)) {
				this.TFRandomWalkEventQueue.scheduleNextEvent(n, moleculeID, time);
			}
			else {
				this.scheduleNextTFRepressionEvent(n, moleculeID, time, reprRate, nextRepressionEvent);
			}
		}
	}*/

	/*public void scheduleNextTFRepressionEvent(Cell n, int moleculeID, double time, double propensity, int nextEvent){
		assert n.TFspecies[n.dbp[moleculeID].speciesID].repressionRate > 0.0;
		double nextTime = Gillespie.computeNextReactionTime(propensity, n.randomGenerator);
		int position = n.dbp[moleculeID].getPosition();
		int size = n.TFspecies[n.dbp[moleculeID].speciesID].sizeTotal;
		int sizeLeft = n.TFspecies[n.dbp[moleculeID].speciesID].repressionLeftSize;
		int sizeRight = n.TFspecies[n.dbp[moleculeID].speciesID].repressionRightSize;
		this.TFRepressionEventQueue.add(new RepressionEvent(time+nextTime, nextEvent, moleculeID,
				position - sizeLeft, position + size + sizeRight));
	}*/

	/** FG
	 * schedules the next TF repression event
	 */
	public void scheduleNextTFRepressionEvent(Cell n, int moleculeID, double time) {
		RepressionEvent re = (RepressionEvent) TFRepressionEventQueue.createNextEvent(n, moleculeID, time);
		TFRepressionEventQueue.scheduleNextEvent(n, moleculeID, re);
	}


	/**
	 * schedules the next TF random walk event
	 */
	public void scheduleNextTFRandomWalkEvent(Cell n, int moleculeID, double time){
		ProteinEvent pe = (ProteinEvent) TFRandomWalkEventQueue.createNextEvent(n, moleculeID, time);
		TFRandomWalkEventQueue.scheduleNextEvent(n, moleculeID, pe);
	}

}
