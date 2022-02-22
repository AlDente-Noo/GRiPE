package event;

import utils.Constants;
import utils.Gillespie;
import utils.Utils;
import environment.Cell;

/**
 * random walk event class using Direct Method
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFRandomWalkEventQueueDM extends TFEventQueue {

	/**
	 * 
	 */
	private static final long serialVersionUID = -8386579868462491249L;

	private ProteinEvent nextEvent;
	
	private double[] avgMoveRate;
	private double avgMoveRateSum;

	private int sectorSize;
	private int numberOfSectors;
	private double[] avgMoveRateSectorSum;

	
	public TFRandomWalkEventQueueDM(Cell n){
		this.avgMoveRateSum=0;
		avgMoveRate = new double[n.dbp.length];
		for(int i=0;i<n.dbp.length;i++){
			avgMoveRate[i] = 0;
		}
		sectorSize = n.ip.EVENT_LIST_SUBGROUP_SIZE.value;
		if(n.ip.EVENT_LIST_SUBGROUP_SIZE.value==0){
			sectorSize =  (int) Math.floor(Math.sqrt(n.dbp.length));
		} else if(n.ip.EVENT_LIST_SUBGROUP_SIZE.value<0 || n.ip.EVENT_LIST_SUBGROUP_SIZE.value > n.dbp.length){
			sectorSize = n.dbp.length;
		}
		numberOfSectors = (int) Math.ceil((double)n.dbp.length/sectorSize);
		avgMoveRateSectorSum = new double[numberOfSectors];
	}
	
	
	/**
	 * adds an event as the next event
	 */
	public void add(Event e){
		nextEvent = (ProteinEvent) e;
	}
	
	/**
	 * peeks the next event
	 */
	public ProteinEvent peek(){
		return nextEvent;
	}
	
	/**
	 * pops the next event
	 */
	public ProteinEvent pop(){
		ProteinEvent pe = nextEvent;
		
		avgMoveRateSum = avgMoveRateSum - this.avgMoveRate[pe.proteinID];
		this.avgMoveRateSectorSum[pe.proteinID/this.sectorSize] -= this.avgMoveRate[pe.proteinID];
		this.avgMoveRate[pe.proteinID] = 0;

		nextEvent = null;
		return pe;
	}
	
	/**
	 * returns the size of this event list 0 or 1
	 */
	public int size(){
		return nextEvent==null?0:1;
	}
	
	/**
	 * checks whether there is any event as scheduled for the next event to be processed
	 */
	public boolean isEmpty(){
		return nextEvent==null;
	}

	/**
	 * updates the event of a bound molecule. Because the DM implementation always checks for the latest move rate and updates it in the list and then redraws a new event we can just schedule the event list
	 */
	public void updateNextEvent(Cell n, int moleculeID, double time){
		ProteinEvent pe = this.createNextEvent(n, moleculeID, time);
		this.scheduleNextEvent(n, moleculeID, pe);
	}

	public ProteinEvent createNextEvent(Cell n, int moleculeID, double time) {
		avgMoveRateSum= avgMoveRateSum - this.avgMoveRate[moleculeID];

		int sectorID = moleculeID/this.sectorSize;
		this.avgMoveRateSectorSum[sectorID] -= this.avgMoveRate[moleculeID];

		if(n.dbp[moleculeID].getPosition()!=Constants.NONE){
			this.avgMoveRate[moleculeID]  = n.dbp[moleculeID].getMoveRate();
			avgMoveRateSum= avgMoveRateSum + this.avgMoveRate[moleculeID];
			this.avgMoveRateSectorSum[sectorID] += this.avgMoveRate[moleculeID];

		} else{
			this.avgMoveRate[moleculeID] = 0;
		}

		if(avgMoveRateSum>0){
			//next reaction
			double nextTime = Gillespie.computeNextReactionTime(avgMoveRateSum, n.randomGenerator);
			//next molecule
			moleculeID = Gillespie.getNextReaction(n.randomGenerator.nextDouble()*avgMoveRateSum, avgMoveRate, this.avgMoveRateSectorSum, this.sectorSize);

			if(n.dbp[moleculeID].getPosition()!=Constants.NONE){
				int position = n.dbp[moleculeID].getPosition();
				int newPosition=position;
				int nextAction = Constants.NONE;
				int speciesID = n.dbp[moleculeID].speciesID;
				int direction = n.dbp[moleculeID].getDirection();

				boolean isHoppingEvent = false;

				//compute the time the TF stays stucked

				double randomNumber=n.randomGenerator.nextDouble()*n.TFspecies[speciesID].slideRightNo;
				if(randomNumber < n.TFspecies[speciesID].jumpNo){
					nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
					newPosition = Constants.NONE;
				} else if(randomNumber < n.TFspecies[speciesID].hopNo){
					isHoppingEvent = true;
					nextAction = Constants.EVENT_TF_RANDOM_WALK_HOP;
					newPosition  = Utils.generateNextNormalDistributedInteger(n.randomGenerator, position,n.TFspecies[speciesID].hopSTDdisplacement);

//					if(newPosition < 0){
//						//newPosition = 0;
//						nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
//						newPosition = Constants.NONE;
//						n.TFspecies[speciesID].countTFHopsOutside++;
//					} else if( newPosition >= n.dna.strand.length-n.dbp[moleculeID].size){
//						//newPosition =  n.dna.strand.length-n.dbp[moleculeID].size;
//						nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
//						newPosition = Constants.NONE;
//						n.TFspecies[speciesID].countTFHopsOutside++;
//					} else if (Math.abs(newPosition-position) > n.TFspecies[speciesID].uncorrelatedDisplacementSize){
//						//the hop size is to big and the hop becomes a jump
//						nextAction = Constants.EVENT_TF_RANDOM_WALK_JUMP;
//						newPosition = Constants.NONE;
//						n.TFspecies[speciesID].countTFforcedJumpsEvents++;
//					} else if(newPosition > position && n.dbp[moleculeID].size > (newPosition-position)){
//						// the size it too small and the molecule slides to  right
//						nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_RIGHT;
//					} else if(newPosition < position && n.dbp[moleculeID].size > (position - newPosition)){
//						// the size it too small and the molecule slides to  left
//						nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_LEFT;
//					}


				} else if(randomNumber < n.dna.TFSlideLeftNo[speciesID][position][direction]){
					nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_LEFT;
					newPosition= position-n.TFspecies[speciesID].stepLeftSize;
				} else if(randomNumber <  n.dna.TFSlideRightNo[speciesID][position][direction]){
					nextAction = Constants.EVENT_TF_RANDOM_WALK_SLIDE_RIGHT;
					newPosition= position+ n.TFspecies[speciesID].stepRightSize;
				}


				return new ProteinEvent(time+nextTime, moleculeID, newPosition, true, nextAction, isHoppingEvent, avgMoveRateSum);

			}
		}
		return new ProteinEvent();
	}

	/**
	 * schedules next random walk event
	 */
	public void scheduleNextEvent(Cell n, int moleculeID, Event e){
		ProteinEvent pe = (ProteinEvent) e;
		if (!pe.isEmpty()) {
			n.dbp[moleculeID].pe = pe;
			n.dbp[moleculeID].re = null;
			this.add(pe);
		}
	}
	
	
	/**
	 * returns the time of the next event for a species 
	 */
	public double getNextEventTime(int moleculeID){
		return 1.0/this.avgMoveRate[moleculeID];
	}
	
}
