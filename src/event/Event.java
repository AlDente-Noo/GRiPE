package event;

import java.io.Serializable;

/**
 * Generic class for event. It contains only the event execution time
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class Event implements Comparable<Event>,  Serializable {

	private static final long serialVersionUID = 394823679537690993L;
	public double time;
	public int nextAction;

	/**
	 * class constructor
	 */
	public Event(double time, int  nextAction){
		this.time=time;
		this.nextAction = nextAction;
	}
	
	/**
	 * events are order by time at which there are scheduled 
	 */
	public int compareTo(Event e){
		if(e.time >time){
			return -1;
		} else if(e.time== time){
			return 0;
		}
		return +1;
	}

	public boolean isEmpty() {
		return time == Double.MAX_VALUE;
	}
}
