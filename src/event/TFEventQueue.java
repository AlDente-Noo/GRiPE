package event;

import java.io.Serializable;

import environment.Cell;

/**
 * interface for random walk event classes
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public abstract class TFEventQueue implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -7798070719703888097L;

	public abstract void add(Event e);
	
	public abstract Event peek();
	
	public abstract Event pop();
	
	
	public abstract int  size();
	public abstract boolean isEmpty();

	public abstract Event createNextEvent(Cell n, int moleculeID, double time);
	
	public abstract void scheduleNextEvent(Cell n, int moleculeID, Event e);

	public abstract void updateNextEvent(Cell n, int moleculeID, double time);

	
}
