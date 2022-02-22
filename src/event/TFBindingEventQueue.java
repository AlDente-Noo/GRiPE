package event;

import java.io.Serializable;

import environment.Cell;

/**
 * this class holds all TF binding events
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFBindingEventQueue  implements Serializable {

	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 2270248942901610094L;
	protected double[] proteinBindingPropensity;// the propensity that ribosomes bind to various mRNA species
	protected double proteinBindingPropensitySum;
	private ProteinEvent bindingEvent;

	public TFBindingEventQueue(Cell n){
		//binding event
		this.bindingEvent = null;
		
		this.proteinBindingPropensitySum = 0;
		this.proteinBindingPropensity = new double[n.getNoOfDBPspecies()];
		for(int i=0;i<proteinBindingPropensity.length;i++){
			this.proteinBindingPropensity[i]=computePropensity(i,n);
			this.proteinBindingPropensitySum+=this.proteinBindingPropensity[i];
		}


	}
	
	/**
	 * returns the next binding event
	 * @return
	 */
	public ProteinEvent peek(){
		return bindingEvent;
	}
	
	/**
	 * returns the next binding event and removes it from the list
	 * @return
	 */
	public ProteinEvent pop(){
		ProteinEvent pe=bindingEvent;
		bindingEvent = null;
		return pe;
	}
	
	/**
	 * replaces the current binding event with a new one
	 * @param pe the new event
	 */
	public void add(ProteinEvent pe){
		bindingEvent=pe;
		//System.out.println(re);
	}
	
	
	/**
	 * returns true if the binding event is null
	 * @return
	 */
	public boolean isEmpty(){
		return bindingEvent==null;
	}
	
	/**
	 * computes the sum of all propensities as in Gillespie algorithm
	 * @return
	 */
	private double propensitySum(){
		double result=0;
		for(double d:proteinBindingPropensity){
			result+=d;
		}
		return result;
	}
	
	/**
	 * returns the propensity sum
	 * @return
	 */
	public double getPropensitySum(){
		return this.propensitySum();
	}
	
	/**
	 * returns the array of all propensities
	 * @return
	 */
	public double[] getPropensities(){
		return this.proteinBindingPropensity;
	}
	
	/**
	 * deletes current protein binding event
	 */
	public void clear(){
		this.bindingEvent = null;
	}


	/**
	 * updates propensities once a TF abundance has changed
	 * @param TFspeciesID the ID of the TF that changed.
	 * @param n
	 */
	public void updateProteinBindingPropensities(int TFspeciesID, Cell n){
		if(n.isInDebugMode()){
			n.printDebugInfo("Smart update of propensities for TF "+n.TFspecies[TFspeciesID].name);
		}
		//remove old propensity from propensity sum
		proteinBindingPropensitySum-=proteinBindingPropensity[TFspeciesID];
		//compute new propensity
		proteinBindingPropensity[TFspeciesID]=computePropensity(TFspeciesID,n);

		//add it to the sum
		proteinBindingPropensitySum+=proteinBindingPropensity[TFspeciesID];
		
		//when the propensities are updated the binding is deleted in order to force its regeneration in the simulator
		this.clear();
	}

	

	
	/**
	 * computes the propensity of binding a TF to the DNA
	 * @param TFspecieID the id of the species for which the propensity is computed
	 * @param n
	 * @return
	 */
	public double computePropensity(int TFspecieID, Cell n){
		return n.TFspecies[TFspecieID].assocRate * n.freeTFmolecules.get(TFspecieID).size();
	}
}
