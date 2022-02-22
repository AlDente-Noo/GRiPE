package event;

import java.io.Serializable;

import utils.Constants;
import environment.Cell;

/**
 * class that contains the list of binding events (3D diffusion) 
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFBindingEventQueueDNAOccupancy extends TFBindingEventQueue  implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 4894832130940997417L;


	public TFBindingEventQueueDNAOccupancy(Cell n) {
		super(n);
		this.updateProteinBindingPropensities(Constants.NONE, n);
	}
	
	

	/**
	 * updates propensities once a TF abundance has changed
	 * @param TFspeciesID the ID of the TF that changed.
	 * @param n
	 */
	public void updateProteinBindingPropensities(int speciesID, Cell n){
		if(n.isInDebugMode()){
			n.printDebugInfo("Full update of propensities for TF binding");
		}

		this.proteinBindingPropensitySum = 0;
		for(int i=0;i<proteinBindingPropensity.length;i++){
			this.proteinBindingPropensity[i]=computePropensity(i,n);
			this.proteinBindingPropensitySum+=this.proteinBindingPropensity[i];
		}
		
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
		//if(n.isInDebugMode()){
		//	n.printDebugInfo("propensity for species "+TFspecieID+" is "+ (n.TFspecies[TFspecieID].assocRate * n.freeTFmolecules.get(TFspecieID).size() * n.dna.effectiveTFavailabilitySum[TFspecieID]/n.dna.effectiveTFavailabilityMaxSum[TFspecieID])+"; "+n.freeTFmolecules.get(TFspecieID).size()+" free copies;  "+ n.dna.effectiveTFavailabilitySum[TFspecieID]+" positions available out of "+n.dna.effectiveTFavailabilityMaxSum[TFspecieID]);
		//}
		return n.TFspecies[TFspecieID].assocRate * n.freeTFmolecules.get(TFspecieID).size() * n.dna.effectiveTFavailabilitySum[TFspecieID]/n.dna.effectiveTFavailabilityMaxSum[TFspecieID];
	}

}
