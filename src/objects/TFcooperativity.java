package objects;


import java.io.Serializable;

import utils.Constants;

/**
 * class that describes cooperativity
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFcooperativity  implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 4619689557618327801L;
	public int species0ID;
	public int species1ID;
	public int type;
	public int direction0; //the orientation of the original TF species
	public int direction1; // the orientation of the new species specified by the speciesID
	public DNAregion region0;
	public DNAregion region1;
	public double dimerisationProbability;
	public double affinityIncrease;
	public boolean isReversible;
	
	public boolean isAdditive;
	public boolean isFixed;
	
	public TFcooperativity(int species0ID, int species1ID, int type, int direction0, int direction1, DNAregion region0, DNAregion region1, double dimerisationProbability, double affinityIncrease, boolean isReversible, boolean isAdditive, boolean isFixed){
		this.species0ID = species0ID;
		this.species1ID = species1ID;
		this.type = type;
		this.direction0 = direction0; 
		this.direction1 = direction1;
		this.region0 = region0;
		this.region1 = region1;
	
		this.dimerisationProbability = dimerisationProbability;
		this.affinityIncrease = affinityIncrease;
		this.isReversible = isReversible;
		this.isAdditive = isAdditive;
		this.isFixed = isFixed;
	}	

	
	/**
	 * swaps data between species 0 and species 1 
	 */
	public void swap(){
		
		int buffer;
		
		buffer = species0ID;
		this.species0ID = species1ID;
		species1ID = buffer;
		
		buffer = direction0; 
		this.direction0 = direction1;
		direction1 = buffer;
		
		DNAregion bufferRegion;
		bufferRegion = new DNAregion(region0);
		region0.loadRegion(region1);
		region1.loadRegion(bufferRegion);
	}
	
	/**
	 * constructor that copies another object
	 * @param coop
	 */
	public TFcooperativity( TFcooperativity coop){
		this.species0ID = coop.species0ID;
		this.species1ID = coop.species1ID;
		this.type = coop.type;
		this.direction0 = coop.direction0; 
		this.direction1 = coop.direction1;
		this.region0 = new DNAregion(coop.region0);
		this.region1 = new DNAregion(coop.region1);
	
		this.dimerisationProbability = coop.dimerisationProbability;
		this.affinityIncrease = coop.affinityIncrease;
		this.isAdditive = coop.isAdditive;
		this.isFixed = coop.isFixed;
	}	

	/**
	 * computes the move rate after is affected by cooperativity
	 * @param moveRate0 the move rate of the current 
	 * @param moveRate1
	 * @param factor
	 * @return
	 */
	public double computeMoveRate(double moveRate0, double moveRate1, double factor){
		return moveRate0*factor;
	}
		
	/**
	 * generates a string that describes the coop
	 */
	public String toString(){
		String str="";
		str+="["+this.species0ID;
		if(this.direction0!=Constants.NONE){
			str+="("+this.direction0+")";
		}
		str+="->"+this.species1ID;
		if(this.direction1!=Constants.NONE){
			str+="("+this.direction1+")";
		}
		str+="]";
		
		if(this.type == Constants.TF_COOPERATIVITY_TYPE_DNA){
			str+="[type=DNA]["+this.region0+"->"+this.region1+"]";
			if(this.isReversible){
				str+="[reversible]";
			} else{
				str+="[ireversible]";
			}
			str+="[affinityIncrease="+this.affinityIncrease+" (multiplicative)]";
		} else{
			str+="[type=direct]";
			str+="[dimerisationProbability="+this.dimerisationProbability+"]";
			if(this.isReversible){
				str+="[reversible]";
			} else{
				str+="[ireversible]";
			}
		
			str+="[affinityIncrease="+this.affinityIncrease;

			if(this.isAdditive){
				if(this.isFixed){
					str+="(additive-fixed)";
				} else{
					str+="(additive-non-fixed)";
				}
			} else{
				str+="(multiplicative)";
			}
			str+="]";
		}
		
		
		return str;
	}
}
