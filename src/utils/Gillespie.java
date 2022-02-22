package utils;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;


/**
 * this class implements the methods of the Gillespie Direct Method SSA
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class Gillespie {
	/**
	 * computes the next reaction time according to Gillespie algorithm
	 * @param propensitySum the sum of all propensities
	 * @param generator  a random number generator
	 * @return the time of the next reaction
	 */
	public static double computeNextReactionTime(double propensitySum, Random generator){
		return (1/propensitySum)*Math.log(1/generator.nextDouble());
	}

	/**
	 * computes the next reaction time according to Gillespie algorithm
	 * @param waitingTimeSum the sum of all waiting times
	 * @param generator  a random number generator
	 * @return the time of the next reaction
	 */
	public static double computeNextReactionTimeFromWaitingTime(double waitingTimeSum, Random generator){
		return (waitingTimeSum)*Math.log(1/generator.nextDouble());
	}
	
	/**
	 * generate the next reaction to take place
	 * @param value the propensity sum multiplied by a random number
	 * @return the index of the next mRNA species likely to attract a ribosome
	 */
	public static  int getNextReaction(double value,double[] ribosomeBindingPropensity){
		double sum=0;

		for(int i=0;i<ribosomeBindingPropensity.length;i++){
			if(ribosomeBindingPropensity[i]>0){
				sum+=ribosomeBindingPropensity[i];
				if(sum>=value){
					return i;
				}
			}
		}
		
		//double stored error return the last value
		if(ribosomeBindingPropensity.length>0){
			return ribosomeBindingPropensity.length-1;
		}

		throw new RuntimeException("Error with computing next reaction sum:"+value+"\nPropensities: ");
	}
	
	
	/**
	 * generate the next reaction to take place
	 * @param value the propensity sum multiplied by a random number
	 * @return the index of the next mRNA species likely to attract a ribosome
	 */
	public static  int getNextReaction(double value, int[] reactionPropensity){
		long sum=0;
		for(int i=0;i<reactionPropensity.length;i++){
			if(reactionPropensity[i]>0){
				sum+=reactionPropensity[i];
				if(sum>=value){
					return i;
				}
			}
		}
		

		throw new RuntimeException("Error with computing next reaction sum:"+value+"\nPropensities: "+sum);
	}
	

	
	
	/**
	 * generate the next reaction to take place
	 * @param value the propensity sum multiplied by a random number
	 * @return the index of the next mRNA species likely to attract a ribosome
	 */
	public static  int getNextReaction(double value, boolean[] reactionPropensity){
		long sum=0;
		
		for(int i=0;i<reactionPropensity.length;i++){
			if(reactionPropensity[i]){
				sum++;
				if(sum>=value){
					return i;
				}
			}
		}
	

		throw new RuntimeException("Error with computing next reaction sum:"+value+"\nPropensities: "+sum);
	}
	
	
	
	 /** generate the next position to bind on DNA based on Gillespie algorithm
	 * @param value the propensity sum multiplied by a random number
	 * @return the index of the next mRNA species likely to attract a ribosome
	 */
	public static  int getNextReaction(double value, double[] propensity, double[] sectorSum, int sectorSize){
		double sum=0;
		
		for(int i=0;i<sectorSum.length;i++){
			if(sectorSum[i]+sum >= value){	
				for(int j=i*sectorSize;j< (int) Math.min(propensity.length, (i+1)*sectorSize); j++){
					sum+=propensity[j];
					if(sum>=value){
							return j;
						}
					
				}
			} 
			sum+=sectorSum[i];	
		}
		//double stored error return the last value
		if(propensity.length>0){
			return propensity.length-1;
		}

		throw new RuntimeException("Error with computing next reaction sum:"+value+"\nPropensities: "+sum);
	}
	
	
	
	/**
	 * generate the next position to bind on DNA based on Gillespie algorithm
	 * @param value the propensity sum multiplied by a random number
	 * @return the index of the next mRNA species likely to attract a ribosome
	 */
	public static  int getNextPositionToBind(double value, boolean[] availability, int[] sectorSum, int sectorSize){
		long sum=0;
		
		for(int i=0;i<sectorSum.length;i++){
			if(sectorSum[i]+sum >= value){	
				for(int j=i*sectorSize;j< (int) Math.min(availability.length, (i+1)*sectorSize); j++){
					if(availability[j]){
						sum++;
						if(sum>=value){
							return j;
						}
					}
				}
			} 
			sum+=sectorSum[i];	
		}
		
		throw new RuntimeException("Error with computing next reaction sum:"+value+"\nPropensities: "+sum);
	}
	
	
	
	
	/**
	 * generate the next reaction to take place
	 * @param value the propensity sum multiplied by a random number
	 * @return the index of the next mRNA species likely to attract a ribosome
	 */
	public static int getNextReaction(double propensitySum,HashMap<Integer,Double> propensity, double reactionRate){
		double sum=0;
		Integer key;
		Iterator<Integer> it = propensity.keySet().iterator();
	    while (it.hasNext()) {
	    	key =  it.next();
	    	sum+=propensity.get(key);
			if(sum>=propensitySum){
				return key;
			}
	    }

		throw new RuntimeException("Error with computing next reaction:");
	}
	

}
