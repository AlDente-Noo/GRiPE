package objects;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import environment.Cell;

import utils.CellUtils;
import utils.Constants;
import utils.Utils;

/**
 * class that contains the description of TF species
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFspecies  implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -8664545127198799155L;
	public int id; // the id of the TF species
	public String name; // the id of the TF species
	
	public byte[] dbd; // the sequence regonise at DNA binding domain
	public PFM pfm;
	public String landscapeFile;
	public int landscapePosCol;
	public int landscapeAffinityColLR;
	public int landscapeAffinityColRL;
	public int landscapeEscapeLines;
	public String seqsFile;
	public double seqsDefaultValue;
	public int seqsEscapeLines;
	public String dbdFile;
	
	
	public int copyNumber; // number of molecules of this species
	public double es; // specific binding energy per nucleotide
	public int sizeLeft;
	public int sizeRight;
	public int sizeTotal;
	public double assocRate;
	public DNAregion initialDrop;
	public DNAregion relInitialDrop;
	public double maxAffinity;
	public double timeBoundAvg;
	public double timeBoundVar;
	public double[] timeReachedTargetSites; 
	public boolean isCognate;
	//public ArrayList<TargetSite> ts; // target sites on DNA
	
	public double unBindingProbability;
	public double slideLeftProbability;
	public double slideRightProbability;
	public double jumpingProbability;
	public double hopSTDdisplacement;
	public double specificWaitingTime;
	public int stepLeftSize;
	public int stepRightSize;
	public int uncorrelatedDisplacementSize;
	public boolean stallsIfBlocked;
	public double collisionUnbindingProbability;	
	public double affinityLandscapeRoughness;
	public double repressionRate;
	public double derepressionRate;
	public int repressionLeftSize;
	public int repressionRightSize;
	
	public double jumpNo;
	public double hopNo;
	public double slideLeftNo;
	public double slideRightNo;

	public double preboundProportion;
	public boolean preboundToHighestAffinity;

	public long countTFBindingEvents;
	public long countTFUnbindingEvents;
	public long countTFSlideLeftEvents;
	public long countTFSlideRightEvents;
	public long countTFHoppingEvents;
	public long countTFRepressionEvents;
	public long countTFDerepressionEvents;

	public long countTFforcedJumpsEvents;
	public long countTFHopsOutside;

	public double slidingEventsPerBinding;
	public double slidingLengthPerBinding;
	public double observedSlidingLengthPerBinding;
	public double residenceTimePerBinding;

	public boolean isCooperative;
	
	public boolean hasDNAbasedCooperativity;
	public boolean hasDirectCooperativity;

	public ArrayList<TFcooperativity> TFcoop;
	public int[][] isCooperativeSite;
	
	private HashMap<Integer,Integer> directCooperativeSpecies;
	
	public ArrayList<Integer> slidingLength;
	public ArrayList<Integer> observedSlidingLength;
	public ArrayList<Integer> slidingEvents;	
	
	
	public boolean isBiasedRandomWalk;
	public boolean isTwoStateRandomWalk;
	
	
	public boolean isImmobile;
	
	//public double moveRateThreshold;
	
	/**
	 * class constructor for the random tfs
	 * @param DNAstrand
	 * @param pos
	 * @param dbdLength
	 * @param copyNumber
	 * @param ens
	 * @param es
	 */
	public TFspecies(int id, byte[] DNAstrand, int pos, int dbdLength, int copyNumber, double es, DNAregion dnaRegion, int sizeLeft, int sizeRight, double assocRate, DNAregion initialDrop, boolean isCognate, double unBindingProbability, double slideLeftProbability, double slideRightProbability, double jumpingProbability, double hopSTDdisplacement, double specificWaitingTime, int stepLeftSize, int stepRightSize, int uncorrelatedDisplacementSize, boolean stallsIfBlocked, double collisionUnbindingProbability, double affinityLandscapeRoughness, double preboundProportion, boolean preboundToHighestAffinity, boolean isImmobile, int TFdirections, boolean isBiasedRandomWalk, boolean isTwoStateRandomWalk, double repressionRate, double derepressionRate, int repressionLenLeft, int repressionLenRight){
		this.id = id;
		name="TF"+id;
		this.es = es;
		this.copyNumber = copyNumber;
		pfm=null;
		this.landscapeFile ="";
		this.seqsFile = "";
		dbd = new byte[dbdLength];
		maxAffinity = 0;
		for(int i = pos; i < pos + dbdLength && i < DNAstrand.length;i++){
			dbd[i-pos] = DNAstrand[i];
			if(dbd[i-pos]!=CellUtils.bpANYID){
				maxAffinity+=es;
			}
		}
		this.sizeLeft = sizeLeft;
		this.sizeRight = sizeRight;
		this.sizeTotal = sizeLeft + sizeRight+dbdLength;
		this.assocRate = assocRate;
		this.timeBoundAvg = 0;
		this.timeBoundVar = 0;
		
		/*ts=new ArrayList<TargetSite>();
		ts.add(new TargetSite("",pos,pos+dbdLength, dnaRegion, id, this.sizeTotal, DNAstrand.length,TFdirections));*/

		this.initialDrop = initialDrop;
		this.isCognate = isCognate;
		
		this.unBindingProbability=unBindingProbability;
		this.slideLeftProbability=slideLeftProbability;
		this.slideRightProbability=slideRightProbability;
		this.jumpingProbability=jumpingProbability;
		this.repressionRate = repressionRate;
		this.derepressionRate = derepressionRate;
		this.repressionLeftSize=repressionLenLeft;
		this.repressionRightSize=repressionLenRight;
		this.hopSTDdisplacement=hopSTDdisplacement;
		this.specificWaitingTime=specificWaitingTime;
		this.stepLeftSize=stepLeftSize;
		this.stepRightSize=stepRightSize;
		this.uncorrelatedDisplacementSize=uncorrelatedDisplacementSize;
		this.stallsIfBlocked =stallsIfBlocked;
		this.collisionUnbindingProbability=collisionUnbindingProbability;	
		this.affinityLandscapeRoughness = affinityLandscapeRoughness;
		this.preboundProportion = preboundProportion;
		this.preboundToHighestAffinity = preboundToHighestAffinity;
		jumpNo = this.unBindingProbability*this.jumpingProbability;
		hopNo = this.unBindingProbability;
		slideLeftNo = hopNo+this.slideLeftProbability;
		slideRightNo = slideLeftNo + this.slideRightProbability;
		
	
		this.countTFBindingEvents =0;
		this.countTFUnbindingEvents =0;		
		this.countTFSlideLeftEvents =0;
		this.countTFSlideRightEvents =0;
		this.countTFHoppingEvents =0;
		this.countTFforcedJumpsEvents=0;
		this.countTFHopsOutside=0;
		this.countTFRepressionEvents=0;
		
		this.isCooperative = false;
		hasDNAbasedCooperativity = false;
		hasDirectCooperativity = false;
		TFcoop = new ArrayList<TFcooperativity>();
		directCooperativeSpecies = new HashMap<Integer, Integer>();
		slidingLength = new ArrayList<Integer>();
		slidingEvents = new ArrayList<Integer>();
		observedSlidingLength = new ArrayList<Integer>();
		//moveRateThreshold = Constants.NONE;
		
		this.isImmobile=isImmobile;
		
		this.isBiasedRandomWalk = isBiasedRandomWalk;
		this.isTwoStateRandomWalk = isTwoStateRandomWalk;
		this.dbdFile="";
		
	}
	
	
	//, ArrayList<TargetSite> ts
	/**
	 * class constructor
	 * @param id the id of the TF
	 * @param name the name
	 * @param dbd the dbd sequence
	 * @param copyNumber the copyNumber
	 * @param ens the non specific energy
	 * @param es the specific energy
	 * @param ts the list of target sites
	 */
	public TFspecies(DNAregion dnaRegion, int id, String name, String dbd, int copyNumber, double es, int sizeLeft, int sizeRight, double assocRate, DNAregion initialDrop, boolean isCognate, double unBindingProbability, double slideLeftProbability, double slideRightProbability, double jumpingProbability, double hopSTDdisplacement, double specificWaitingTime, int stepLeftSize, int stepRightSize, int uncorrelatedDisplacementSize, boolean stallsIfBlocked, double collisionUnbindingProbability, double affinityLandscapeRoughness, double preboundProportion, boolean preboundToHighestAffinity, boolean isImmobile, boolean isBiasedRandomWalk, boolean isTwoStateRandomWalk, double repressionRate, double derepressionRate, int repressionLenLeft, int repressionLenRight, Cell n){
		this.id = id; 
		this.name=name; 
		this.isCognate = false;

		this.dbd = new byte[0];
		int sizeInBP = 0;
		pfm=null;
		this.landscapeFile ="";
		this.seqsFile = "";
		this.dbdFile="";

		
		this.sizeLeft = sizeLeft;
		this.sizeRight = sizeRight;
		
		sizeInBP = parseDBD(n,dbd);
		
		this.copyNumber = copyNumber; 
		this.es=es; 
		

		this.sizeTotal = sizeLeft + sizeRight+ sizeInBP;
		this.assocRate = assocRate;
		this.timeBoundAvg = 0;
		this.timeBoundVar = 0;
			
		/*this.ts=new ArrayList<TargetSite>(); // target sites on DNA
		for(TargetSite t:ts){
			this.ts.add(t);
		}*/
		
		//initial drop
		relInitialDrop=new DNAregion("",initialDrop.start-dnaRegion.start, initialDrop.end-dnaRegion.start);
		if(relInitialDrop.start > 0 && relInitialDrop.end-1 < dnaRegion.size()){
			this.initialDrop=initialDrop;
		} else{
			this.initialDrop = new DNAregion("",0, dnaRegion.size());
			this.relInitialDrop = new DNAregion("",0, dnaRegion.size());
		}
	
		
		this.unBindingProbability=unBindingProbability;
		this.slideLeftProbability=slideLeftProbability;
		this.slideRightProbability=slideRightProbability;
		this.jumpingProbability=jumpingProbability;
		this.repressionRate = repressionRate;
		this.derepressionRate = derepressionRate;
		this.repressionLeftSize=repressionLenLeft;
		this.repressionRightSize=repressionLenRight;
		this.hopSTDdisplacement=hopSTDdisplacement;
		this.specificWaitingTime=specificWaitingTime;
		this.stepLeftSize=stepLeftSize;
		this.stepRightSize=stepRightSize;
		this.uncorrelatedDisplacementSize=uncorrelatedDisplacementSize;
		this.stallsIfBlocked =stallsIfBlocked;
		this.collisionUnbindingProbability=collisionUnbindingProbability;	
		this.affinityLandscapeRoughness = affinityLandscapeRoughness;
		this.preboundProportion = preboundProportion;
		this.preboundToHighestAffinity = preboundToHighestAffinity;
		jumpNo = this.unBindingProbability*this.jumpingProbability;
		hopNo = this.unBindingProbability;
		slideLeftNo = hopNo+this.slideLeftProbability;
		slideRightNo = slideLeftNo + this.slideRightProbability;
		
		//System.out.println("jumpNo="+jumpNo+"; hopNo="+hopNo+"; slideLeftNo="+slideLeftNo+"; slideRightNo="+slideRightNo);
		
		
		this.countTFBindingEvents =0;
		this.countTFUnbindingEvents =0;		
		this.countTFSlideLeftEvents =0;
		this.countTFSlideRightEvents =0;
		this.countTFHoppingEvents =0;
		this.countTFforcedJumpsEvents=0;
		this.countTFHopsOutside=0;
		this.countTFRepressionEvents=0;
		
		this.isCooperative = true;
		TFcoop = new ArrayList<TFcooperativity>();
		directCooperativeSpecies = new HashMap<Integer, Integer>();

		this.isImmobile = isImmobile;
		
		slidingLength = new ArrayList<Integer>();
		observedSlidingLength = new ArrayList<Integer>();
		slidingEvents = new ArrayList<Integer>();
		//moveRateThreshold = Constants.NONE;
		this.isBiasedRandomWalk = isBiasedRandomWalk;
		this.isTwoStateRandomWalk = isTwoStateRandomWalk;
		

	}

	public boolean isRepressor() {
		return this.repressionRate > 0.0;
	}

	/**
	 * parse a String into DBD
	 * @param n the cell
	 * @param dbd the DNA binding domain text 
	 * @return
	 */
	public int parseDBD(Cell n, String dbd){
		int sizeInBP=0;
		if(dbd.startsWith(Constants.DBD_TYPE_SEQ)){
			ArrayList<Byte> bufferDBD;
			dbd = dbd.replaceAll(Constants.DBD_TYPE_SEQ, "");
			bufferDBD = CellUtils.getSequenceIDs(dbd);
			
			this.dbd =  new byte[bufferDBD.size()];
			this.maxAffinity  = 0;
			for(int i=0;i<bufferDBD.size();i++){
				this.dbd[i]=bufferDBD.get(i);
				if(this.dbd[i]!=CellUtils.bpANYID){
					maxAffinity+=es;
				}
			}
			sizeInBP = this.dbd.length;
			if(this.dbd.length>0){
				this.isCognate = true;
			}
		}  else if(dbd.startsWith(Constants.DBD_TYPE_PFM) || dbd.startsWith(Constants.DBD_TYPE_PWM)){
			pfm = new PFM(dbd,n);
			if(this.pfm.isCorrect){
				this.isCognate = true;
			}
			sizeInBP = this.pfm.motifSize;
		}  else if(dbd.startsWith(Constants.DBD_TYPE_LANDSCAPE)){
			String buffer = dbd.replaceAll(Constants.DBD_TYPE_LANDSCAPE, "");
			this.dbdFile = dbd;
			if(buffer!=null && !buffer.isEmpty() && buffer.contains(Constants.DBD_TYPE_SEPARATOR)){
				
				String[] bufferStr = buffer.split(Constants.DBD_TYPE_SEPARATOR);
				
				if(bufferStr.length>= 6){
					

					//is cognate
					this.isCognate = Utils.parseBoolean(bufferStr[1], true);

					//DBD size in bp
					sizeInBP = Utils.parseInteger(bufferStr[0], Constants.NONE);
					if(sizeInBP<0 || (sizeInBP==0 && sizeLeft==0 && sizeRight==0)){
						n.stopSimulation("Could not load the affinity landscape "+dbd+": no motif size");
					}
					
					
					if(bufferStr.length== 7){
						this.landscapePosCol = Utils.parseInteger(bufferStr[2],Constants.NONE);
						this.landscapeAffinityColLR = Utils.parseInteger(bufferStr[3],Constants.NONE);
						this.landscapeAffinityColRL = Utils.parseInteger(bufferStr[4],Constants.NONE);
						this.landscapeFile = bufferStr[5];
						this.landscapeEscapeLines = Utils.parseInteger(bufferStr[6],0);
						this.landscapePosCol--;
						if(this.landscapePosCol <0){
							n.stopSimulation("Could not load the affinity landscape "+dbd+": no column for DNA position");

						}
					} else{
						this.landscapePosCol = Constants.NONE;
						this.landscapeAffinityColLR = Utils.parseInteger(bufferStr[2],Constants.NONE);
						this.landscapeAffinityColRL = Utils.parseInteger(bufferStr[3],Constants.NONE);
						this.landscapeFile = bufferStr[4];
						this.landscapeEscapeLines = Utils.parseInteger(bufferStr[5],0);
					}
					this.landscapeAffinityColLR --;
					this.landscapeAffinityColRL--;

					if(this.landscapeAffinityColLR <0){
						n.stopSimulation("Could not load the affinity landscape "+dbd+": no column for affinity in 3'-5' direction");

					}

					if(this.landscapeAffinityColRL <0){
						n.stopSimulation("Could not load the affinity landscape "+dbd+": no column for affinity in 5'-3' direction");
					}
					
					if(this.landscapeFile==null || this.landscapeFile.isEmpty()){
						n.stopSimulation("Could not load the affinity landscape "+dbd+": no file containing the landscape");

					} else{
						File f=new File(this.landscapeFile );
						if(!f.exists() || f.isDirectory() || !f.canRead() || f.length()==0){
							n.stopSimulation("Could not load the affinity landscape "+dbd+": cannot find the landscape file: "+landscapeFile);

						} 
					}
					
					
				} 
			
				
			} else{
				n.stopSimulation("Could not load the affinity landscape "+dbd+": could not split the text by "+Constants.DBD_TYPE_SEPARATOR);
			} 
			
			//if(!isLandscapeSpecificationCorrect){
			//	n.stopSimulation("Could not load the affinity landscape: "+dbd);
			//}
			
		}  else if(dbd.startsWith(Constants.DBD_TYPE_SEQS)){
			//System.out.println("load sequences files");
			String buffer = dbd.replaceAll(Constants.DBD_TYPE_SEQS, "");
			this.dbdFile = dbd;

			if(buffer!=null && !buffer.isEmpty() && buffer.contains(Constants.DBD_TYPE_SEPARATOR)){
				
				String[] bufferStr = buffer.split(Constants.DBD_TYPE_SEPARATOR);
				
				if(bufferStr.length== 5){
					

					//is cognate
					this.isCognate = Utils.parseBoolean(bufferStr[1], true);

					//DBD size in bp
					sizeInBP = Utils.parseInteger(bufferStr[0], Constants.NONE);
					if(sizeInBP<0 || (sizeInBP==0 && sizeLeft==0 && sizeRight==0)){
						n.stopSimulation("Could not load the affinity landscape "+dbd+": no motif size");
					}
					
					
	
					this.seqsDefaultValue = Utils.parseDouble(bufferStr[2],Constants.NONE);
					this.seqsFile = bufferStr[3];
					this.seqsEscapeLines = Utils.parseInteger(bufferStr[4],0);
					
					if(this.seqsDefaultValue <0){
						n.stopSimulation("Could not load the sequence file "+dbd+": no default value");

					}

					
					if(this.seqsFile==null || this.seqsFile.isEmpty()){
						n.stopSimulation("Could not load the sequence file "+dbd+": no file containing the sequences");

					} else{
						File f=new File(this.seqsFile );
						if(!f.exists() || f.isDirectory() || !f.canRead() || f.length()==0){
							n.stopSimulation("Could not load the sequence file "+dbd+": cannot find the sequence file: "+seqsFile);

						} 
					}
					
					
				} 
			
				
			} else{
				n.stopSimulation("Could not load the sequence file"+dbd+": could not split the text by "+Constants.DBD_TYPE_SEPARATOR);
			} 
		}
		return sizeInBP;
		
	}
	
	
	/**
	 * @return returns the string of the current TF species
	 */
	public String toString(Cell n, boolean reduced){
		
		//"\"name\", \"DBD\",
		String str="\""+ name+"\", \"";
		if(pfm!=null && pfm.isCorrect && pfm.motifSize > 0){
			str+=pfm.toString(n.dna.bpFreq);
		} else if(dbd.length > 0){
			str+="SEQ:";
			for(int i=0;i<dbd.length;i++){
				str+=CellUtils.bps.bps[dbd[i]];
			}
		} else{
			str+=this.dbdFile;
		}
		
		//\"ES\", \"COPYNUMBER\", \"SIZELEFT\", \"SIZERIGHT\", \"ASSOCRATE\", \"INITIALDROP\", \"UNBINDINGPROBABILITY\", \"SLIDELEFTPROBABILITY\", \"SLIDERIGHTPROBABILITY\", \"JUMPINGPROBABILITY\", \"HOPSTDDISPLACEMENT\", \"SPECIFICWAITINGTIME\", \"STEPLEFTSIZE\", \"STEPRIGHTSIZE\", \"UNCORRELATEDDISPLACEMENTSIZE\", \"STALLSIFBLOCKED\", \"COLLISIONUNBINDPROBABILITY\", \"AFFINITYLANDSCAPEROUGHNESS\", \"PREBOUNDPROPORTION\", \"PREBOUNDTOHIGHESTAFFINITY\"
		str+="\","+es+", "+copyNumber+", "+this.sizeLeft+", "+this.sizeRight+", "+this.assocRate+", \""+this.initialDrop+"\"" + ", " +    unBindingProbability + ", " +    slideLeftProbability + ", " +    slideRightProbability + ", " +    jumpingProbability + ", " +    hopSTDdisplacement + ", " +    specificWaitingTime + ", " +    stepLeftSize + ", " +    stepRightSize + ", " +    uncorrelatedDisplacementSize + ", " + stallsIfBlocked + ", " +    collisionUnbindingProbability+ ", " +     affinityLandscapeRoughness+", "+preboundProportion+ ", "+preboundToHighestAffinity+ ", "+isImmobile;
		
		//"ISBIASEDRANDOMWALK\", \"ISTWOSTATERANDOMWALK\",
		str+= ", "+this.isBiasedRandomWalk +", "+ this.isTwoStateRandomWalk;

		//"REPRESSIONRATE", "DEREPRESSIONRATE", "REPRLENLEFT", "REPRLENRIGHT"
		str += ", " + this.repressionRate + ", " + this.derepressionRate + ", " + this.repressionLeftSize + ", " + this.repressionRightSize;

		if(!reduced){
			//"eventsBindingTotal\", \"eventsUnbindingTotal\", \"eventsSlideLeftTotal\", \"eventsSlideRightTotal\", \"eventsSlideTotal\", \"eventsHoppingTotal\", \"eventsForcedJumps\", \"eventsHopOutsideDNA\", \"collisionsCount\",		
			str+=  ", "+ (this.countTFBindingEvents) + ", " + (this.countTFUnbindingEvents)+ ", " + (this.countTFSlideLeftEvents) + ", " + (this.countTFSlideRightEvents) + ", " + (this.countTFSlideLeftEvents+this.countTFSlideRightEvents) + ", " + (this.countTFHoppingEvents)+ ", " + (this.countTFforcedJumpsEvents)+ ", " + (this.countTFHopsOutside)+ ", " + (n.dna.collisionsCountTotal);
	
			//     \"sizeTotal\"" , \"isCognate\"
			str+=  ", "+this.sizeTotal+", " + this.isCognate;
			
			
			//, timeBoundAvg \"residenceTimePerBinding\", \"slidingEventsPerBinding\", \"slidingLengthPerBinding\", \"observedSlidingLengthPerBinding\"";
			if(this.slidingLength!=null && !slidingLength.isEmpty()){
				//System.out.println(slidingLength);
				str+=", "+this.timeBoundAvg+", "+residenceTimePerBinding+", "+Utils.computeMean(slidingEvents)+", "+Utils.computeMean(this.slidingLength)+", "+Utils.computeMean(this.observedSlidingLength);
			} else {
				str+=", "+this.timeBoundAvg+", "+residenceTimePerBinding+", "+slidingEventsPerBinding+", "+this.slidingLengthPerBinding+", "+this.observedSlidingLengthPerBinding;
			}	
			str+=", \""+this.getCoopString()+"\"";

			str += ", " + this.countTFRepressionEvents + ", " + this.countTFDerepressionEvents;
		}
		return str;
	}

	/**
	 * @return returns the string header of the info of the current TF species
	 */
	//public String headerToString(){
	//	return "#id, name, copyNumber, es, dbd, sizeLeft, sizeRight, sizeTotal, assocRate, timeBoundAvg, residenceTimePerBinding, slidingEventsPerBinding, slidingLengthPerBinding, observedSlidingLength, initialDrop, isCognate,  unBindingProbability,  slideLeftProbability,  slideRightProbability,  jumpingProbability,  hopSTDdisplacement,  specificWaitingTime,  stepLeftSize,  stepRightSize,  uncorrelatedDisplacementSize,  stallsIfBlocked,  collisionUnbindingProbability, affinityLandscapeRoughness, preboundProportion, eventsBindingTotal, eventsUnbindingTotal, eventsSlideLeftTotal, eventsSlideRightTotal, eventsSlideTotal, eventsHoppingTotal, eventsForcedJumps, eventsHopOutsideDNA, collisionsCount, cooperativity, isBiasedRandomWalk, isTwoStateRandomWalk";
	//}
	
	public String headerToString(boolean reduced){
		String str = "\"name\", \"DBD\", \"ES\", \"COPYNUMBER\", \"SIZELEFT\", \"SIZERIGHT\", \"ASSOCRATE\", \"INITIALDROP\", \"UNBINDINGPROBABILITY\", \"SLIDELEFTPROBABILITY\", \"SLIDERIGHTPROBABILITY\", \"JUMPINGPROBABILITY\", \"HOPSTDDISPLACEMENT\", \"SPECIFICWAITINGTIME\", \"STEPLEFTSIZE\", \"STEPRIGHTSIZE\", \"UNCORRELATEDDISPLACEMENTSIZE\", \"STALLSIFBLOCKED\", \"COLLISIONUNBINDPROBABILITY\", \"AFFINITYLANDSCAPEROUGHNESS\", \"PREBOUNDPROPORTION\", \"PREBOUNDTOHIGHESTAFFINITY\", \"TFISIMMOBILE\", \"ISBIASEDRANDOMWALK\", \"ISTWOSTATERANDOMWALK\", \"REPRESSIONRATE\", \"DEREPRESSIONRATE\", \"REPRLENLEFT\", \"REPRLENRIGHT\"";
		if(!reduced){
			str+= ", \"eventsBindingTotal\", \"eventsUnbindingTotal\", \"eventsSlideLeftTotal\", \"eventsSlideRightTotal\", \"eventsSlideTotal\", \"eventsHoppingTotal\", \"eventsForcedJumps\", \"eventsHopOutsideDNA\", \"collisionsCount\", \"sizeTotal\", \"isCognate\", \"timeBoundAvg\", \"residenceTimePerBinding\", \"slidingEventsPerBinding\", \"slidingLengthPerBinding\", \"observedSlidingLengthPerBinding\", \"cooperativity\"";
			str+= ", \"eventsRepression\", \"eventsDerepression\"";
		}
		return str;
	}

	


	
	/**
	 * adds TF coopearativity
	 * @param TFcoop
	 */
	public void addCooperativity(TFcooperativity coop, int DNAlength, int TFdirections){
		this.TFcoop.add(coop);
		
		this.isCooperative = true;
		if(!this.hasDNAbasedCooperativity && coop.type == Constants.TF_COOPERATIVITY_TYPE_DNA){
			this.hasDNAbasedCooperativity = true;
		}
		if(!this.hasDirectCooperativity && coop.type == Constants.TF_COOPERATIVITY_TYPE_DIRECT){
			this.hasDirectCooperativity = true;
		}		
		

		
		
		//if there is a cooperative site then mark this on the array
		if(coop.type ==  Constants.TF_COOPERATIVITY_TYPE_DNA){
			
			if(this.isCooperativeSite == null){
				intiateCooperativeSite(DNAlength, TFdirections);
			}
			setSiteAsCooperative(coop.region0, coop.direction0, this.TFcoop.size()-1);
		}
		
		//if there is a cooperative site then mark this on the array
		if(coop.type ==  Constants.TF_COOPERATIVITY_TYPE_DIRECT){
			directCooperativeSpecies.put(coop.species1ID, this.TFcoop.size()-1);
		}
		
	}
	
	/**
	 * 
	 * @param position
	 * @return
	 */
	//public boolean isAtCooperativeSite(int position){
	//	for 
	//}
	
	/**
	 * initialise the cooperative site
	 */
	private void intiateCooperativeSite(int DNAlength, int TFdirections){
		isCooperativeSite = new int[DNAlength][TFdirections];
		for(int i=0;i<DNAlength;i++){
			for(int j=0;j<TFdirections;j++){
				this.isCooperativeSite[i][j] = Constants.NONE;
			}
		}
	}
	
	
	/**
	 * sets a site as cooperative
	 * @param region
	 */
	private void setSiteAsCooperative(DNAregion region, int direction, int coopID){
		int start=(int) Math.max(0,region.start);
		int end=(int) Math.min(1,region.end);
		end = Math.max(start+1, end);
		
		int startDir = 0, endDir=this.isCooperativeSite[0].length;
		if(this.TFcoop.get(coopID).direction0!=Constants.NONE){
			startDir = this.TFcoop.get(coopID).direction0;
			endDir = this.TFcoop.get(coopID).direction0+1;
		}
			
		for(int i=start;i<end;i++){
			for(int j=startDir;j<endDir;j++){
				this.isCooperativeSite[i][j]=coopID;
				
				//System.out.println("for TF species "+this.id+" position "+i+" in direction "+j+" was marked as generating coop");
				
			}
		}
	}
	
	
	/**
	 * returns true if this site is marked as being cooperative
	 * @param position
	 * @return
	 */
	public boolean isCooperativeSite(int position, int dir){	
		return this.isCooperativeSite[position][dir]!=Constants.NONE;
	}
	
	/**
	 * returns the species that is affected by binding this molecule to region0
	 * @param position
	 * @return
	 */
	public TFcooperativity getCooperativity(int position, int direction){
		return this.TFcoop.get(this.isCooperativeSite[position][direction]);
	}

	/**
	 * prints the sliding lengths to a file
	 * @param filename
	 */
	public void printSlidingLengths(String path, String filename){
		try {
		   BufferedWriter out = new BufferedWriter(new FileWriter(path+filename));
		    out.write("\"no\", \"slidingLength\", \"slidingEvents\" \n");
			StringBuffer strBuf =  new StringBuffer();
			
			for(int i=0;i< this.slidingLength.size();i++){
				strBuf.delete(0,strBuf.length());				
					strBuf.append(i);
					strBuf.append(", ");
					strBuf.append(this.slidingLength.get(i));
					strBuf.append(", ");
					strBuf.append(this.slidingEvents.get(i));
				out.write(strBuf.toString());
				out.newLine();
			}
		    
		    out.close();
		} catch (IOException e) {
		}
		
		
	}
	
	/**
	 * prints the sliding lengths to a file
	 * @param filename
	 */
	public void printObservedSlidingLengths(String path, String filename){
		try {
		   BufferedWriter out = new BufferedWriter(new FileWriter(path+filename));
			StringBuffer strBuf =  new StringBuffer();
			
		    out.write("\"no\", \"observedSlidingLength\" \n");

			
			for(int i=0;i<this.observedSlidingLength.size();i++){
				strBuf.delete(0,strBuf.length());			
				strBuf.append(i);
				strBuf.append(", ");
				strBuf.append(this.observedSlidingLength.get(i));
				out.write(strBuf.toString());
				out.newLine();
			}
		    
		    out.close();
		} catch (IOException e) {
		}
		
		
	}
	
	/**
	 * returns the id of the cooperativity for a right neighbour
	 * @param speciesID1
	 * @param direction0
	 * @param direction1
	 * @return
	 */
	public TFcooperativity getDirectCooperativityRight(int speciesID1, int direction0, int direction1){
		TFcooperativity result = null;
		
		if(directCooperativeSpecies.containsKey(speciesID1)){
			int coopID=directCooperativeSpecies.get(speciesID1);
			if((TFcoop.get(coopID).direction0==Constants.NONE || TFcoop.get(coopID).direction0==direction0) && (TFcoop.get(coopID).direction1==Constants.NONE || TFcoop.get(coopID).direction1==direction1)){
				result=TFcoop.get(coopID);
			}
		}
		
		return result;
	}
	
	
	/**
	 * returns the id of the cooperativity for a left neighbour
	 * @param speciesID1
	 * @param direction0
	 * @param direction1
	 * @return
	 */
	public TFcooperativity getDirectCooperativityLeft(int speciesID1, int direction0, int direction1){
		TFcooperativity result = null;
		
		if(directCooperativeSpecies.containsKey(speciesID1)){
			int coopID=directCooperativeSpecies.get(speciesID1);
			if((TFcoop.get(coopID).direction0==Constants.NONE || TFcoop.get(coopID).direction0==direction1) && (TFcoop.get(coopID).direction1==Constants.NONE || TFcoop.get(coopID).direction1==direction0)){
				result=TFcoop.get(coopID);
			}
		}
		
		return result;
	}
	
	/**generates the string description of the cooperativity 
	 * 
	 * @return
	 */
	private String getCoopString(){
		String str="";
		if(this.TFcoop==null || this.TFcoop.isEmpty()){
			str="none";
		} else{
			for(int i=0;i<this.TFcoop.size();i++){
				str+="; "+this.TFcoop.get(i).toString();
			}
			str=str.substring(2);
		}
		return str;
	}
	
	
	/*public boolean updateTStimeReached(int position, int direction, double time){
		boolean found=false;
		for(int i=0;i<ts.size();i++){
			for(int j=ts.get(i).relStart;j<ts.get(i).relEnd;j++){	
				if(j==position){
					ts.get(i).timeReached[j-ts.get(i).relStart][direction] += time;
					ts.get(i).countReached[j-ts.get(i).relStart][direction]++;
					found = true;				
				}
			}
		}
		
		return found;
	}*/
	
}

