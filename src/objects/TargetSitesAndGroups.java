package objects;

import java.io.Serializable;
import java.util.ArrayList;

import environment.Cell;

import utils.Constants;

/**
 * contains a list of target sites and a list of the target sites groups
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TargetSitesAndGroups   implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -8051357555423486023L;
	public ArrayList<TargetSite> ts;
	public ArrayList<TargetSitesGroup> tsg;
	public boolean[] occupancy;
	
	/**
	 * class constructor
	 */
	public TargetSitesAndGroups(){
		ts = new ArrayList<TargetSite>();
		tsg = new  ArrayList<TargetSitesGroup>();
	}
	
	/**
	 * class constructor
	 */
	public TargetSitesAndGroups(String text){
		ts = new ArrayList<TargetSite>();
		tsg = new  ArrayList<TargetSitesGroup>();
	}
	
	/**
	 * this is true if no TS is created
	 * @return
	 */
	public boolean isEmpty(){
		return ts.isEmpty();
	}
	
	
	/**
	 * adds a target site i
	 * if the target site does not exists is added to the list and its ID is returned otherwise the ID of the existing target site is returned;
	 * @param ts the target site
	 */
	public int addTargetSite(TargetSite ts){
		int result=Constants.NONE;
		
		for(int i=0;i<this.ts.size() && result==Constants.NONE; i++){
			if(this.ts.get(i).eqauls(ts)){
				result=i;
			}
		}
		
		if(result==Constants.NONE){
			ts.targetSiteID=this.ts.size();
			this.ts.add(ts);
			occupancy = new boolean[this.ts.size()];
			for(int i=0;i<occupancy.length;i++){
				occupancy[i] = false;
			}
			result=this.ts.size()-1;
		}
		
		return result;
		
	}
	
	
	/**
	 * generate the text string
	 */
	
	public String toString(){
		String text="";
		
		/*for(int i=0; i< ts.size();i++){
			text+="\""+ts.get(i).toString()+"\", "+ts.get(i).firstTimeReached+", "+ts.get(i).timesReached+", "+ts.get(i).timeOccupied+"\n";
		}*/
		
		for(int i=0; i< tsg.size();i++){
			text+="\""+tsg.get(i).toString()+"\", "+tsg.get(i).firstTimeReached+", "+tsg.get(i).timesReached+", "+tsg.get(i).timeOccupied+"\n";
		}
		
		return text;
	}
	
	
	/** 
	 * adds a new group
	 * @param n the cell
	 * @param pos the position on the DNA of the target site
	 * @param TFid the ID of the TF
	 */
	public void addGroup(Cell n, int pos, int TFid){
		int newID = this.addTargetSite(new TargetSite(n, this.ts.size(), "", pos, pos+1, n.dna.region, TFid, n.TFspecies[TFid].sizeTotal, n.dna.strand.length, n.TFreadingDirection));
		TargetSitesGroup bufferTSG = new TargetSitesGroup(tsg.size(),ts.get(newID).toString());
		bufferTSG.addTargetSite(newID);
		bufferTSG.generateRPN(newID+"");
		bufferTSG.text = this.ts.get(newID).toString();
		tsg.add(bufferTSG);
	}
	
	/**
	 * updates statistics about a target site
	 * @param tsID the ID of the target site
	 * @param time the current simulation time
	 * @param lastTime the last time the molecule made an action
	 * @param position the lcurrent position of the TF

	 */
	public void updateTargetSiteStatistics(int tsID, double time, boolean bound, double timeBound){

		this.occupancy[tsID] = bound;
		
		//check for each group the current TF belongs to the state
		boolean evaluateGroup;
		for(int i:ts.get(tsID).group){
			evaluateGroup = tsg.get(i).evaluateRPNTree(occupancy);
			if(evaluateGroup && !tsg.get(i).isOccupied){
				tsg.get(i).updateTimesReachedStatistics(time);
			} 
			
			if(tsg.get(i).isOccupied){
				tsg.get(i).updateOccupancyStatistics(time);
			}
			
			tsg.get(i).isOccupied = evaluateGroup;
			tsg.get(i).lastTimeUpdate = time;
		}
		
	}
	
	
	/**
	 * verifies if there are target sites to be reached
	 * @return true if there are unreached site and false otherwise
	 */
	public boolean areTargetSitesToBeReached(){
		boolean result=false;
		
		for(int tsID=0;tsID<ts.size() && !result;tsID++){
			for(int i=0;i<ts.get(tsID).group.size() && !result; i++){
				if(tsg.get(i).timesReached == 0){
					result=true;
				}
			}
		}
		
		return result;
		
	}
	
	/**
	 * returns an String with the names of all target site groups 
	 * @return
	 */
	public String getTargetSiteGroupsString(){
		
		String text="";
		for(int i=0; i< tsg.size();i++){
			if(i>0){
				text+=", ";
			}
			text+="\""+tsg.get(i).text+"\"";
		}
		
		return text;
	}
	
	/**
	 * returns an String with the occupancy of all target site groups 
	 * @return
	 */
	public String getTargetSiteGroupsOccupancyString(){
		
		String text="";
		for(int i=0; i< this.tsg.size();i++){
			if(i>0){
				text+=", ";
			}
			if(this.tsg.get(i).isOccupied){
				text+=" 1";
			} else{
				text+=" 0";
			}
		}
		
		return text;
	}
}
