package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import objects.TargetSite;
import objects.TargetSitesAndGroups;
import objects.TargetSitesGroup;
import environment.Cell;

public class TSfileParser {

	public TargetSitesAndGroups tsg;
	private int TSid;
	private int TSGid;
	
	public TSfileParser(String filename, Cell n){
	
	
		//TARGETSITES
		tsg=new TargetSitesAndGroups();			
		BufferedReader reader = null;
		TSid = 0;
		TSGid = 0;
        try{
            reader = new BufferedReader(new FileReader(filename));
            String text = null;            
            while ((text = reader.readLine()) != null){
	            	text=text.trim();
	            	if(!text.isEmpty()){
	            		extractTargetSitesAndGroups(n,text);		
	            	}
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            
        		try{
                if (reader != null){
                    reader.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
	
			
	}
	
	
	/**
	 * extracts the target site
	 * @param n the cell object
	 * @param str the text to be parsed
	 * @param dnaRegion the DNA region of the subsequence analysed
	 * @param TFid the ID of the TF species
	 * @param TFsize the size of the TF 
	 * @param DNAsize the size of the DNA
	 * @param TFdirections the direction of the TS
	 * @return
	 */
	public TargetSite extractTargetSite(Cell n, int targetSiteID, String str){
		TargetSite bufferTS = new TargetSite(n,targetSiteID, str, "", Constants.NONE, Constants.NONE);
		if(bufferTS.region.start!=Constants.NONE && bufferTS.region.end!=Constants.NONE ){
			n.printDebugInfo("target site " + bufferTS+" loaded");
		} else{
			n.printDebugInfo("error on parsing target site " +bufferTS);
			bufferTS = null;
		}
		
		return bufferTS;
	}


	
	/**
	 * parse one line from the TS file 
	 * @param n
	 * @param text
	 * @return
	 */
	public void extractTargetSitesAndGroups(Cell n, String text){
		TargetSitesGroup bufferTSG;
		bufferTSG = new TargetSitesGroup(TSGid, text);
		TargetSite bufferTS;
		//sitesStr = new ArrayList<String>();
		String[] sitesStr;
		String bufferLogicExpression;
		int newID;

		sitesStr = text.trim().split("[ \t()]+");
		bufferLogicExpression=text;
		for(String str:sitesStr){
			str = str.trim();
			if(!str.isEmpty() &&  !CellUtils.tgsl.isOperator(str)){
				bufferTS = extractTargetSite(n,TSid, str);
				if(bufferTS!=null){
					newID = tsg.addTargetSite(bufferTS);
					bufferTSG.addTargetSite(newID);
					tsg.ts.get(newID).group.add(TSGid);
					if(newID==TSid){
						TSid++;
					}
					bufferLogicExpression = bufferLogicExpression.replaceAll(str, newID+"");
				}
			}
		}
		bufferLogicExpression = CellUtils.tgsl.replaceOperators(bufferLogicExpression);
		
		bufferTSG.generateRPN(bufferLogicExpression);
		/*boolean[] occupancy = new boolean[2];
		occupancy[0]=false; occupancy[1]=false;System.out.println("(false AND false) OR (NOT false) = "+bufferTSG.evaluateRPNTree(occupancy) );
		occupancy[0]=false; occupancy[1]=true;System.out.println("(false AND true) OR (NOT false) = "+bufferTSG.evaluateRPNTree(occupancy) );
		occupancy[0]=true; occupancy[1]=false;System.out.println("(true AND false) OR (NOT true) = "+bufferTSG.evaluateRPNTree(occupancy) );
		occupancy[0]=true; occupancy[1]=true;System.out.println("(true AND true) OR (NOT true) = "+bufferTSG.evaluateRPNTree(occupancy) );*/
		//System.out.println("logic:"+bufferLogicExpression);
		tsg.tsg.add(bufferTSG);
		TSGid++;
	}
	
}
