package utils;

import java.util.ArrayList;

import environment.Cell;

import objects.DNAregion;
import objects.TFcooperativity;
import objects.TFcooperativityDirectAdditive;
import objects.TFcooperativityDirectAdditiveFixed;
import objects.TFcooperativityDirectMultiplicative;

/**
 * class that parses a TF cooperativity file
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFcooperativityFileParser {
	public ArrayList<TFcooperativity> data;
	public boolean loaded;
	public boolean hasVariousRowSizes;
	public boolean hasAllHeaders;
	public boolean parsed;
	
	
	
	public TFcooperativityFileParser(Cell n, String filename, int defaultDirectionOf0, int defaultDirectionOf1, DNAregion defaultDNAregion, double defaultDimerisationProbability, DNAregion region){
		parsed = false;
		data = new  ArrayList<TFcooperativity>();
		
		
		CSVparser csv = new CSVparser(filename, true, true);
		

		hasAllHeaders = false;
		// manage to parse csv
		this.loaded = csv.loaded;
		this.hasVariousRowSizes = csv.hasVariousRowSizes;

		if(csv.loaded && !csv.hasVariousRowSizes){
			hasAllHeaders = true;
			// check header
			for(String key: Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER){
				if(!csv.header.containsKey(key)){
					hasAllHeaders = false;
					n.printDebugInfo("Error when parsing the TF cooperativity csv file "+filename+": missing column "+key);
				}
			}
			

			
			if(hasAllHeaders){
				if(!csv.data.isEmpty()){
					parsed=true;
				}
				
				String cellContent;
				
				int species0ID;
				int species1ID;
				int type;
				int directionOf0;
				int directionOf1;
				DNAregion DNAregion0;
				DNAregion DNAregion1;
				double dimerisationProbability;
				double affinityIncrease;
				boolean isReversible = true;
				int rowID=0;
				
				boolean isAdditive;
				boolean isFixed;
				
				for(ArrayList<String> buffer: csv.data){
					//init params
					species0ID=Constants.NONE;
					species1ID=Constants.NONE;
					type=Constants.NONE;
					directionOf0=defaultDirectionOf0;
					directionOf1=defaultDirectionOf1;
					DNAregion0 = defaultDNAregion;
					DNAregion1 = defaultDNAregion; 
					dimerisationProbability = defaultDimerisationProbability;
					affinityIncrease = Constants.NONE;
					isReversible = true;
					isAdditive = true;
					isFixed=false;
					
					rowID++;
					
					//species 0
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[0])); 
					species0ID = n.getTFspeciesID(cellContent);

					//species 1
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[1])); 
					species1ID = n.getTFspeciesID(cellContent);
					
					//type
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[2])); 
					type = Utils.parseInteger(cellContent, type);
					
					//sideOf0
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[3])); 
					directionOf0 = Utils.parseInteger(cellContent, directionOf0);
					
					//sideOf1
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[4])); 
					directionOf1 = Utils.parseInteger(cellContent, directionOf1);
							
					
					//DNA REGION 0
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[5])); 
					DNAregion0 = new DNAregion(cellContent, "", Constants.NONE, Constants.NONE,false,false);
					if(DNAregion0.start==Constants.NONE || DNAregion0.end==Constants.NONE ){
						n.printDebugInfo("warning when parsing the TF cooperativity csv file "+filename+": default value for DNA region 0 was used." );
						DNAregion0 = defaultDNAregion;
					}
					
					//DNA REGION 1
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[6])); 
					DNAregion1 = new DNAregion(cellContent, "", Constants.NONE, Constants.NONE,false,false);
					if(DNAregion1.start==Constants.NONE || DNAregion1.end==Constants.NONE ){
						n.printDebugInfo("warning when parsing the TF cooperativity csv file "+filename+": default value for DNA region 1 was used." );
						DNAregion1 = defaultDNAregion;
					}
					
					
					//dimmerisationRate
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[7])); 
					dimerisationProbability = Utils.parseDouble(cellContent, dimerisationProbability);
					
					
					//undimmerisationRate
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[8])); 
					affinityIncrease = Utils.parseDouble(cellContent, affinityIncrease);
						
					
					//undimmerisationRate
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[9])); 
					isReversible = Utils.parseBoolean(cellContent, isReversible);
								
					//undimmerisationRate
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[10])); 
					isAdditive = Utils.parseBoolean(cellContent, isAdditive);
					
					//undimmerisationRate
					cellContent =buffer.get(csv.header.get(Constants.PARSER_TF_COOPERATIVITY_CSV_FILE_HEADER[11])); 
					isFixed = Utils.parseBoolean(cellContent, isFixed);
					
					
					if(species0ID==Constants.NONE){
						n.printDebugInfo("Error on loading row "+rowID+" of TF cooperativity file. Species 0 is undefined.");
					} else if(species1ID==Constants.NONE){
						n.printDebugInfo("Error on loading row "+rowID+" of TF cooperativity file. Species 1 is undefined.");
					} else if(type==Constants.NONE){
						n.printDebugInfo("Error on loading row "+rowID+" of TF cooperativity file. Cooperativity type is undefined.");
					} else if(type == Constants.TF_COOPERATIVITY_TYPE_DNA && !DNAregion0.isSpecific(defaultDNAregion) && !DNAregion1.isSpecific(defaultDNAregion) && region.includes(DNAregion0) && region.includes(DNAregion1)){ 
						n.printDebugInfo("Error on loading row "+rowID+" of TF cooperativity file. For DNA based cooperativity specific DNA regions need to be defined.");
					} else if(type == Constants.TF_COOPERATIVITY_TYPE_DNA && affinityIncrease <= 0){ 
						n.printDebugInfo("Error on loading row "+rowID+" of TF cooperativity file. For DNA based cooperativity you need to specify the affinity increase.");
					//} else if(type == Constants.TF_COOPERATIVITY_TYPE_DIRECT &&( (dimmerisationRate <= 0 || undimmerisationRate<=0)  )){
					} else if(type == Constants.TF_COOPERATIVITY_TYPE_DIRECT &&( dimerisationProbability <= 0 || affinityIncrease <= 0)){ 
						n.printDebugInfo("Error on loading row "+rowID+" of TF cooperativity file. For direct TF-TF cooperativity you need to specify the dimerisation and undimeristion rate.");
					}  else if(type == Constants.TF_COOPERATIVITY_TYPE_DNA && (affinityIncrease==1)){ 
						n.printDebugInfo("Row "+rowID+" of TF cooperativity file was not load. An affinity increase of 1 means no cooperativity.");
					} else{
						
						 TFcooperativity bufferCoop;
						 if(type==Constants.TF_COOPERATIVITY_TYPE_DIRECT){
							 if(isAdditive){
								 if(isFixed){
									 bufferCoop = new TFcooperativityDirectAdditiveFixed(species0ID, species1ID, type, directionOf0, directionOf1, DNAregion0, DNAregion1, dimerisationProbability, affinityIncrease, isReversible, isAdditive, isFixed);
								 } else{
									 bufferCoop = new TFcooperativityDirectAdditive(species0ID, species1ID, type, directionOf0, directionOf1, DNAregion0, DNAregion1, dimerisationProbability, affinityIncrease, isReversible, isAdditive, isFixed);

								 }
							 } else{
								 bufferCoop = new TFcooperativityDirectMultiplicative(species0ID, species1ID, type, directionOf0, directionOf1, DNAregion0, DNAregion1, dimerisationProbability, affinityIncrease, isReversible, isAdditive, isFixed);
							 }
						 } else{
							 bufferCoop = new TFcooperativity(species0ID, species1ID, type, directionOf0, directionOf1, DNAregion0, DNAregion1, dimerisationProbability, affinityIncrease, isReversible, isAdditive, isFixed);
						 }
						
						data.add(bufferCoop);
					}
				}
			} 
			
		} else{
			n.printDebugInfo("could not load TF cooperativity csv file!");
		}
		
	}
}
