package objects;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import environment.Cell;
import utils.CellUtils;
import utils.Constants;
import utils.Utils;

/**
 * class that imports a Position Frequency Matrix 
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class PFM  implements Serializable{

    private static final long serialVersionUID = 1372756026390102207L;
    public ArrayList<ArrayList<Double>> pfm;
    public ArrayList<ArrayList<Double>> normPFM;

    public int[] nucleotidePosition;
    public int motifSize;
    public boolean isCorrect;

    /**
     * class constructor
     * an example of a correct PFM specification is A=[3, 1, 5, 7, 3, 6, 4, 7, 1, 9, 8, 5, 4, 2]; C=[1, 1, 2, 0, 0, 0, 1, 0, 8, 0, 0, 0, 0, 3]; G=[4, 1, 1, 1, 0, 0, 4, 1, 0, 0, 1, 3, 0, 2]; T=[1, 6, 1, 1, 6, 3, 0, 1, 0, 0, 0, 1, 5, 2]
     * @param str string with PFM
     */
    public PFM(String str, Cell n){
        if(str.startsWith(Constants.DBD_TYPE_PWM)){
            str = str.replaceAll(Constants.DBD_TYPE_PWM, "").trim();
        }
        else {
            n.stopSimulation("PFM (PWM) string has to begin with '"+Constants.DBD_TYPE_PWM+"'");
        }
        pfm = new ArrayList<ArrayList<Double>>();
        motifSize = Constants.NONE;
        ArrayList<Double> bufferPFM;
        isCorrect = true;
        nucleotidePosition = new int[CellUtils.bps.numberOfBP];
        Arrays.fill(nucleotidePosition, Constants.NONE);

        if(str.contains(Constants.PFM_NUCLEOTIDE_SEPARATOR)){
            String[] bufferNucleotide, bufferNucleotideContainer;
            int nucleotideID, matrixID;
            String nucleotideKey, nucleotidePFM;

            bufferNucleotide = str.split(Constants.PFM_NUCLEOTIDE_SEPARATOR);

            for(String buffer:bufferNucleotide){
                if(buffer.contains(Constants.PFM_NUCLEOTIDE_ASSIGNMENT)){
                    bufferNucleotideContainer = buffer.split(Constants.PFM_NUCLEOTIDE_ASSIGNMENT);
                    if(bufferNucleotideContainer.length == 2){
                        nucleotideKey = bufferNucleotideContainer[0].trim();
                        if(nucleotideKey.length()==1 && CellUtils.bps.bpsID.containsKey(nucleotideKey)){
                            nucleotideID = CellUtils.bps.bpsID.get(nucleotideKey);

                            nucleotidePFM = bufferNucleotideContainer[1].trim();
                            nucleotidePFM = nucleotidePFM.substring(1);
                            nucleotidePFM = nucleotidePFM.substring(0,nucleotidePFM.length()-1);

                            if(nucleotidePFM.contains(Constants.PFM_NUCLEOTIDE_CONTAINER_SEPRATOR)){
                                bufferPFM = new ArrayList<Double>();
                                for(String bufferNucleotidePosPFM:nucleotidePFM.split(Constants.PFM_NUCLEOTIDE_CONTAINER_SEPRATOR)){
                                    bufferPFM.add(Utils.parseDouble(bufferNucleotidePosPFM.trim(), 0));
                                }
                                if(motifSize==Constants.NONE){
                                    motifSize=bufferPFM.size();
                                }

                                if(bufferPFM.size()!=motifSize){
                                    isCorrect =false;
                                    n.printDebugInfo("pfm " +str+ " various size rows pfm");
                                }

                                pfm.add(bufferPFM);
                                matrixID = pfm.size()-1;
                                this.nucleotidePosition[nucleotideID] = matrixID;
                            }
                        }
                    }

                }
            }
        }
        else {
            n.stopSimulation("nucleotides in PFM (PWM) string has to be separated with '"+Constants.PFM_NUCLEOTIDE_SEPARATOR+"'");
        }

        //not all fill with zero
        if(isCorrect && pfm.size()!=CellUtils.bps.numberOfBP){
            int id;
            for(String buffer:CellUtils.bps.bps){
                if(CellUtils.bps.bpsID.containsKey(buffer)){
                    id = CellUtils.bps.bpsID.get(buffer);
                    if(id!=CellUtils.bps.getANYID() && nucleotidePosition[id]==Constants.NONE){
                        bufferPFM = new ArrayList<Double>();
                        for(int i=0;i<motifSize;i++){
                            bufferPFM.add(0.0);
                        }
                        pfm.add(bufferPFM);
                        nucleotidePosition[id] =pfm.size()-1;
                    }
                }
            }
        }

        this.normPFM = this.pfm;
    }


    /**
     * FG
     * returns the PWM-score of the nucleotide
     * @param nucleotide current nucleotide
     * @param position in the motif
     */
    public double getScorePFM(byte nucleotide, int position){
        return normPFM.get(this.nucleotidePosition[nucleotide]).get(position);
    }

    /**
     * AD
     * returns max PWM score for a position on the DNA
     * FG: this function is needed to normalize PWM-scores, so that the strongest site has PWM-score = 0
     * and weaker sites have negative scores. The affinity of a site is
     */
    public double getMaxScorePFM(int position){
        double max = 0;
        for (int value : this.nucleotidePosition) {
            if (normPFM.get(value).get(position) > max) {
                max = normPFM.get(value).get(position);
            }
        }
        return max;
    }


    /**
     * generates a string description of the class
     * @return
     */
    public String toString(double[] bpFreq){
        StringBuffer str= new StringBuffer("");
        str.append(Constants.DBD_TYPE_PWM);
        for(int i=0;i<this.pfm.size();i++){
            str.append(CellUtils.bps.bps[i]);
            str.append("=[");
            for (int j=0;j<motifSize; j++){
                str.append(this.normPFM.get(i).get(j));
                if(j<motifSize-1){
                    str.append(", ");
                }
            }
            if(i<this.pfm.size()-1){
                str.append("]; ");
            } else{
                str.append("]");
            }
        }
        return str.toString();
    }


}
