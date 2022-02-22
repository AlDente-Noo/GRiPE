package objects;

import java.io.Serializable;

import utils.CellUtils;

public class DNAsequence implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1329373046613469228L;
	public byte[] seq;
	
	/**
	 * class constructor
	 * @param seq
	 */
	public DNAsequence(byte[] seq){
		this.seq = new byte[seq.length];
		for(int i=0; i< seq.length;i++){
			this.seq[i]=seq[i];
		}
	}
	
	/**
	 * class constructor
	 * @param seq
	 */
	public DNAsequence(String seqStr){
		byte[] seq= CellUtils.getSeqIDs(seqStr);
		this.seq = new byte[seq.length];
		for(int i=0; i< seq.length;i++){
			this.seq[i]=seq[i];
		}
	}
	
	
	
	/**
	 * override equals to make this a custom class for hasmap keys
	 * @param other the object to compare with
	 * @return
	 */
	 @Override public boolean equals(Object other){
	        if (!(other instanceof DNAsequence)) {
	            return false;
	        }
	        DNAsequence seq = (DNAsequence) other;
	        
	       /// System.out.println(CellUtils.sequenceToString(this.seq) +" ? "+ CellUtils.sequenceToString(seq.seq)+" = "+CellUtils.areSequencesEqual(this.seq, seq.seq));
	        
	        return CellUtils.areSequencesEqual(this.seq, seq.seq);
	 }
		
	
	
		/**
		 * override hashCode to make this a custom class for hasmap keys
		 * @param seq
		 * @return
		 */
	  @Override public int hashCode() {
	        int hash = 0;
	        hash=CellUtils.sequenceToString(this.seq).hashCode();
	        return hash;
	    }

	
}
