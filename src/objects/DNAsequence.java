package objects;

import utils.CellUtils;

import java.io.Serializable;

public class DNAsequence implements Serializable {

    private static final long serialVersionUID = 1329373046613469228L;
    public byte[] seq;

    /**
     * class constructor
     */
    public DNAsequence(byte[] seq) {
        this.seq = new byte[seq.length];
		System.arraycopy(seq, 0, this.seq, 0, seq.length);
    }

    /**
     * class constructor
     */
    public DNAsequence(String seqStr) {
        byte[] seq = CellUtils.getSeqIDs(seqStr);
        this.seq = new byte[seq.length];
		System.arraycopy(seq, 0, this.seq, 0, seq.length);
    }


    /**
     * override equals to make this a custom class for hasmap keys
     *
     * @param other the object to compare with
     */
    @Override
    public boolean equals(Object other) {
        if (!(other instanceof DNAsequence)) {
            return false;
        }
        DNAsequence seq = (DNAsequence) other;
        return CellUtils.areSequencesEqual(this.seq, seq.seq);
    }

    /**
     * override hashCode to make this a custom class for hashmap keys
     */
    @Override
    public int hashCode() {
        return CellUtils.sequenceToString(this.seq).hashCode();
    }

}
