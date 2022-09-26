package utils;

import objects.BasePairs;
import objects.DNAsequence;
import objects.PFM;
import objects.TargetSitesGroupLogic;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
 * this class contains static fields and methods that are used for the biology part of the application
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
public class CellUtils {

    public static final BasePairs bps = new BasePairs();
    public static final TargetSitesGroupLogic tgsl = new TargetSitesGroupLogic();

    public static final byte bpANYID = bps.getANYID();

    /**
     * generate random sequence of nucleotides ( can be used for DBD of TFs or DNA sequence the probability of
     * nucleotides are uniformly distributed
     *
     * @param generator the random number generator
     * @param length    the length in nucleotides
     * @return a string of bytes represnting the nucleotide sequence
     */
    public static byte[] generateRandomDNASequence(Random generator, int length, double proportionOfA,
                                                   double proportionOfT, double proportionOfC, double proportionOfG) {
        byte[] seq = new byte[length];
        double buffer;

        double quotaOfA = proportionOfA;
        double quotaOfT = quotaOfA + proportionOfT;
        double quotaOfC = quotaOfT + proportionOfC;
        double quotaOfG = quotaOfC + proportionOfG;

        byte current;

        for (int i = 0; i < length; i++) {
            buffer = Utils.generateNextDouble(generator, 0, quotaOfG);

            if (buffer < quotaOfA) {
                current = bps.bpsID.get("A");
            } else if (buffer < quotaOfT) {
                current = bps.bpsID.get("T");
            } else if (buffer < quotaOfC) {
                current = bps.bpsID.get("C");
            } else {
                current = bps.bpsID.get("G");
            }

            seq[i] = current;
        }

        return seq;
    }


    /**
     * computes the affinities between a TF and DNA
     *
     * @param strand    the DNA strand
     * @param TFseq     the recognise DNA sequence
     * @param sizeLeft  The size on the left of the DBD that the TF occupies on the DNA.
     * @param sizeTotal The total number of bp  that the TF occupies on the DNA.
     * @param es        the specific energy
     * @return the affinity vector
     */
    public static double[] computeTFAffinities(Random generator, byte[] strand, byte[] TFseq, int sizeLeft,
                                               int sizeTotal, double es, int direction, double roughness) {
        double[] affinities = new double[strand.length];
        sizeTotal--;

        if (TFseq != null && TFseq.length > 0) {
            if (direction == 0) {
                for (int i = 0; i < strand.length - sizeTotal; i++) {
                    affinities[i] = computeTFAffinityLR(strand, i, TFseq, sizeLeft, es);
                }
            }
            if (direction == 1) {
                for (int i = 0; i < strand.length - sizeTotal; i++) {
                    affinities[i] = computeTFAffinityRL(strand, i, TFseq, sizeLeft, es);
                }
            }
        } else {
            for (int i = 0; i < strand.length - sizeTotal; i++) {
                affinities[i] = Utils.generateNextNormalDistributedDouble(generator, es, roughness, 0);
            }
        }

        for (int i = strand.length - sizeTotal; i < strand.length; i++) {
            affinities[i] = Constants.NONE;
        }

        return affinities;
    }


    /**
     * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 two ways from 5'
     * to 3' and from 3' to 5'
     *
     * @param DNAseq   the DNA sequence
     * @param DNApos   the index of the starting bp of the position on the DNA where the affinity is computed
     * @param TFseq    the recognise DNA sequence
     * @param sizeLeft TF size to the left of its motif
     * @param es       the specific energy
     * @return the value of the affinity
     */
    public static double computeTFAffinity2Way(byte[] DNAseq, int DNApos, byte[] TFseq, int sizeLeft, double es) {
        double sumLR = 0;
        double sumRL = 0;

        byte[] revSeq = getReversedComplementSequences(TFseq);

        for (int i = 0; i < TFseq.length; i++) {
            if (TFseq[i] != bpANYID && DNAseq[DNApos + i + sizeLeft] != TFseq[i]) {
                sumLR += es;
            }

            if (revSeq[i] != bpANYID && DNAseq[DNApos + i + sizeLeft] != revSeq[i]) {
                sumRL += es;
            }
        }
        return Math.min(sumLR, sumRL);
    }


    /**
     * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 one way from 5'
     * to 3'
     *
     * @param DNAseq   the DNA sequence
     * @param DNApos   the index of the starting bp of the position on the DNA where the affinity is computed
     * @param TFseq    the recognise DNA sequence
     * @param sizeLeft TF size to the left of its motif
     * @param es       the specific energy (mismatch energy in kT units)
     * @return the value of the affinity
     */
    public static double computeTFAffinityLR(byte[] DNAseq, int DNApos, byte[] TFseq, int sizeLeft, double es) {
        double sumLR = 0;
        for (int i = 0; i < TFseq.length; i++) {
            if (TFseq[i] != bpANYID && DNAseq[DNApos + i + sizeLeft] != TFseq[i]) {
                sumLR += es;
            }
        }
        return sumLR;
    }


    /**
     * computes the affinity between a TF and the DNA at a specific position using the Gerland 2002 two ways from 3'
     * to 5'
     *
     * @param DNAseq   the DNA sequence
     * @param DNApos   the index of the starting bp of the position on the DNA where the affinity is computed
     * @param TFseq    the recognise DNA sequence
     * @param sizeLeft TF size to the left of its motif
     * @param es       the specific energy (mismatch energy in kT units)
     * @return the value of the affinity
     */
    public static double computeTFAffinityRL(byte[] DNAseq, int DNApos, byte[] TFseq, int sizeLeft, double es) {
        double sumRL = 0;
        byte[] revSeq = getReversedComplementSequences(TFseq);
        for (int i = 0; i < TFseq.length; i++) {
            if (revSeq[i] != bpANYID && DNAseq[DNApos + i + sizeLeft] != revSeq[i]) {
                sumRL += es; // TODO: change to -es ??
            }
        }
        return sumRL;
    }


    /**
     * computes the affinity between a TF and the DNA at a specific position using a known sequence affinity 3' -> 5'
     *
     * @param DNAseq          the DNA sequence
     * @param DNApos          the index of the starting bp of the position on the DNA where the affinity is computed
     * @param seqsAffinities  the list of sequences and their affinties
     * @param sizeLeft        the left size
     * @param sizeMotif       the size of the motif
     * @param defaultAffinity the default affinity if not found in the list
     * @return the value of the affinity
     */
    public static double computeTFAffinityLR(byte[] DNAseq, int DNApos, HashMap<DNAsequence, Double> seqsAffinities,
                                             int sizeLeft, int sizeMotif, double defaultAffinity) {
        double sumLR = defaultAffinity;

        byte[] seq = new byte[sizeMotif - sizeLeft];

        for (int i = 0; i < sizeMotif; i++) {
            seq[i] = DNAseq[i + DNApos + sizeLeft];
        }

        DNAsequence buffer = new DNAsequence(seq);

        if (seqsAffinities.containsKey(buffer)) {
            sumLR = seqsAffinities.get(buffer);
        }

        return sumLR;
    }


    /**
     * computes the affinity between a TF and the DNA at a specific position using a known sequence affinity 5' -> 3'
     *
     * @param DNAseq          the DNA sequence
     * @param DNApos          the index of the starting bp of the position on the DNA where the affinity is computed
     * @param seqsAffinities  the list of sequences and their affinties
     * @param sizeLeft        the left size
     * @param sizeMotif       the size of the motif
     * @param defaultAffinity the default affinity if not found in the list
     * @return the value of the affinity
     */
    public static double computeTFAffinityRL(byte[] DNAseq, int DNApos, HashMap<DNAsequence, Double> seqsAffinities,
                                             int sizeLeft, int sizeMotif, double defaultAffinity) {
        double sumLR = defaultAffinity;

        byte[] seq = new byte[sizeMotif - sizeLeft];

        for (int i = 0; i < sizeMotif; i++) {
            seq[i] = DNAseq[DNApos + sizeLeft + (sizeMotif - 1) - i];
        }

        DNAsequence buffer = new DNAsequence(seq);

        if (seqsAffinities.containsKey(buffer)) {
            sumLR = seqsAffinities.get(buffer);
        }

        return sumLR;
    }


    /**
     * computes the affinities between a TF and DNA (using PWM or randomly)
     *
     * @param generator random number generator
     * @param strand    the DNA strand
     * @param pfm       PFM (PWM)
     * @param sizeLeft  The size on the left of the DBD that the TF occupies on the DNA.
     * @param sizeTotal The total number of bp  that the TF occupies on the DNA.
     * @param es        the specific energy
     * @param direction the direction (5'->3' or 3'->5')
     * @param roughness the affinity landscape roughness
     * @return the affinity vector
     */
    public static double[] computeTFAffinities(Random generator, byte[] strand, PFM pfm, int sizeLeft, int sizeTotal,
                                               double es, int direction, double roughness) {
        double[] affinities = new double[strand.length];
        sizeTotal--;

        if (pfm != null && pfm.isCorrect && pfm.motifSize > 0) {
            if (direction == 0) {
                for (int i = 0; i < strand.length - sizeTotal; i++) {
                    affinities[i] = computeTFAffinityLR(strand, i, pfm, sizeLeft, es, false);
                }
            }
            if (direction == 1) {
                for (int i = 0; i < strand.length - sizeTotal; i++) {
                    affinities[i] = computeTFAffinityRL(strand, i, pfm, sizeLeft, es, false);
                }
            }
        } else {
            for (int i = 0; i < strand.length - sizeTotal; i++) {
                affinities[i] = Utils.generateNextNormalDistributedDouble(generator, es, roughness, 0);
            }
        }

        for (int i = strand.length - sizeTotal; i < strand.length; i++) {
            affinities[i] = Constants.NONE;
        }

        return affinities;
    }


    /**
     * computes the affinities between a TF and DNA from user-specified affinities
     *
     * @param strand    the DNA strand
     * @param sizeLeft  The size on the left of the DBD that the TF occupies on the DNA.
     * @param sizeTotal The total number of bp  that the TF occupies on the DNA.
     * @return the affinity vector
     */
    public static double[] computeTFAffinities(byte[] strand, HashMap<DNAsequence, Double> seqsAffinities,
                                               double defaultAffinity, int sizeLeft, int sizeMotif, int sizeTotal,
                                               int direction) {
        double[] affinities = new double[strand.length];
        sizeTotal--;

        if (seqsAffinities != null && seqsAffinities.size() > 0) {
            if (direction == 0) {
                for (int i = 0; i < strand.length - sizeTotal; i++) {
                    affinities[i] = computeTFAffinityLR(strand, i, seqsAffinities, sizeLeft, sizeMotif,
                            defaultAffinity);
                }
            }
            if (direction == 1) {
                for (int i = 0; i < strand.length - sizeTotal; i++) {
                    affinities[i] = computeTFAffinityRL(strand, i, seqsAffinities, sizeLeft, sizeMotif,
                            defaultAffinity);
                }
            }
        } else {
            for (int i = 0; i < strand.length - sizeTotal; i++) {
                affinities[i] = defaultAffinity;
            }
        }

        for (int i = strand.length - sizeTotal; i < strand.length; i++) {
            affinities[i] = Constants.NONE;
        }

        return affinities;
    }

    /**
     * computes the PWM (PFM) affinity between a TF and the DNA at a specific position for the 5'->3' strand
     *
     * @param DNAseq   the DNA sequence
     * @param DNApos   the index of the starting bp of the position on the DNA where the affinity is computed
     * @param pfm      PFM
     * @param sizeLeft TF size to the left of its motif
     * @param es       the specific energy (converts PWM-score to energy in kT units: E = -es * PWM-score)
     * @param isPWMscore if true returns the score, otherwise returns normalized score
     * @return the value of the affinity
     */
    public static double computeTFAffinityLR(byte[] DNAseq, int DNApos, PFM pfm, int sizeLeft, double es,
                                             boolean isPWMscore) {
        double sumLR = 0, sumMax = 0;
        for (int i = 0; i < pfm.motifSize; i++) {
            sumLR += es * pfm.getScorePFM(DNAseq[DNApos + i + sizeLeft], i);
            sumMax += es * pfm.getMaxScorePFM(i);
        }
        return isPWMscore ? sumLR : (sumLR - sumMax);
    }


    /**
     * computes the PWM (PFM) affinity between a TF and the DNA at a specific position for the 3'->5' strand
     *
     * @param DNAseq   the DNA sequence
     * @param DNApos   the index of the starting bp of the position on the DNA where the affinity is computed
     * @param pfm      PFM
     * @param sizeLeft TF size to the left of its motif
     * @param es       the specific energy (converts PWM-score to energy in kT units)
     * @param isPWMscore if true returns the score, otherwise returns normalized score
     * @return the value of the affinity
     */
    public static double computeTFAffinityRL(byte[] DNAseq, int DNApos, PFM pfm, int sizeLeft, double es,
                                             boolean isPWMscore) {
        double sumRL = 0, sumMax = 0;
        byte[] revComplement = getReversedComplementSequences(DNAseq, DNApos + sizeLeft, pfm.motifSize);
        for (int i = 0; i < pfm.motifSize; i++) {
            sumRL += es * pfm.getScorePFM(revComplement[i], i);
            sumMax += es * pfm.getMaxScorePFM(i);
        }
        return isPWMscore ? sumRL : (sumRL - sumMax);
    }


    /**
     * computes the rate at which a bound TF will take a decision whether to move or not from current position.
     */
    public static double computeAvgMoveRate(double specificWaitingTime, double bindingEnergy) {
        return 1.0 / (specificWaitingTime * Math.exp(-bindingEnergy));
    }

    /**
     * converts waiting time into binding energy
     */
    public static double computeBindingEnergy(double specificWaitingTime, double TFwaitingTime) {
        // FG: TFwaitingTime = specificWaitingTime * exp(-bindingEnergy)
        return -Math.log(TFwaitingTime / specificWaitingTime);
    }


    /**
     * converts a string of DNA sequence into the corresponding vector of bytes
     */
    public static ArrayList<Byte> getSequenceIDs(String strand) {

        ArrayList<Byte> buffer = new ArrayList<Byte>();

        strand = strand.trim();
        strand = strand.toUpperCase();

        String letter;

        for (int i = 0; i < strand.length(); i++) {
            letter = strand.substring(i, i + 1);
            if (bps.bpsID.containsKey(letter)) {
                buffer.add(bps.bpsID.get(letter));
            }
        }

        return buffer;
    }

    /**
     * converts a string of DNA sequence into the corresponding vector of bytes
     */
    public static byte[] getSeqIDs(String strand) {

        byte[] buffer = new byte[strand.length()];

        strand = strand.trim();
        strand = strand.toUpperCase();

        String letter;

        for (int i = 0; i < strand.length(); i++) {
            buffer[i] = Constants.NONE;
            letter = strand.substring(i, i + 1);
            if (bps.bpsID.containsKey(letter)) {
                buffer[i] = bps.bpsID.get(letter);
            }
        }

        return buffer;
    }

    /**
     * returns the literal string of the byte specified sequence
     */
    public static String sequenceToString(byte[] seq) {
        StringBuilder result = new StringBuilder();
        if (seq != null && seq.length > 0) {
            for (byte b : seq) {
                if (b >= 0 && b < bps.bps.length) {
                    result.append(bps.bps[b]);
                }
            }
        }
        return result.toString();
    }


    /**
     * returns the reversed complement of a sequence
     */
    public static ArrayList<Byte> getReversedComplementSequences(ArrayList<Byte> seq) {
        ArrayList<Byte> revSeq = new ArrayList<Byte>();
        for (int i = seq.size() - 1; i >= 0; i--) {
            revSeq.add(bps.getComplement(seq.get(i)));
        }
        return revSeq;
    }


    /**
     * returns the reversed complement of a sequence
     */
    public static byte[] getReversedComplementSequences(byte[] seq) {
        byte[] revSeq = new byte[seq.length];
        int j = 0;
        for (int i = seq.length - 1; i >= 0; i--) {
            revSeq[j] = bps.getComplement(seq[i]);
            j++;
        }
        return revSeq;
    }


    /**
     * returns the reversed complement of a subsequence
     *
     * @param seq    the DNA sequence
     * @param start  the start position
     * @param length the length of the sequence
     */
    public static byte[] getReversedComplementSequences(byte[] seq, int start, int length) {
        if (length <= 0 || ((start + length) >= seq.length)) {
            length = seq.length;
            start = 0;
        }
        byte[] revSeq = new byte[length];
        int j = 0;
        for (int i = (length + start) - 1; i >= start; i--) {
            revSeq[j] = bps.getComplement(seq[i]);
            j++;
        }
        return revSeq;
    }


    /**
     * compare two sequences
     */
    public static boolean areSequencesEqual(byte[] seq1, byte[] seq2) {
        boolean result = false;
        if (seq1.length == seq2.length && seq1.length > 0) {
            result = true;
            for (int i = 0; i < seq1.length && result; i++) {
                if (seq1[i] != seq2[i] && seq1[i] != bpANYID && seq2[i] != bpANYID) {
                    result = false;
                }
            }
        }
        return result;
    }


    /**
     * print a DNA sequence to a file
     *
     * @param path        the path were to save the file
     * @param filename    the filename
     * @param description the description
     * @param seq         the DNA sequence
     */
    public static void printSequence(String path, String filename, String description, byte[] seq) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(path, filename)));
            out.write(">");
            out.write(description);
            out.newLine();
            int sectorSize = 100;
            int sectors = (int) Math.ceil((double) seq.length / sectorSize);
            StringBuilder buffer;
            for (int i = 0; i < sectors; i++) {
                buffer = new StringBuilder();
                for (int j = i * sectorSize; j < Math.min((i + 1) * sectorSize, seq.length); j++) {
                    buffer.append(CellUtils.bps.bps[seq[j]]);
                }
                out.write(buffer.toString());
                out.newLine();
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * concatenates two sequences and returns the result
     */
    public static byte[] concatenateDNAseq(byte[] seq1, byte[] seq2) {
        byte[] result;
        result = new byte[seq1.length + seq2.length];
        System.arraycopy(seq1, 0, result, 0, seq1.length);
        System.arraycopy(seq2, 0, result, seq1.length, seq2.length);
        return result;
    }


    /**
     * creates a copy of a DNA sequence and returns it
     */
    public static byte[] copySequence(byte[] seq) {
        byte[] result = null;
        if (seq != null && seq.length > 0) {
            result = new byte[seq.length];
            System.arraycopy(seq, 0, result, 0, seq.length);
        }
        return result;
    }

    /**
     * replaces the DNA sequence in the strand from position pos with the subsequence subseq
     */
    public static byte[] replaceDNAseq(byte[] strand, byte[] subSeq, int pos) {
        byte[] result = copySequence(strand);
        if (pos > 0 && pos < strand.length) {
            for (int i = pos; i < Math.min(pos + subSeq.length, strand.length); i++) {
                if (subSeq[i - pos] != bpANYID) {
                    result[i] = subSeq[i - pos];
                }
            }
        }
        return result;
    }


    /**
     * generates a DNA strand
     *
     * @param length length of the generated strand
     */
    public static byte[] generateEmptyDNAStrand(int length) {
        byte[] result = new byte[length];
        for (int i = 0; i < length; i++) {
            result[i] = bpANYID;
        }
        return result;
    }

    /**
     * generates a DNA seq that contains only one bp
     *
     * @param length length of the generated strand
     * @param bp     the only base pair of the strand
     */
    public static byte[] generateDNAStrand(int length, String bp) {
        byte[] result = new byte[length];
        byte bpID = bps.bpsID.get(bp);
        for (int i = 0; i < length; i++) {
            result[i] = bpID;
        }

        return result;
    }

    /**
     * extracts a subsequence form a sequence
     *
     * @param seq   the original sequence
     * @param start the start
     * @param end   the end (exclusive)
     * @return the subsequence
     */
    public static byte[] extractSubSequence(byte[] seq, int start, int end) {
        int s = Math.max(0, start);
        int e = Math.min(seq.length, end);
        byte[] result = null;
        if (e > s) {
            result = new byte[e - s];
            System.arraycopy(seq, s, result, 0, e - s);
        }
        return result;
    }
}
