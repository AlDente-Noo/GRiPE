package utils;

public class RepressionData {
    private final double time;
    private final int repressedLength;
    private final double repressedRepScore;
    private final double repressedActScore;

    public RepressionData(double time, int repressedLength, double repressedRepScore, double repressedActScore) {
        this.time = time;
        this.repressedLength = repressedLength;
        this.repressedRepScore = repressedRepScore;
        this.repressedActScore = repressedActScore;
    }

    public double getTime() {
        return time;
    }

    public int getRepressedLength() {
        return repressedLength;
    }

    public double getRepressedRepScore() {
        return repressedRepScore;
    }

    public double getRepressedActScore() {
        return repressedActScore;
    }
}
