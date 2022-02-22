package objects;


/**
 * class that extends TF cooperativity to be multiplicative 
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFcooperativityDirectMultiplicative  extends TFcooperativity{

	/**
	 * 
	 */
	private static final long serialVersionUID = 237093256997502380L;

	/**
	 * class constructor
	 * @param species0id
	 * @param species1id
	 * @param type
	 * @param direction0
	 * @param direction1
	 * @param region0
	 * @param region1
	 * @param dimerisationProbability
	 * @param affinityIncrease
	 * @param isReversible
	 */
	public TFcooperativityDirectMultiplicative(int species0id, int species1id,
			int type, int direction0, int direction1, DNAregion region0,
			DNAregion region1, double dimerisationProbability,
			double affinityIncrease, boolean isReversible, boolean isAdditive, boolean isFixed) {
		super(species0id, species1id, type, direction0, direction1, region0, region1,
				dimerisationProbability, affinityIncrease, isReversible, isAdditive, isFixed);
	}

	/**
	 * class constructor
	 * @param coop
	 */
	public TFcooperativityDirectMultiplicative( TFcooperativity coop){
		super (coop);
	}
	
	/**
	 * method that computes 
	 */
	public double computeMoveRate(double moveRate0, double moveRate1, double factor){
		return moveRate0/factor;
	}
	
}
