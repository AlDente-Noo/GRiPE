package objects;

/**
 * class that extends TF cooperativity class to be fixed additive 
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class TFcooperativityDirectAdditiveFixed   extends TFcooperativity{

	/**
	 * 
	 */
	private static final long serialVersionUID = -6043187560043870350L;

	public TFcooperativityDirectAdditiveFixed(int species0id, int species1id,
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
	public TFcooperativityDirectAdditiveFixed( TFcooperativity coop){
		super (coop);
	}
	
	/**
	 * method that computes 
	 */
	public double computeMoveRate(double moveRate0, double moveRate1, double factor){
		return (moveRate0 * factor)/(moveRate0 + factor);
	}
}
