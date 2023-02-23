package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Random;

/**
 * a class with useful general methods to write into a file or generate random numbers
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class Utils {

	/**
	 * generates a random double number between two values. when the two values are equal it returns default value
	 * @param generator random generator
	 * @param min inclusive
	 * @param max exclusive
	 * @return random integer from between min and max
	 */
	public static double generateNextDouble(Random generator, double min, double max){
		double result=min;
		if(min<max){
			result = min + generator.nextDouble()*(max - min);
		}
		return result;
	}
	
	/**
	 * generates a random integer number between two values. when the two values are equal it returns default value
	 * @param generator random generator
	 * @param min inclusive
	 * @param max exclusive
	 * @return random integer from between min and max
	 */
	public static int generateNextInteger(Random generator, int min, int max){
		int result=min;
		if(min<max){
			result = min+generator.nextInt(max-min);
		}
		return result;
	}


	/**
	 * returns a normal distributed integer with a specific mean and  standard deviation
	 */
	public static int generateNextNormalDistributedInteger(Random generator, double mean, double stddev){
		return  (int) Math.round((generator.nextGaussian()*stddev+mean));
	}
	

	/**
	 * returns a normal distributed double with a specific mean and stddev and min value
	 */
	public static double generateNextNormalDistributedDouble(Random generator, double mean, double stddev, double min){
		double result =   (generator.nextGaussian()*stddev+mean);
		return Math.max(result, min);
	}
	
	/**
	 * checks if a text is double and if so it gets the number
	 */
	public static double parseDouble(String str, double none){	
		double result = none;
		try{
			result = Double.parseDouble(str);
		} catch(NumberFormatException e){
			e.printStackTrace();
		}
		return result;
	}
	

	/**
	 * checks if a text is int and if so it gets the number
	 */
	public static int parseInteger(String str, int none){	
		int result = none;
		try{
			result = Integer.parseInt(str);
		} catch(NumberFormatException e){
			e.printStackTrace();
		}
		
		return result;
	}
	

	/**
	 * checks if a text is int and if so it gets the number
	 */
	public static long parseLong(String str, long none){	
		long result = none;
		try{
			result = Long.parseLong(str);
		} catch(NumberFormatException e){
			e.printStackTrace();
		}
		
		return result;
	}
	

	/**
	 * checks if a text is boolean and if so it gets the value
	 */
	public static boolean parseBoolean(String str, boolean none){	
		boolean result = none;
		try{
			result = Boolean.parseBoolean(str);
		} catch(Exception e){
			e.printStackTrace();
		}
		
		return result;
	}
	
	/**
	 * rounds a double to a 2 decimal number
	 * @param d the double number
	 * @return the 2 decimal double number
	 */
	public static double roundTwoDecimals(double d) {
		DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.getDefault());
		otherSymbols.setDecimalSeparator('.');
		DecimalFormat twoDForm = new DecimalFormat("#.##", otherSymbols);
		return Double.parseDouble(twoDForm.format(d));
	}
	
	
	/**
	 * computes the sum of a int vector
	 */
	public static long computeSum(int[] vector){
		long sum=0;
		for(int v:vector){
			sum+=v;
		}
		return sum;
	}


	/**
	 * computes the average value of a vector
	 */
	public static <N> double computeMean(ArrayList<N> vector){
		double avg = 0;
		for(N v: vector){
			avg += (double) v;
		}
		if(vector.size()>0){
			avg = avg / vector.size();
		}
		return avg;
	}

	/**
	 * breaks the command line into two strings the parameter and the value
	 * @param text the command line
	 * @return a 2 value array of strings 0-the parameter 1- the value
	 */
	public static String[] extractParameterFromCommandLine(String text, String assignmentChar){
		String[] result= new String[2];
		int firstPos, lastPos;
		result[0]="";
		result[1]="";
		text = text.trim();
		if(text.contains(assignmentChar)){
			firstPos = text.indexOf(assignmentChar);
			lastPos = text.lastIndexOf(assignmentChar);
			if(firstPos >=0 && firstPos == lastPos){
				result[0]=text.substring(0,firstPos);
				//remove spaces
				result[0] = result[0].trim();
				result[1]=text.substring(firstPos+1,text.length()-1);
				//remove spaces
				result[1] = result[1].trim();
				//remove "
				result[1] = result[1].replace("\"", "");
			} else{
				System.out.println("error in parameters line: \""+text+"\" (more than one equal signs)");
			}
		} else{
			System.out.println("error in parameters line: \""+text+"\" (no equal sign)");
		}

		
		return result;
	}
	
	
	/**
	 * return all lines from a file into an array list of strings
	 */
	public static ArrayList<String> readLinesFromFile(String filename) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(filename));
        ArrayList<String> lines = new  ArrayList<String>();
        String line;
	    try {
	        while ((line = is.readLine()) !=null) {
	        		if(!line.isEmpty()){
	        			lines.add(line);
	        		}
	        }
	    } finally {
	        is.close();
	    }
        return lines;
	}

	/**
	 * parse a string into an array of doubles
	 * @param line the string to be parsed
	 * @param delimiter the delimiter between the numbers
	 * @param defaultValue the default value to substitute an error in parse
	 * @return an array of doubles
	 */
	public static double[] parseCSVline(String line, String delimiter, double defaultValue){
		double[] result=null;
		String[] cells = line.split(delimiter);
		
		if(cells.length>0){
			result=new double[cells.length]; 
			for(int i=0;i<cells.length;i++){
				result[i] = parseDouble(cells[i], defaultValue);
			}
		}
			
        return result;
	}
}
