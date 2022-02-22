package utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

/**
 * a class with useful general methods to write into a file or generate random numbers
 * @author n.r.zabet@gen.cam.ac.uk
 *
 */
public class Utils {

	
	private static boolean hasNextNextGaussian=false;
	private static double nextNextGaussian=0;

	
	

	
	/**
	 * generates a random double number between two values. when the two values are equal it returns default value
	 * @param generator
	 * @param min inclusive
	 * @param max exclusive
	 * @return
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
	 * @param generator
	 * @param min inclusive
	 * @param max exclusive
	 * @return
	 */
	public static int generateNextInteger(Random generator, int min, int max){
		int result=min;
		if(min<max){
			result = min+generator.nextInt(max-min);
		}
		return result;
	}


	/**
	 * returns a normal distributed integer with a specific mean and  stdandard deviation
	 * @param generator
	 * @param mean
	 * @param stddev
	 * @return
	 */
	public static int generateNextNormalDistributedInteger(Random generator, double mean, double stddev){

		return  (int) Math.round((generator.nextGaussian()*stddev+mean));
	}
	

	/**
	 * returns a normal distributed double with a specific mean and stddev and min value
	 * @param generator
	 * @param mean
	 * @param stddev
	 * @return
	 */
	public static double generateNextNormalDistributedDouble(Random generator, double mean, double stddev, double min){
		double result =   (generator.nextGaussian()*stddev+mean);
		
		return result>=min?result:min;
	}	

	/**
	 * returns a normal distributed double with a specific mean and variance
	 * @param generator
	 * @param mean
	 * @param variance
	 * @return
	 */
	/*public static double generateNextNormalDistributedDouble(Random generator, double mean, double variance, double min){
		double result = generator.nextGaussian()*variance+mean;
		return  result>=min?result:min;
	}*/
	
	
	
	 public static double generateNextGaussian(Random generator, double mean, double variance) {
		     double v1, v2, s;
		     if (hasNextNextGaussian) {
		    	 	hasNextNextGaussian = false;
		         return nextNextGaussian;
		       } else {
			     do {
			       v1 = 2 * generator.nextDouble() - 1;   // between -1.0 and 1.0
			       v2 = 2 * generator.nextDouble() - 1;   // between -1.0 and 1.0
			       s = v1 * v1 + v2 * v2;
			     } while (s >= 1 || s == 0);
			     double mult=Math.sqrt(-2 * Math.log(s) / s);
			     hasNextNextGaussian = true;
			     nextNextGaussian = mean+variance*v2*mult;
			     return mean+variance*v1 * mult;
		       }
		 }
	

		/**
		 * method to generate a Poisson random number
		 * @param generator the random number generator
		 * @param mean the mean
		 * @return
		 */
		public static int generateNextPoissonNumber(Random generator, double mean){
			double l = Math.exp(-mean);
			int k=0;
			double p=1;
			while (p>l){
				k++;
				p=p*generator.nextDouble();
			}
			return k-1;
		}

	
	
	/**
	 * checks if a text is double and if so it gets the number
	 * @param str
	 * @param none
	 * @return
	 */
	public static double parseDouble(String str, double none){	
		double result = none;
		try{
			result = Double.parseDouble(str);
		} catch(NumberFormatException e){
			
		}
		return result;
	}
	

	/**
	 * checks if a text is int and if so it gets the number
	 * @param str
	 * @param none
	 * @return
	 */
	public static int parseInteger(String str, int none){	
		int result = none;
		try{
			result = Integer.parseInt(str);
		} catch(NumberFormatException e){
			
		}
		
		return result;
	}
	

	/**
	 * checks if a text is int and if so it gets the number
	 * @param str
	 * @param none
	 * @return
	 */
	public static long parseLong(String str, long none){	
		long result = none;
		try{
			result = Long.parseLong(str);
		} catch(NumberFormatException e){
			
		}
		
		return result;
	}
	

	/**
	 * checks if a text is boolean and if so it gets the value
	 * @param str
	 * @param none
	 * @return
	 */
	public static boolean parseBoolean(String str, boolean none){	
		boolean result = none;
		try{
			result = Boolean.parseBoolean(str);
		} catch(Exception e){
			
		}
		
		return result;
	}
	
	/**
	 * generates a String from a list of Strings 
	 * @param v the list of strings
	 * @return
	 */
	public static String arrayToString(String[] v){
		String result="[";
		
		for(String s: v){
			result+=s+", ";
		}
		
		result = result.substring(0, result.length()-2);
		result+="]";
		return result;
	}
	
	/**
	 * rounds a double to a 2 decimal number
	 * @param d the double number
	 * @return the 2 decimal double number
	 */
	public static double roundTwoDecimals(double d) {
        	DecimalFormat twoDForm = new DecimalFormat("#.##");
		return Double.valueOf(twoDForm.format(d));
	}
	
	
	
	public static double[] normalizeVector(double[] vector){
		double[] normalizedVector= new double[vector.length];
		
		//compute min max value
		double sum =0;
		for(int i=0;i<vector.length;i++){
			sum+=vector[i];
		}

		for(int i=0;i<vector.length;i++){
			normalizedVector[i] = (double)(vector[i]/sum);

		}
		
		
		return normalizedVector;
	}
	

	/**
	 * log 2 normalization of the vector
	 * @param vector
	 * @return
	 */
	public static double[] normalizeLog2Vector(double[] vector){
		double[] normalizedVector= new double[vector.length];
		
		//compute min max value
		for(int i=0;i<vector.length;i++){
			normalizedVector[i] = (Math.log(vector[i])/Math.log(2));

		}
		
		
		
		return normalizedVector;
	}
	
	/**
	 * computes the sum of a double vector
	 * @param vector
	 * @return
	 */
	public static double computeSum(double[] vector){
		double sum=0;
		for(double v:vector){
			sum+=v;
		}
		
		return sum;
	}
	
	/**
	 * computes the sum of a double vector
	 * @param vector
	 * @return
	 */
	public static double sum(ArrayList<Double> vector){
		double sum=0;
		for(double v:vector){
			sum+=v;
		}
		
		return sum;
	}

	
	
	/**
	 * computes the sum of a int vector
	 * @param vector
	 * @return
	 */
	public static long computeSum(int[] vector){
		long sum=0;
		for(int v:vector){
			sum+=v;
		}
		
		return sum;
	}


	/**
	 * computes the sum of a long vector
	 * @param vector
	 * @return
	 */
	public static long computeSum(long[] vector){
		long sum=0;
		for(long v:vector){
			sum+=v;
		}
		
		return sum;
	}

	
	/**
	 * computes the sum of a double vector
	 * @param vector
	 * @return
	 */
	public static double computeSum(ArrayList<Integer> vector){
		double sum=0;
		for(double v:vector){
			sum+=v;
		}
		
		return sum;
	}
	

	/**
	 * computes the sum of a boolean vector
	 * @param vector
	 * @return
	 */
	public static long computeSum(boolean[] vector){
		long sum=0;
		for(boolean v:vector){
			if(v){
				sum++;
			}
		}
		
		return sum;
	}
	
	/**
	 * gets a maximum value in a vector
	 * @param vector
	 * @return
	 */
	public static double getMax(double[] vector){
		double result = Double.MIN_VALUE;
		
		for(double d: vector){
			if(d>result){
				result=d;
			}
		}
		return result;
	}
	
	/**
	 * gets a maximum value in a vector
	 * @param vector
	 * @return
	 */
	public static int getMax(int[] vector){
		int result = Integer.MIN_VALUE;
		
		for(int d: vector){
			if(d>result){
				result=d;
			}
		}
		return result;
	}
	
	/**
	 * gets a maximum value in a vector
	 * @param vector
	 * @return
	 */
	public static double getMax(ArrayList<Double> vector){
		double result = Double.MIN_VALUE;
		
		for(double d: vector){
			if(d>result){
				result=d;
			}
		}
		return result;
	}

	/**
	 * gets a maximum value in a vector
	 * @param vector
	 * @return
	 */
	public static long getMax(long[] vector){
		long result = Long.MIN_VALUE;
		
		for(long d: vector){
			if(d>result){
				result=d;
			}
		}
		return result;
	}
	
	
	/**
	 * gets a minimum value in a vector
	 * @param vector
	 * @return
	 */
	public static double getMin(double[] vector){
		double result = Double.MAX_VALUE;
		
		for(double d: vector){
			if(d<result){
				result=d;
			}
		}
		return result;
	}
	
	/**
	 * gets a minimum value in a vector
	 * @param vector
	 * @return
	 */
	public static int getMin(int[] vector){
		int result = Integer.MAX_VALUE;
		
		for(int d: vector){
			if(d<result){
				result=d;
			}
		}
		return result;
	}
	

	/**
	 * gets a minimum value in a vector
	 * @param vector
	 * @return
	 */
	public static long getMin(long[] vector){
		long result = Long.MAX_VALUE;
		
		for(long d: vector){
			if(d<result){
				result=d;
			}
		}
		return result;
	}
	
	
	/**
	 * results the procentage a protein is within one interval
	 * @param intStart
	 * @param intEnd
	 * @param position
	 * @param size
	 * @return
	 */
	public static double procentageInInterval(int intStart, int intEnd, int position, int size){
		double result=0;
		
		if(intEnd-intStart<size || (position>intStart && position + size< intEnd)){
			result=1.0;
		}
		
		
		return result;
	}
	
	
	
	/**
	 * computes the average value of a vector
	 * @param vector
	 * @return
	 */
	public static double computeAverage(ArrayList<Double> vector){
		double avg = 0;

		for(double v: vector){
			avg+=v;
		}
		
		if(vector.size()>0){
			avg = avg/vector.size();
		}
		
		return avg;
	}
	
	
	/**
	 * computes the average value of a vector
	 * @param vector
	 * @return
	 */
	public static double computeMean(ArrayList<Integer> vector){
		double avg = 0;

		for(int v: vector){
			avg+=v;
		}
		
		if(vector.size()>0){
			avg = avg/vector.size();
		}
		
		return avg;
	}
	
	/**
	 * computes the variance of a vector
	 * @param vector
	 * @return
	 */
	public static double computeVariance(ArrayList<Double>  vector, double mean){
		double var = 0;

		for(double v: vector){
			var+=(v-mean)*(v-mean);
		}
		
		if(vector.size()>0){
			var = var/vector.size();
		}
		
		return var;
	}
	
	/**
	 * comptes the log base 2 of a number
	 * @param num
	 * @return
	 */
	public static double log2(double num){
		return (Math.log(num)/Math.log(2));
	} 

	/**
	 * comptes the log base 10 of a number
	 * @param num
	 * @return
	 */
	public static double log10(double num){
		return (Math.log(num)/Math.log(10));
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
	 * @param filename
	 * @return
	 * @throws IOException
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
	 * @throws IOException
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
	
	
	/**
	 * checks whether a value is found in an array
	 * @param array
	 * @param value
	 * @return
	 */
	public static boolean containsValue(double[] array, double value){
		for(double val:array){
			if(val==value){
				return true;
			}
		}
		return false;
		
	}
	
	/**
	 * generates a log scale array
	 * @param min
	 * @param max
	 * @return
	 */
	public static ArrayList<Integer> generateLog10ScaleInetegerArray(int min, int max){
		ArrayList<Integer> result=new ArrayList<Integer>();
		
		int minbase = (int) Math.floor(log10(min));
		int maxbase = (int) Math.floor(log10(max));

		for(int i=minbase;i<=maxbase;i++){
			for(int j=1;j<10;j++){
				result.add((int)Math.pow(10, i)*j);
			}
			if(i==maxbase){
				result.add((int)Math.pow(10, i)*10);
			}
		}
		return result;
	} 

	/**
	 * generates a log scale for a double value
	 * @param min
	 * @param max
	 * @param value
	 * @return
	 */
	public static ArrayList<Double> multiplyArray( ArrayList<Integer> scale, double value, int digits){
		ArrayList<Double> result=new ArrayList<Double>();
		for(int i=0; i<scale.size();i++){
			result.add(Math.round(value*scale.get(i)*Math.pow(10, digits))/Math.pow(10, digits));
		}
		
		return result;
	} 
	
}
