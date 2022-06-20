package simulator;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import environment.Cell;

import objects.DNA;
import objects.DNAregion;
import objects.TFspecies;

import utils.CellUtils;
import utils.TFfileParser;
import utils.Utils;

public class GenerateDNAsequence {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		
		Random rand= new Random();		
		String path="biodata";
		int K12DNAlength=4639675;
		
		
		
		ArrayList<Integer> fullCognateCopyNumber=Utils.generateLog10ScaleInetegerArray(1,1000);
		
			
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(1);fullCognateCopyNumber.add(5);fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(50);fullCognateCopyNumber.add(100);fullCognateCopyNumber.add(500);fullCognateCopyNumber.add(1000);
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(50);fullCognateCopyNumber.add(100);

		
		System.out.println(fullCognateCopyNumber);
		//cognateCopyNumber.add(1);		
		//cognateCopyNumber.add(1);cognateCopyNumber.add(10);cognateCopyNumber.add(50);cognateCopyNumber.add(100);cognateCopyNumber.add(1000);
		ArrayList<Integer> nonCognateCopyNumber=new ArrayList<Integer>();nonCognateCopyNumber.add(0);nonCognateCopyNumber.add(10000);nonCognateCopyNumber.add(30000);nonCognateCopyNumber.add(50000);nonCognateCopyNumber.add(70000);
		ArrayList<Double> fullCognateWaitingTime=Utils.multiplyArray(fullCognateCopyNumber, 0.033,3);
		fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(0.066);fullCognateWaitingTime.add(0.33);fullCognateWaitingTime.add(0.66);
		System.out.println(fullCognateWaitingTime);

		
		//cognateWaitingTime.add(0.003302781);cognateWaitingTime.add(0.03302781);cognateWaitingTime.add(0.3302781);cognateWaitingTime.add(3.302781);cognateWaitingTime.add(33.02781);
		int avgCognateCopyNumber = 50;
		double avgWaitingTime = 0.33;

		//the length of the DNA segment
		int DNAlength = 20000;
		
		//base pair composition of the DNA segment
		double DNAproportionOfA=0.246,  DNAproportionOfT=0.246,  DNAproportionOfC=0.254,  DNAproportionOfG=0.254;
		byte[] DNAsequence = CellUtils.generateRandomDNASequence(rand, DNAlength,  DNAproportionOfA,  DNAproportionOfT,  DNAproportionOfC,  DNAproportionOfG);
		DNAsequence = CellUtils.generateDNAStrand(DNAlength,"A");
		DNAregion dnaRegion = new DNAregion("K12",0,DNAlength); 
		
		DNA loadDNA;
		boolean loadedFasta=false;
		File f;
		/*f = new File(path,  "randomDNAsequence.fasta");
		if(f.exists()){
			loadDNA = FastaFileParser.fileParser(path+"/"+ "randomDNAsequence.fasta");
			if(loadDNA!=null && loadDNA.strand!=null && loadDNA.strand.length>0){
				DNAsequence = CellUtils.copySequence(loadDNA.strand);
				loadedFasta = true;
			}
			
		}

		
		if(!loadedFasta){
			CellUtils.printSequence(path, "randomDNAsequence.fasta", DNAsequenceDescription, DNAsequence);
		}*/
		
		String DNAsequenceDescription = "K12:0.."+DNAsequence.length+"; subsequence = K12:0.."+DNAsequence.length+"; copy = 1; boundary = \"reflexive\" ";

		
		//the length of the TF motif
		int TFlength = 20;
		double stopTime=3000;
		int replicates=100;
		int sizeLeft, sizeRight;
		
		TFlength = 1;
		
		//base pair composition of the TF motif
		double TFproportionOfA=71610,  TFproportionOfT=73293,  TFproportionOfC=50158,  TFproportionOfG=50207, total;
		total = TFproportionOfA+TFproportionOfT+TFproportionOfC+TFproportionOfG;
		TFproportionOfA/=total;TFproportionOfC/=total;TFproportionOfG/=total;TFproportionOfT/=total;

		byte[] TF1seq = CellUtils.generateRandomDNASequence(rand, TFlength,  TFproportionOfA,  TFproportionOfT,  TFproportionOfC,  TFproportionOfG);
		byte[] TF2seq = CellUtils.generateRandomDNASequence(rand, TFlength,  TFproportionOfA,  TFproportionOfT,  TFproportionOfC,  TFproportionOfG);

		TF1seq =  CellUtils.generateDNAStrand(TFlength,"C");
		TF2seq =  CellUtils.generateDNAStrand(TFlength,"G");
		
		/*loadedFasta=false;
		f = new File(path,  "randomTF1motifSequence.fasta");
		if(f.exists()){
			loadDNA = FastaFileParser.fileParser(path+"/"+ "randomTF1motifSequence.fasta");
			if(loadDNA!=null && loadDNA.strand!=null && loadDNA.strand.length>0){
				TF1seq = CellUtils.copySequence(loadDNA.strand);
				loadedFasta=true;
			}
		}
		if(!loadedFasta){
			CellUtils.printSequence(path, "randomTF1motifSequence.fasta", "0..20; A="+TFproportionOfA+", T="+TFproportionOfT+", C="+TFproportionOfC+", G="+TFproportionOfG, TF1seq);
		}
		
		
		loadedFasta=false;
		f = new File(path,  "randomTF2motifSequence.fasta");
		if(f.exists()){
			loadDNA = FastaFileParser.fileParser(path+"/"+ "randomTF2motifSequence.fasta");
			if(loadDNA!=null && loadDNA.strand!=null && loadDNA.strand.length>0){
				TF2seq = CellUtils.copySequence(loadDNA.strand);
				loadedFasta=true;
			}
		} 
		if(!loadedFasta){
			CellUtils.printSequence(path, "randomTF2motifSequence.fasta", "0..20; A="+TFproportionOfA+", T="+TFproportionOfT+", C="+TFproportionOfC+", G="+TFproportionOfG, TF2seq);
		}*/
		
		//CellUtils.printSequence(path, "randomTF3motifSequence.fasta", "0..20; A="+TFproportionOfA+", T="+TFproportionOfT+", C="+TFproportionOfC+", G="+TFproportionOfG, TF3seq);
		//CellUtils.printSequence(path, "randomTF4motifSequence.fasta", "0..20; A="+TFproportionOfA+", T="+TFproportionOfT+", C="+TFproportionOfC+", G="+TFproportionOfG, TF4seq);

		
		Cell n = new Cell("params/params.grp", null, true);
		ArrayList<Integer> TSposition1,TSposition2,TSposition3,TSposition4;
		
		int insertPos;
		
		//1a
		byte[] seq_1TF1_1TF2;
		byte[] DNAsequence_1TF1_1TF2;
		int distance;
		String dist="";
		
		
		//1a 2 sites in a uniform distribution landscape and 
		
		 
		sizeLeft= 10; sizeRight=10;
		replicates = 5;
		int seqNo = 0;
		distance = 100;
		
		int motifLength = TF1seq.length +sizeLeft+sizeRight;

		byte[] seq_TF1 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF1seq),CellUtils.generateEmptyDNAStrand(sizeRight));
		byte[] seq_TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF2seq),CellUtils.generateEmptyDNAStrand(sizeRight));
		
		seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(seq_TF1,seq_TF2);
		dist=distance+"s";
		if(distance <0){
			dist=(-distance)+"o";
			seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_TF1,0,seq_TF1.length + distance),seq_TF2);
		} else if(distance > 0){
			seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_TF1, CellUtils.generateEmptyDNAStrand(distance)),seq_TF2);
		} 

		insertPos = (int) Math.floor((DNAsequence.length-motifLength+1)/2 - seq_1TF1_1TF2.length/2);

		
		DNAsequence_1TF1_1TF2 = CellUtils.replaceDNAseq(DNAsequence,seq_1TF1_1TF2,insertPos);
		
		TSposition1=new ArrayList<Integer>();TSposition2=new ArrayList<Integer>();TSposition3=new ArrayList<Integer>();TSposition4=new ArrayList<Integer>();
		TSposition1.add(insertPos);TSposition2.add(insertPos+motifLength+distance);
		
		CellUtils.printSequence(path, "randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", DNAsequenceDescription, DNAsequence_1TF1_1TF2);
		generateFiles(n, rand, "biodata/TF_1lacI_1nc.csv", dnaRegion, K12DNAlength,"randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", fullCognateCopyNumber,nonCognateCopyNumber, fullCognateWaitingTime, avgWaitingTime, avgCognateCopyNumber, TSposition1,TSposition2, TF1seq , TF2seq, sizeLeft, sizeRight, "biodata", "TFfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "TSfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "params/params_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".grp", "run_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".sh", "rand_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"", stopTime, replicates, loadedFasta);		
		
		
		
		
		/*ArrayList<Integer> distances =new ArrayList<Integer>(); distances.add(-5); distances.add(0); distances.add(5); distances.add(50); distances.add(100);
		distances =new ArrayList<Integer>(); distances.add(-2);
		avgCognateCopyNumber = 10;
		avgWaitingTime = 0.33;
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(1);fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(100);
		fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(0.033);fullCognateWaitingTime.add(0.33);fullCognateWaitingTime.add(3.3);

		for(int i: distances){
			
			sizeLeft= 10; sizeRight=10;
			replicates = 5;
			int seqNo = 0;
			distance = i;
			
			int motifLength = TF1seq.length +sizeLeft+sizeRight;

			byte[] seq_TF1 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF1seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			byte[] seq_TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF2seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			
			seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(seq_TF1,seq_TF2);
			dist=distance+"s";
			if(distance <0){
				dist=(-distance)+"o";
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_TF1,0,seq_TF1.length + distance),seq_TF2);
			} else if(distance > 0){
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_TF1, CellUtils.generateEmptyDNAStrand(distance)),seq_TF2);
			} 

			insertPos = (int) Math.floor((DNAsequence.length-motifLength+1)/2 - seq_1TF1_1TF2.length/2);

			
			DNAsequence_1TF1_1TF2 = CellUtils.replaceDNAseq(DNAsequence,seq_1TF1_1TF2,insertPos);
			
			TSposition1=new ArrayList<Integer>();TSposition2=new ArrayList<Integer>();TSposition3=new ArrayList<Integer>();TSposition4=new ArrayList<Integer>();
			TSposition1.add(insertPos);TSposition2.add(insertPos+motifLength+distance);
			
			CellUtils.printSequence(path, "randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", DNAsequenceDescription, DNAsequence_1TF1_1TF2);
			generateFiles(n, rand, "biodata/TF_1lacI_1nc.csv", dnaRegion, K12DNAlength,"randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", fullCognateCopyNumber,nonCognateCopyNumber, fullCognateWaitingTime, avgWaitingTime, avgCognateCopyNumber, TSposition1,TSposition2, TF1seq , TF2seq, sizeLeft, sizeRight, "biodata", "TFfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "TSfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "params/params_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".grp", "run_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".sh", "rand_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"", stopTime, replicates);		
			
		}*/
		
		
		
		/*f = new File(path,  "randomDNAsequence.fasta");
		if(f.exists()){
			loadDNA = FastaFileParser.fileParser(path+"/"+ "randomDNAsequence.fasta");
			if(loadDNA!=null && loadDNA.strand!=null && loadDNA.strand.length>0){
				DNAsequence = CellUtils.copySequence(loadDNA.strand);
				loadedFasta = true;
			}
			
		}
		DNAsequenceDescription = "K12:0.."+DNAsequence.length+"; subsequence = K12:0.."+DNAsequence.length+"; copy = 1; boundary = \"reflexive\" ";

		
		if(!loadedFasta){
			CellUtils.printSequence(path, "randomDNAsequence.fasta", DNAsequenceDescription, DNAsequence);
		}
		loadedFasta=false;
		f = new File(path,  "randomTF1motifSequence.fasta");
		if(f.exists()){
			loadDNA = FastaFileParser.fileParser(path+"/"+ "randomTF1motifSequence.fasta");
			if(loadDNA!=null && loadDNA.strand!=null && loadDNA.strand.length>0){
				TF1seq = CellUtils.copySequence(loadDNA.strand);
				loadedFasta=true;
			}
		}
		if(!loadedFasta){
			CellUtils.printSequence(path, "randomTF1motifSequence.fasta", "0..20; A="+TFproportionOfA+", T="+TFproportionOfT+", C="+TFproportionOfC+", G="+TFproportionOfG, TF1seq);
		}
		
		
		loadedFasta=false;
		f = new File(path,  "randomTF2motifSequence.fasta");
		if(f.exists()){
			loadDNA = FastaFileParser.fileParser(path+"/"+ "randomTF2motifSequence.fasta");
			if(loadDNA!=null && loadDNA.strand!=null && loadDNA.strand.length>0){
				TF2seq = CellUtils.copySequence(loadDNA.strand);
				loadedFasta=true;
			}
		} 
		if(!loadedFasta){
			CellUtils.printSequence(path, "randomTF2motifSequence.fasta", "0..20; A="+TFproportionOfA+", T="+TFproportionOfT+", C="+TFproportionOfC+", G="+TFproportionOfG, TF2seq);
		}
		
		ArrayList<Integer> distances =new ArrayList<Integer>(); distances.add(-2); distances.add(0); distances.add(5); distances.add(50); distances.add(100);
		avgCognateCopyNumber = 10;
		avgWaitingTime = 0.33;
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(1);fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(100);
		fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(0.033);fullCognateWaitingTime.add(0.33);fullCognateWaitingTime.add(3.3);

		for(int i: distances){
			
			sizeLeft= 1; sizeRight=1;
			replicates = 5;
			int seqNo = 1;
			distance = i;
			
			int motifLength = TF1seq.length +sizeLeft+sizeRight;

			byte[] seq_TF1 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF1seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			byte[] seq_TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF2seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			
			seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(seq_TF1,seq_TF2);
			dist=distance+"s";
			if(distance <0){
				dist=(-distance)+"o";
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_TF1,0,seq_TF1.length + distance),seq_TF2);
			} else if(distance > 0){
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_TF1, CellUtils.generateEmptyDNAStrand(distance)),seq_TF2);
			} 

			insertPos = (int) Math.floor((DNAsequence.length-motifLength+1)/2 - seq_1TF1_1TF2.length/2);

			
			DNAsequence_1TF1_1TF2 = CellUtils.replaceDNAseq(DNAsequence,seq_1TF1_1TF2,insertPos);
			
			TSposition1=new ArrayList<Integer>();TSposition2=new ArrayList<Integer>();TSposition3=new ArrayList<Integer>();TSposition4=new ArrayList<Integer>();
			TSposition1.add(insertPos);TSposition2.add(insertPos+motifLength+distance);
			
			CellUtils.printSequence(path, "randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", DNAsequenceDescription, DNAsequence_1TF1_1TF2);
			generateFiles(n, rand, "biodata/TF_1lacI_1nc.csv", dnaRegion, K12DNAlength,"randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", fullCognateCopyNumber,nonCognateCopyNumber, fullCognateWaitingTime, avgWaitingTime, avgCognateCopyNumber, TSposition1,TSposition2, TF1seq , TF2seq, sizeLeft, sizeRight, "biodata", "TFfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "TSfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "params/params_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".grp", "run_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".sh", "rand_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"", stopTime, replicates);		
		}*/
	
		
		/*ArrayList<Integer> distances =new ArrayList<Integer>(); distances.add(-5); distances.add(0); distances.add(5); distances.add(50); distances.add(100);
		avgCognateCopyNumber = 10;
		avgWaitingTime = 0.33;
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(1);fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(100);
		fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(0.033);fullCognateWaitingTime.add(0.33);fullCognateWaitingTime.add(3.3);
		stopTime = 3000*282;
		for(int i: distances){
			
			sizeLeft= 10; sizeRight=10;
			replicates = 5;
			int seqNo = 0;
			distance = i;
			
			int motifLength = TF1seq.length +sizeLeft+sizeRight;

			byte[] seq_TF1 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF1seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			byte[] seq_TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF2seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			
			seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(seq_TF1,seq_TF2);
			dist=distance+"s";
			if(distance <0){
				dist=(-distance)+"o";
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_TF1,0,seq_TF1.length + distance),seq_TF2);
			} else if(distance > 0){
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_TF1, CellUtils.generateEmptyDNAStrand(distance)),seq_TF2);
			} 

			insertPos = (int) Math.floor((DNAsequence.length-motifLength+1)/2 - seq_1TF1_1TF2.length/2);

			
			DNAsequence_1TF1_1TF2 = CellUtils.replaceDNAseq(DNAsequence,seq_1TF1_1TF2,insertPos);
			
			TSposition1=new ArrayList<Integer>();TSposition2=new ArrayList<Integer>();TSposition3=new ArrayList<Integer>();TSposition4=new ArrayList<Integer>();
			TSposition1.add(insertPos);TSposition2.add(insertPos+motifLength+distance);
			
			CellUtils.printSequence(path, "randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", DNAsequenceDescription, DNAsequence_1TF1_1TF2);
			generateFiles(n, rand, "biodata/TF_1lacI_1nc_3d.csv", dnaRegion, K12DNAlength,"randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", fullCognateCopyNumber,nonCognateCopyNumber, fullCognateWaitingTime, avgWaitingTime, avgCognateCopyNumber, TSposition1,TSposition2, TF1seq , TF2seq, sizeLeft, sizeRight, "biodata", "TFfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "TSfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "params/params_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".grp", "run_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".sh", "rand_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"", stopTime, replicates);		
			
		}*/
		
		
		//3D only strong
		/*ArrayList<Integer> distances =new ArrayList<Integer>(); distances.add(-5); distances.add(0); distances.add(5); distances.add(50); distances.add(100);
		avgCognateCopyNumber = 10;
		avgWaitingTime = 14.85;
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(1);fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(100);
		fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(1.485);fullCognateWaitingTime.add(14.85);fullCognateWaitingTime.add(148.5);
		stopTime = 18.6;
		for(int i: distances){
			
			sizeLeft= 10; sizeRight=10;
			replicates = 5;
			int seqNo = 0;
			distance = i;
			
			int motifLength = TF1seq.length +sizeLeft+sizeRight;

			byte[] seq_TF1 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF1seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			byte[] seq_TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF2seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			
			seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(seq_TF1,seq_TF2);
			dist=distance+"s";
			if(distance <0){
				dist=(-distance)+"o";
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_TF1,0,seq_TF1.length + distance),seq_TF2);
			} else if(distance > 0){
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_TF1, CellUtils.generateEmptyDNAStrand(distance)),seq_TF2);
			} 

			insertPos = (int) Math.floor((DNAsequence.length-motifLength+1)/2 - seq_1TF1_1TF2.length/2);

			
			DNAsequence_1TF1_1TF2 = CellUtils.replaceDNAseq(DNAsequence,seq_1TF1_1TF2,insertPos);
			
			TSposition1=new ArrayList<Integer>();TSposition2=new ArrayList<Integer>();TSposition3=new ArrayList<Integer>();TSposition4=new ArrayList<Integer>();
			TSposition1.add(insertPos);TSposition2.add(insertPos+motifLength+distance);
			
			CellUtils.printSequence(path, "randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+"_strong.fasta", DNAsequenceDescription, DNAsequence_1TF1_1TF2);
			generateFiles(n, rand, "biodata/TF_1lacI_1nc_3d.csv", dnaRegion, K12DNAlength,"randomDNAsequence_1TF1_1TF2_"+dist+"_seq"+seqNo+"_strong.fasta", fullCognateCopyNumber,nonCognateCopyNumber, fullCognateWaitingTime, avgWaitingTime, avgCognateCopyNumber, TSposition1,TSposition2, TF1seq , TF2seq, sizeLeft, sizeRight, "biodata", "TFfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"_strong.csv", "TSfile_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"_strong.csv", "params/params_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"_strong.grp", "run_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"_strong.sh", "rand_1TF1_1TF2_0nc_"+dist+"_seq"+seqNo+"_strong", stopTime, replicates);		
			
		}*/
		
		//site duplicates
		/*ArrayList<Integer> distances =new ArrayList<Integer>(); distances.add(-5); distances.add(0); distances.add(5); distances.add(50); distances.add(100);
		
		distances =new ArrayList<Integer>(); distances.add(10); distances.add(20); distances.add(40); distances.add(75);
		avgCognateCopyNumber = 10;
		//fullCognateCopyNumber.add(1);fullCognateCopyNumber.add(10);avgWaitingTime = 0.33;
		//fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(1);fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(100);
		//fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(0.033);fullCognateWaitingTime.add(0.33);fullCognateWaitingTime.add(3.3);
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(10);
		fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(0.033);
		avgWaitingTime = 0.33;
		//fullCognateWaitingTime.add(0.033);fullCognateWaitingTime.add(3.3);
		
		for(int i: distances){
			
			sizeLeft= 10; sizeRight=10;
			replicates = 400;
			int seqNo = 0;
			distance = i;
			
			int motifLength = TF2seq.length +sizeLeft+sizeRight;

			byte[] seq_TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF2seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			
			byte[] seq_1TF2_1TF2 = CellUtils.concatenateDNAseq(seq_TF2,seq_TF2);
			dist=distance+"s";
			if(distance <0){
				dist=(-distance)+"o";
				seq_1TF2_1TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_TF2,0,seq_TF2.length + distance),seq_TF2);
			} else if(distance > 0){
				seq_1TF2_1TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_TF2, CellUtils.generateEmptyDNAStrand(distance)),seq_TF2);
			} 

			insertPos = (int) Math.floor((DNAsequence.length-motifLength+1)/2 - seq_1TF2_1TF2.length/2);

			
			DNAsequence_1TF1_1TF2 = CellUtils.replaceDNAseq(DNAsequence,seq_1TF2_1TF2,insertPos);
			
			TSposition1=new ArrayList<Integer>();TSposition2=new ArrayList<Integer>();TSposition3=new ArrayList<Integer>();TSposition4=new ArrayList<Integer>();
			TSposition2.add(insertPos);TSposition2.add(insertPos+motifLength+distance);
			
			CellUtils.printSequence(path, "randomDNAsequence_2TF2_"+dist+"_seq"+seqNo+".fasta", DNAsequenceDescription, DNAsequence_1TF1_1TF2);
			generateFiles(n, rand, "biodata/TF_1lacI_1nc.csv", dnaRegion, K12DNAlength,"randomDNAsequence_2TF2_"+dist+"_seq"+seqNo+".fasta", fullCognateCopyNumber,nonCognateCopyNumber, fullCognateWaitingTime, avgWaitingTime, avgCognateCopyNumber, TSposition1,TSposition2, TF1seq , TF2seq, sizeLeft, sizeRight, "biodata", "TFfile_2TF2_0nc_"+dist+"_seq"+seqNo+".csv", "TSfile_2TF2_0nc_"+dist+"_seq"+seqNo+".csv", "params/params_2TF2_0nc_"+dist+"_seq"+seqNo+".grp", "run_2TF2_0nc_"+dist+"_seq"+seqNo+".sh", "rand_2TF2_0nc_"+dist+"_seq"+seqNo+"", stopTime, replicates,false);		
			
		}*/
		
		
		/*ArrayList<Integer> distances =new ArrayList<Integer>(); distances.add(-5); distances.add(0); distances.add(5); distances.add(50); distances.add(100);
		avgCognateCopyNumber = 10;
		avgWaitingTime = 0.33;
		fullCognateCopyNumber = new ArrayList<Integer>();fullCognateCopyNumber.add(5);fullCognateCopyNumber.add(10);fullCognateCopyNumber.add(100);
		fullCognateWaitingTime = new ArrayList<Double>();fullCognateWaitingTime.add(0.033);fullCognateWaitingTime.add(0.33);fullCognateWaitingTime.add(3.3);

		for(int i: distances){
			
			sizeLeft= 10; sizeRight=10;
			replicates = 5;
			int seqNo = 0;
			distance = i;
			
			int motifLength = TF1seq.length +sizeLeft+sizeRight;

			byte[] seq_TF1 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF1seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			byte[] seq_TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(CellUtils.generateEmptyDNAStrand(sizeLeft), TF2seq),CellUtils.generateEmptyDNAStrand(sizeRight));
			
			seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(seq_TF2,seq_TF1);
			dist=distance+"s";
			if(distance <0){
				dist=(-distance)+"o";
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_TF1,0,seq_TF1.length + distance),seq_TF2);
			} else if(distance > 0){
				seq_1TF1_1TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_TF1, CellUtils.generateEmptyDNAStrand(distance)),seq_TF2);
			} 

			byte[] seq_1TF1_2TF2 = CellUtils.concatenateDNAseq(seq_1TF1_1TF2,seq_TF1);
			dist=distance+"s";
			if(distance <0){
				dist=(-distance)+"o";
				seq_1TF1_2TF2 = CellUtils.concatenateDNAseq(CellUtils.extractSubSequence(seq_1TF1_1TF2,0,seq_1TF1_1TF2.length + distance),seq_TF1);
			} else if(distance > 0){
				seq_1TF1_2TF2 = CellUtils.concatenateDNAseq(CellUtils.concatenateDNAseq(seq_1TF1_1TF2, CellUtils.generateEmptyDNAStrand(distance)),seq_TF1);
			} 
			
			insertPos = (int) Math.floor((DNAsequence.length-motifLength+1)/2 - seq_1TF1_2TF2.length/2);

			
			byte[] DNAsequence_1TF1_2TF2 = CellUtils.replaceDNAseq(DNAsequence,seq_1TF1_2TF2,insertPos);
			
			TSposition1=new ArrayList<Integer>();TSposition2=new ArrayList<Integer>();TSposition3=new ArrayList<Integer>();TSposition4=new ArrayList<Integer>();
			TSposition2.add(insertPos+motifLength+distance);TSposition1.add(insertPos);TSposition1.add(insertPos+2*(motifLength+distance));
			
			CellUtils.printSequence(path, "randomDNAsequence_2TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", DNAsequenceDescription, DNAsequence_1TF1_2TF2);
			generateFiles(n, rand, "biodata/TF_1lacI_1nc.csv", dnaRegion, K12DNAlength,"randomDNAsequence_2TF1_1TF2_"+dist+"_seq"+seqNo+".fasta", fullCognateCopyNumber,nonCognateCopyNumber, fullCognateWaitingTime, avgWaitingTime, avgCognateCopyNumber, TSposition1,TSposition2, TF1seq , TF2seq, sizeLeft, sizeRight, "biodata", "TFfile_2TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "TSfile_2TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".csv", "params/params_2TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".grp", "run_2TF1_1TF2_0nc_"+dist+"_seq"+seqNo+".sh", "rand_1TF1_2TF2_0nc_"+dist+"_seq"+seqNo+"", stopTime, replicates,true);		
			
		}*/
		
	}
	
	
	public static void generateFiles(Cell n, Random randomGenerator, String TFfile, DNAregion dnaRegion, int maxDNAsize, String DNAsequenceFile, ArrayList<Integer> cognateCopyNumber, ArrayList<Integer> nonCognateCopyNumber,ArrayList<Double> cognateWaitingTime,	double avgCognateWaitingTime,	int avgCognateCopyNumber, ArrayList<Integer> TSposition1, ArrayList<Integer> TSposition2, byte[] TF1seq , byte[] TF2seq, int sizeLeft, int sizeRight, String path, String outputTFfile, String bufferTSfile, String paramsFile, String scriptFile, String resultFilename, double stopTime, int replicates, boolean full) throws IOException{
		
		BufferedWriter out ;
		n.ip.TS_FILE.value=path+"/"+bufferTSfile;
		n.ip.STOP_TIME.value=stopTime;
		
		
		//TSfile
		String bufferTSstr, bufferTSgroupStr;
		out = new BufferedWriter(new FileWriter(new File(path,bufferTSfile)));
		
		String allBuffer="";
		
		bufferTSgroupStr="";
		if(TSposition1!=null && !TSposition1.isEmpty()){
			for(int i=0;i<TSposition1.size();i++){
				bufferTSstr="TF1:K12:"+TSposition1.get(i)+".."+(TSposition1.get(i)+1)+":-1";	
				out.write(bufferTSstr);
				out.newLine();
				bufferTSgroupStr += i>0?" AND ": "";
				bufferTSgroupStr+=bufferTSstr;
				if(i>0){
					out.write(bufferTSgroupStr);
					out.newLine();
					out.write(bufferTSgroupStr.replaceAll("AND", "OR"));
					out.newLine();
				}
			}
		}
		allBuffer+=bufferTSgroupStr;
		
		
		bufferTSgroupStr="";
		if(TSposition2!=null && !TSposition2.isEmpty()){
			for(int i=0;i<TSposition2.size();i++){
				bufferTSstr="TF2:K12:"+TSposition2.get(i)+".."+(TSposition2.get(i)+1)+":-1";	
				out.write(bufferTSstr);
				out.newLine();
				bufferTSgroupStr += i>0?" AND ": "";
				bufferTSgroupStr +=bufferTSstr;
				if(i>0){
					out.write(bufferTSgroupStr);
					out.newLine();
					out.write(bufferTSgroupStr.replaceAll("AND", "OR"));
					out.newLine();
				}
			}
		}
			
		if(!allBuffer.isEmpty() && !bufferTSgroupStr.isEmpty()){
			out.write(allBuffer+" AND "+ bufferTSgroupStr);
			out.newLine();
			out.write(allBuffer.replaceAll("AND", "OR")+" OR "+bufferTSgroupStr.replaceAll("AND", "OR"));
			out.newLine();
		}
		
		out.flush();
		out.close();
		
		int nonCognateID=0, cognateID=0;
		
		//TF file
		TFfileParser TFparser=new TFfileParser(n, TFfile, Utils.generateNextInteger(randomGenerator, n.ip.TF_COPY_NUMBER_MIN.value, n.ip.TF_COPY_NUMBER_MAX.value), n.ip.TF_ES.value, n.ip.TF_SIZE_LEFT.value, n.ip.TF_SIZE_RIGHT.value, n.ip.TF_ASSOC_RATE.value, n.dna.strand.length, n.ip.TF_UNBINDING_PROBABILITY.value, n.ip.TF_SLIDE_LEFT_PROBABILITY.value, n.ip.TF_SLIDE_RIGHT_PROBABILITY.value, n.ip.TF_JUMPING_PROBABILITY.value,
				n.ip.TF_HOP_STD_DISPLACEMENT.value, n.ip.TF_SPECIFIC_WAITING_TIME.value, n.ip.TF_STEP_LEFT_SIZE.value, n.ip.TF_STEP_RIGHT_SIZE.value, n.ip.TF_UNCORRELATED_DISPLACEMENT_SIZE.value,
				n.ip.TF_STALLS_IF_BLOCKED.value, n.ip.TF_COLLISION_UNBIND_PROBABILITY.value, n.ip.TF_AFFINITY_LANDSCAPE_ROUGHNESS.value, n.ip.TF_PREBOUND_PROPORTION.value, n.ip.TF_PREBOUND_TO_HIGHEST_AFFINITY.value, n.ip.TF_IS_IMMOBILE.value,dnaRegion, n.ip.IS_BIASED_RANDOM_WALK.value,n.ip.IS_TWO_STATE_RANDOM_WALK.value,
				n.ip.TF_REPRESSION_RATE.value, n.ip.TF_DEREPRESSION_RATE.value, n.ip.TF_REPR_LEN_LEFT.value, n.ip.TF_REPR_LEN_RIGHT.value);
		if(TFparser.parsed){
			// load TF species
			if(TFparser.data==null || TFparser.data.isEmpty()){
				System.out.println("Could not load the default TF file");
				System.exit(0);
			}
				
			n.TFspecies = new TFspecies[TFparser.data.size()];
			for(int i=0;i< TFparser.data.size();i++){
				n.TFspecies[i] =  TFparser.data.get(i);
				n.moleculesCopyNumber+=n.TFspecies[i].copyNumber;
				if(n.TFspecies[i].name.equals("lacI")){
					n.TFspecies[i].name="TF1";
					cognateID = i;
					n.TFspecies[cognateID].sizeLeft = sizeLeft;
					n.TFspecies[cognateID].sizeRight = sizeRight;
					System.out.println(n.TFspecies[i]);
					
				} else{
					nonCognateID = i;
				}
			}
			
			double lambda = (double)dnaRegion.size()/(double)maxDNAsize;
			double f = 0.9;
			n.TFspecies[cognateID].assocRate *=(lambda-f*lambda)/(1-f*lambda); 
			n.TFspecies[cognateID].name="TF1";

			n.TFspecies[nonCognateID].copyNumber = 0;
			n.TFspecies[nonCognateID].copyNumber *= (int)Math.round(dnaRegion.size()/maxDNAsize);

			String bufferTFfile = outputTFfile, bufferParamsFile, bufferScriptFile;
			
			
			//generateTFfiles
			for(double specificWaitingTime: cognateWaitingTime){

				for(int copyNumber: cognateCopyNumber){
				
					if(full || (copyNumber == avgCognateCopyNumber || specificWaitingTime == avgCognateWaitingTime)){
					
						bufferTFfile = outputTFfile.replaceAll(".csv", "_"+copyNumber+"TF2cn_"+specificWaitingTime+"TF2swt.csv");
						
						out = new BufferedWriter(new FileWriter(new File(path,bufferTFfile)));
			
						out.write(n.TFspecies[0].headerToString(true));
			    			out.newLine();
			
						if(n.TFspecies[nonCognateID].copyNumber>0){
							out.write(n.TFspecies[nonCognateID].toString(n, true));
			    				out.newLine();
						}
						
						if(!TSposition1.isEmpty()){
							n.TFspecies[cognateID].name="TF1";
							n.TFspecies[cognateID].copyNumber = avgCognateCopyNumber;
							n.TFspecies[cognateID].specificWaitingTime = avgCognateWaitingTime;
							n.TFspecies[cognateID].dbd = CellUtils.copySequence(TF1seq);
							out.write(n.TFspecies[cognateID].toString(n, true));
							out.newLine();
						}						
		
						
						if(!TSposition2.isEmpty()){
							n.TFspecies[cognateID].copyNumber = copyNumber;
							n.TFspecies[cognateID].specificWaitingTime = specificWaitingTime;
							n.TFspecies[cognateID].name="TF2";
							n.TFspecies[cognateID].dbd = CellUtils.copySequence(TF2seq);
							out.write(n.TFspecies[cognateID].toString(n, true));
							out.newLine();
						}
						
						
						out.flush();
						out.close();
						
						
						bufferParamsFile = paramsFile.replaceAll(".grp", "_"+copyNumber+"TF2cn_"+specificWaitingTime+"TF2swt.grp");
						n.ip.TF_FILE.value=path+"/"+bufferTFfile;
						n.ip.OUTPUT_FILENAME.value=resultFilename+ "_"+copyNumber+"TF2cn_"+specificWaitingTime+"TF2swt";
						n.ip.OUTPUT_FOLDER.value="results";
	
						n.ip.DNA_SEQUENCE_FILE.value=path+"/"+DNAsequenceFile;
						n.ip.exportParameterFile(bufferParamsFile);
						
						
						bufferScriptFile = scriptFile.replaceAll(".sh", "_"+copyNumber+"TF2cn_"+specificWaitingTime+"TF2swt.sh");
						generateScriptFile(bufferScriptFile,bufferParamsFile,replicates);
					}
				}
			}
			
			
			
			
			
    		
		}	
		
	}
	
	
	/**
	 * generates a bash script file to run the simulations for a certain parameter file and for a number of replicates times
	 * @param filename
	 * @param paramsFile
	 * @param replicates
	 * @throws IOException
	 */
	public static void generateScriptFile(String filename, String paramsFile, int replicates) throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(filename)));

		out.write("#!/bin/bash");
		out.newLine();
		out.write("i=\"0\"");
		out.newLine();
		out.write("j=\""+replicates+"\"");
		out.newLine();
		out.write("while [ $i -lt $j ]");
		out.newLine();
		out.write("do");
		out.newLine();		
		out.write("   echo \"sample $i\"");
		out.newLine();	
		out.write("   java -Xms32m -Xmx2048m -jar gripCLI.jar "+paramsFile+" 100 -1 false");
		out.newLine();	
		out.write("   i=$[$i+1]");
		out.newLine();	
		out.write("done");
		out.newLine();	
			
		out.flush();
		out.close();

	}
	

}