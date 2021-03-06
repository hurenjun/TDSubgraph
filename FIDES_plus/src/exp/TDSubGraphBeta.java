package exp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.TreeSet;

import algorithm.BoundedProbing;
import algorithm.ConvertAG;
import algorithm.IntervalSelector;
import algorithm.MergeCCs;
import algorithm.SpanningTree;
import algorithm.StrongPrunning;
import graph.DenseSubGraph;
import graph.Edge;
import graph.Graph;
import graph.MGraph;
import graph.TGraph;

public class TDSubGraphBeta {  // Time-envolving Dense SubGraph 
	
	public static void main(String[] args) throws IOException, InterruptedException {
		String option = args[0], srcGraph = args[1];
		// data-driven Find Dense Subgraph
		if (option.equals("-run")) {	
			int startT = Integer.valueOf(args[2]), endT = Integer.valueOf(args[3]);
			int k = Integer.valueOf(args[4]);
			run(option, srcGraph, startT, endT, k);
		}
		// data-driven Find Dense Subgraph given intervals (interval.txt)
		else if (option.equals("-runIntv")) {
			int startT = Integer.valueOf(args[2]), endT = Integer.valueOf(args[3]);
			int k = Integer.valueOf(args[4]);
			runIntv(option, srcGraph, startT, endT, k);
		}
		// compute dense subgraph on aggregate graphs, given time intervals (interval.txt)
		else if (option.equals("-computeADS")) {
			int startT = Integer.valueOf(args[2]), endT = Integer.valueOf(args[3]);
			computeADS(option, srcGraph, startT, endT);
		}
		// computeADS with different combinations of optimizations
		else if (option.equals("-ADSOPT")) {
			int startT = Integer.valueOf(args[2]), endT = Integer.valueOf(args[3]);
			String intv = args[4];
			ADSOPT(option, srcGraph, startT, endT, intv);
		}
		// compute evolving convergence rate
		else if (option.equals("-ecrate")) {
			int startT = Integer.valueOf(args[2]), endT = Integer.valueOf(args[3]);
			ecrate(option, srcGraph, startT, endT);
		}
	}
	
	static void run(String option, String srcGraph, int startT, int endT, int k) {
		System.out.println("FIDES_plus " + option + " " + srcGraph + " " + 
				startT + " " + endT + " " + k);
		// Step 0: loading temporal graph
		long st = System.currentTimeMillis();
		TGraph tg = new TGraph(srcGraph, startT, endT);
		System.out.println("Loading took:\t" + (System.currentTimeMillis()-st));
		st = System.currentTimeMillis();
		
		// Step 1: selecting candidate time interval
		if (k % 2 == 1) k++;
		Edge[] topK = IntervalSelector.doIntervalSelect(tg, k);
		System.out.println("Interval took:\t" + (System.currentTimeMillis()-st));
		
		// Step 2: ComputeADS
		//st = System.currentTimeMillis();
		double bestResult = -1;
		Edge bestIntv = new Edge(-1, -1);
		TreeSet<Edge> bestEdgeSet = null;
		//topK[k - 1] = new Edge(1040, 1048);
		for (int i = 0; i < k; i++) {
			long time = System.currentTimeMillis();
			//System.out.print(i + "[" + topK[i].getS() + "," + topK[i].getT() + "]\t");
			Edge tunedIntv = new Edge();
			Graph ag = tg.aggregateGraph(topK[i].getS(), topK[i].getT());
			if (ag == null) continue;	// time interval error
			//System.out.print(ag.getPosCDen() + "\t");
			MGraph convAG = ConvertAG.doConvertAg(ag);
			//System.out.print("ConvertAg:" + (System.currentTimeMillis()-time) + "\t");
			time = System.currentTimeMillis();
			
			// using all 3 optimization techniques
			MGraph mergCAG = MergeCCs.doMergeCCsDTrees(convAG);
			//System.out.print("MergeCCs:" + (System.currentTimeMillis()-time) + "\t");
			time = System.currentTimeMillis();
			
			MGraph mstAG = SpanningTree.doMinimumST(mergCAG);
			DenseSubGraph dsg = StrongPrunning.doStrongPrunning(mstAG);
			//System.out.print("MST&SP:" + (System.currentTimeMillis()-time) + "\t");
			time = System.currentTimeMillis();
			
			BoundedProbing.bp(mergCAG, dsg);
			TreeSet<Edge> dsgEdges = new TreeSet<Edge>();
			double finalDSGWeight = SpanningTree.doMinimumST(convAG, BoundedProbing.getUsed(), 
					MergeCCs.getFinalNodeID(), dsgEdges);
			finalDSGWeight = IntervalSelector.tuning(tg, dsgEdges, topK[i], finalDSGWeight, 
					tunedIntv);
			//System.out.print("BP&Rest:" + (System.currentTimeMillis()-time) + "\t");
			//System.out.print("[" + (tunedIntv.getS()) + "," + (tunedIntv.getT()) + "]" 
			//		+ finalDSGWeight + "\t");
			if (finalDSGWeight > bestResult) {
				bestResult = finalDSGWeight;
				bestIntv.copy(tunedIntv);
				bestEdgeSet = new TreeSet<Edge>(dsgEdges);
			}
			
			// (2) using only Strong Pruning + Bounded Probing		
			time = System.currentTimeMillis();
			MGraph mstAG2 = SpanningTree.doMinimumST(convAG);
			DenseSubGraph dsg2 = StrongPrunning.doStrongPrunning(mstAG2);
			//System.out.print("MST&SP:" + (System.currentTimeMillis()-time) + "\t");
			time = System.currentTimeMillis();
			
			BoundedProbing.bp(convAG, dsg2);
			TreeSet<Edge> dsgEdges2 = new TreeSet<Edge>();
			double finalDSGWeight2 = SpanningTree.doMinimumST(convAG, 
					BoundedProbing.getUsed(), dsgEdges2);
			finalDSGWeight2 = IntervalSelector.tuning(tg, dsgEdges2, topK[i], finalDSGWeight2, 
					tunedIntv);
			//System.out.println("BP&Rest:" + (System.currentTimeMillis()-time) + "\t");
			time = System.currentTimeMillis();
			
			//System.out.println("[" + (tunedIntv.getS()) + "," + (tunedIntv.getT()) + "]" 
			//		+ finalDSGWeight2 + "\t");
			if (finalDSGWeight2 > bestResult) {
				bestResult = finalDSGWeight2;
				bestIntv.copy(tunedIntv);
				bestEdgeSet = new TreeSet<Edge>(dsgEdges2);
			}
		}
		System.out.println("[" + (bestIntv.getS() + startT) + "," + 
				(bestIntv.getT() + startT) + "]\t" + bestResult);
		//Graph ag = tg.aggregateGraph(bestIntv.getS(), bestIntv.getT());
		//System.out.println(ag.computeSumOfEdgeW(bestEdgeSet) + "checking");
		System.out.println("Evaluation took:\t" + (System.currentTimeMillis()-st));
		System.out.println();
		// note the id2name mapping of nodes in TGraph tg, if printing the dsgEdges
	}
	
	static void runIntv(String option, String srcGraph, int startT, int endT, int k) 
			throws NumberFormatException, IOException {
		System.out.println("FIDES_plus " + option + " " + srcGraph + " " + 
				startT + " " + endT + " " + k);
		// Step 0: loading temporal graph
		long st = System.currentTimeMillis();
		TGraph tgAnces = new TGraph(srcGraph, startT, endT);
		System.out.println("Loading took:\t" + (System.currentTimeMillis()-st));
		int[] intvL = new int[1000];
		int[] intvR = new int[1000];
		int cnt = 0;
		BufferedReader br = new BufferedReader(new FileReader("interval.txt"));
		String line = null;
		char sepChar = ',';
		while((line = br.readLine()) != null) {
			if (line.length() == 0) continue;
			int sep = line.indexOf(sepChar);
			intvL[cnt] = Integer.valueOf(line.substring(0, sep));
			intvR[cnt] = Integer.valueOf(line.substring(sep + 1));
			cnt++;
		}
		br.close();
		
		for (int intv = 0; intv < cnt; intv++) {
			System.out.println("TDSubGraphBeta " + option + " " + srcGraph + " " + 
					intvL[intv] + " " + intvR[intv] + " " + k);
			TGraph tg = new TGraph(tgAnces, intvL[intv], intvR[intv]);
			
			st = System.currentTimeMillis();
			// Step 1: selecting candidate time interval
			if (k % 2 == 1) k++;
			Edge[] topK = IntervalSelector.doIntervalSelect(tg, k);
			System.out.println("Interval took:\t" + (System.currentTimeMillis()-st));
			
			// Step 2: ComputeADS
			double bestResult = -1;
			Edge bestIntv = new Edge(-1, -1);
			//TreeSet<Edge> bestEdgeSet = null;
			for (int i = 0; i < k; i++) {
				Edge tunedIntv = new Edge();
				Graph ag = tg.aggregateGraph(topK[i].getS(), topK[i].getT());
				if (ag == null) continue;	// time interval error
				MGraph convAG = ConvertAG.doConvertAg(ag);
				
				//using all 3 optimization techniques
				MGraph mergCAG = MergeCCs.doMergeCCsDTrees(convAG);
				MGraph mstAG = SpanningTree.doMinimumST(mergCAG);
				DenseSubGraph dsg = StrongPrunning.doStrongPrunning(mstAG);
				BoundedProbing.bp(mergCAG, dsg);
				TreeSet<Edge> dsgEdges = new TreeSet<Edge>();
				double finalDSGWeight = SpanningTree.doMinimumST(convAG, BoundedProbing.getUsed(), 
						MergeCCs.getFinalNodeID(), dsgEdges);
				finalDSGWeight = IntervalSelector.tuning(tg, dsgEdges, topK[i], finalDSGWeight, 
						tunedIntv);
				if (finalDSGWeight > bestResult) {
					bestResult = finalDSGWeight;
					bestIntv.copy(tunedIntv);
					//bestEdgeSet = new TreeSet<Edge>(dsgEdges);
				}
				// using Strong Pruning + Bounded Probing		
				MGraph mstAG2 = SpanningTree.doMinimumST(convAG);
				DenseSubGraph dsg2 = StrongPrunning.doStrongPrunning(mstAG2);
				BoundedProbing.bp(convAG, dsg2);
				TreeSet<Edge> dsgEdges2 = new TreeSet<Edge>();
				double  finalDSGWeight2 = SpanningTree.doMinimumST(convAG, 
						BoundedProbing.getUsed(), dsgEdges2);
				finalDSGWeight2 = IntervalSelector.tuning(tg, dsgEdges2, topK[i], finalDSGWeight2, 
						tunedIntv);
				if (finalDSGWeight2 > bestResult) {
					bestResult = finalDSGWeight2;
					bestIntv.copy(tunedIntv);
					//bestEdgeSet = new TreeSet<Edge>(dsgEdges2);
				}
			}
			System.out.println("[" + (bestIntv.getS() + intvL[intv]) + "," + 
					(bestIntv.getT() + intvL[intv]) + "]" + bestResult);
			System.out.println("Evaluation took:\t" + (System.currentTimeMillis()-st));
			System.out.println();
		}
	}
	
	static void computeADS(String option, String srcGraph, int startT, int endT) 
			throws NumberFormatException, IOException {
		System.out.println("FIDES_plus " + option + " " + srcGraph + " " + startT + " " + endT);
		// Step 0: loading temporal graph
		long st = System.currentTimeMillis();
		TGraph tg = new TGraph(srcGraph, startT, endT);
		System.out.println("Loading took:\t" + (System.currentTimeMillis()-st));
		
		int[] intvL = new int[1000];
		int[] intvR = new int[1000];
		int cnt = 0;
		BufferedReader br = new BufferedReader(new FileReader("interval.txt"));
		String line = null;
		char sepChar = ',';
		while((line = br.readLine()) != null) {
			if (line.length() == 0) continue;
			int sep = line.indexOf(sepChar);
			intvL[cnt] = Integer.valueOf(line.substring(0, sep));
			intvR[cnt] = Integer.valueOf(line.substring(sep + 1));
			cnt++;
		}
		br.close();
		
		st = System.currentTimeMillis();
		for (int i = 0; i < cnt; i++) {
			long stInner = System.currentTimeMillis();
			Graph ag = tg.aggregateGraph(intvL[i], intvR[i]);
			if (ag != null) {	// time interval error
				MGraph convAG = ConvertAG.doConvertAg(ag);
				// using all optimizations
				MGraph mergCAG = MergeCCs.doMergeCCsDTrees(convAG, tg.getId2name());
				MGraph mstAG = SpanningTree.doMinimumST(mergCAG);
				DenseSubGraph dsg = StrongPrunning.doStrongPrunning(mstAG);
				//long bpTime = System.currentTimeMillis();
				BoundedProbing.bp(mergCAG, dsg);
				//System.out.print((System.currentTimeMillis()-bpTime) + "ms\t");
				TreeSet<Edge> dsgEdges = new TreeSet<Edge>();
				double finalDSGWeight = SpanningTree.doMinimumST(convAG, 
						BoundedProbing.getUsed(), 
						MergeCCs.getFinalNodeID(), dsgEdges);
				
				// using Strong Pruning + Bounded Probing		
				MGraph mstAG2 = SpanningTree.doMinimumST(convAG);
				DenseSubGraph dsg2 = StrongPrunning.doStrongPrunning(mstAG2);
				//bpTime = System.currentTimeMillis();
				BoundedProbing.bp(convAG, dsg2);
				//System.out.print((System.currentTimeMillis()-bpTime) + "ms\t");
				TreeSet<Edge> dsgEdges2 = new TreeSet<Edge>();
				double finalDSGWeight2 = SpanningTree.doMinimumST(convAG, 
						BoundedProbing.getUsed(), dsgEdges2);
				
				System.out.print(intvL[i] + "\t" + intvR[i] + "\t");
				
//				System.out.print(finalDSGWeight + "\t");
//				System.out.print(ag.computeSumOfEdgeW(dsgEdges) + "check\t");
//				System.out.print(finalDSGWeight2 + "\t");
//				System.out.print(ag.computeSumOfEdgeW(dsgEdges2) + "check\t");
				
				if (finalDSGWeight >= finalDSGWeight2)
					System.out.print(finalDSGWeight + "\t");
				else {
					System.out.print(finalDSGWeight2 + "\t");
					//System.out.print("bad_strong_merging\t");
				}
				//System.out.print(ag.getPosEdgeCount() + "\t");
				System.out.println((System.currentTimeMillis()-stInner));
				//System.out.println();
			}
			
			//break;
		}
		//System.out.println("Avg Evaluation took:\t" + (System.currentTimeMillis()-st) / cnt);
		System.out.println();
	}
	
	static void ADSOPT(String option, String srcGraph, int startT, int endT, String intv) 
			throws NumberFormatException, IOException, InterruptedException {
		System.out.println("FIDES_plus " + option + " " + srcGraph + " " + startT + " " + endT);
		// Step 0: loading temporal graph
		long st = System.currentTimeMillis();
		TGraph tg = new TGraph(srcGraph, startT, endT);
		System.out.println("Loading took:\t" + (System.currentTimeMillis()-st));
		
		int[] intvL = new int[1000];
		int[] intvR = new int[1000];
		int cnt = 0;
		BufferedReader br = new BufferedReader(new FileReader(intv));
		String line = null;
		char sepChar = ',';
		while((line = br.readLine()) != null) {
			if (line.length() == 0) continue;
			int sep = line.indexOf(sepChar);
			intvL[cnt] = Integer.valueOf(line.substring(0, sep));
			intvR[cnt] = Integer.valueOf(line.substring(sep + 1));
			cnt++;
		}
		br.close();
		
		st = System.currentTimeMillis();
		for (int i = 0; i < cnt; i++) {
			Graph ag = tg.aggregateGraph(intvL[i], intvR[i]);
			if (ag != null) {	// time interval error
				System.out.print(intvL[i] + "\t" + intvR[i] + "\t");
				
				// OPT1: Strong Pruning (including an MST as pre-processing)
				long stBegin = System.currentTimeMillis();
				MGraph convAG1 = ConvertAG.doConvertAg(ag);					
				MGraph mstAG1 = SpanningTree.doMinimumST(convAG1);
				DenseSubGraph dsg1 = StrongPrunning.doStrongPrunning(mstAG1);
				System.out.print(dsg1.getDSGWeight() + "\t");
				System.out.print((System.currentTimeMillis()-stBegin) + "\t"); 
				//tg.printSelectEdge(dsg1.getDSGEdges(), intvL[i], intvR[i]);
				convAG1 = mstAG1 = null; dsg1 = null;
				Thread.sleep(100);
				
				// OPT2: Strong Pruning + Bounded Probing
				stBegin = System.currentTimeMillis();
				MGraph convAG2 = ConvertAG.doConvertAg(ag);					
				MGraph mstAG2 = SpanningTree.doMinimumST(convAG2);
				DenseSubGraph dsg2 = StrongPrunning.doStrongPrunning(mstAG2);
				double finalDSGWeight2 = BoundedProbing.bp(convAG2, dsg2);
				TreeSet<Edge> dsgEdges2 = new TreeSet<Edge>();
				finalDSGWeight2 = SpanningTree.doMinimumST(convAG2, BoundedProbing.getUsed(), 
						dsgEdges2);
				System.out.print(finalDSGWeight2 + "\t");
				System.out.print((System.currentTimeMillis()-stBegin) + "\t"); 
				//tg.printSelectEdge(BoundedProbing.getDSGEdges(), intvL[i], intvR[i]);
				//tg.printSelectEdge(dsgEdges2, intvL[i], intvR[i]);
				convAG2 = mstAG2 = null; dsg2 = null;
				Thread.sleep(100);
				
				// OPT3: Strong Merging + Strong Pruning + MST
				stBegin = System.currentTimeMillis();
				MGraph convAG3 = ConvertAG.doConvertAg(ag);	
				MGraph mergCAG3 = MergeCCs.doMergeCCsDTrees(convAG3);
				//mergCAG.testGraph(tg.getId2name());
				MGraph mstAG3 = SpanningTree.doMinimumST(mergCAG3);
				DenseSubGraph dsg3 = StrongPrunning.doStrongPrunning(mstAG3);
				TreeSet<Edge> dsgEdges3 = new TreeSet<Edge>();
				double finalDSGWeight3 = SpanningTree.doMinimumST(convAG3, dsg3.getUsed(), 
						MergeCCs.getFinalNodeID(), dsgEdges3);
				System.out.print(finalDSGWeight3 + "\t");
				System.out.print((System.currentTimeMillis()-stBegin) + "\t");	
				//tg.printSelectEdge(dsgEdges3, intvL[i], intvR[i]);
				convAG3 = mergCAG3 = mstAG3 = null; dsg3 = null;
				dsgEdges3 = null;
				Thread.sleep(100);
				
				// OPT4: all
				stBegin = System.currentTimeMillis();
				MGraph convAG4 = ConvertAG.doConvertAg(ag);	
				MGraph mergCAG4 = MergeCCs.doMergeCCsDTrees(convAG4);
				MGraph mstAG4 = SpanningTree.doMinimumST(mergCAG4);
				DenseSubGraph dsg4 = StrongPrunning.doStrongPrunning(mstAG4);
				BoundedProbing.bp(mergCAG4, dsg4);
				TreeSet<Edge> dsgEdges4 = new TreeSet<Edge>();
				double finalDSGWeight4 = SpanningTree.doMinimumST(convAG4, BoundedProbing.getUsed(), 
						MergeCCs.getFinalNodeID(), dsgEdges4);
				System.out.print(finalDSGWeight4 + "\t");
				System.out.println(System.currentTimeMillis()-stBegin);
				//tg.printSelectEdge(dsgEdges4, intvL[i], intvR[i]);
				Thread.sleep(100);
			}
		}
		System.out.println();
	}
	
	static void ecrate(String option, String srcGraph, int startT, int endT) {
		System.out.println("FIDES_plus " + option + " " + srcGraph + " " + startT + " " + endT);
		long st = System.currentTimeMillis();
		TGraph tg = new TGraph(srcGraph, startT, endT);
		System.out.println("ECRate: " + tg.getECRate());
		System.out.println("Evaluation took:\t" + (System.currentTimeMillis()-st));
	}
}
