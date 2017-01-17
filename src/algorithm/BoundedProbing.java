package algorithm;

import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import graph.DenseSubGraph;
import graph.Edge;
import graph.MEdge;
import graph.MGraph;

public class BoundedProbing {
	static MGraph mergCAG;
	static int nNode;
	static int[] nodeInd;
	static int[] edgeT;
	static int[] oriS;
	static int[] oriT;
	static double[] w;
	static double[] nodeW;
	static boolean[] used;
	static MEdge[] seq;
	static boolean[] inSeq;
	static double[] score;
	static boolean[] preVisit, visit;
	static int[][] prePath, path, bestPath;
	static double[] preGain, gain, bestGain;
	static int maxDepth = 4;
	static TreeSet<Edge> DSGEdges;
	static double DSGWeight;
	/*
	 * member variables description:
	 * 1) used[i]: whether node i of expandCAG is used in some expanded path;
	 * 2) seq: sequence of nodes of current expanded path;
	 * 3) maxDepth: maxDepth of expanded path;
	 */
	
	/*
	 * 1) Algorithm Bounded Proving first merges current dense subgraph as a single node (node 0) 
	 *    based on mergCAG, current dense subgraph is a sub-tree of mergCAG. 
	 * 2) Then it tries to expand the dense subgraph by expanded paths. 
	 * 3) A valid expanded path is a sequence of nodes, starting at 0, and the sum of 
	 *    node weights except node 0 is greater than the sum of edge weights;
	 * 4) For any valid expanded path, select the node with maximum score, here score is 
	 *    defined as sum of node weight except node 0 minus sum of edge weight;
	 * 5) The depth of expanded path is controlled by maxDepth; 
	 */
	public static double bp(final MGraph _mergCAG, final DenseSubGraph dsg) {
		mergCAG = _mergCAG;
		DSGEdges = new TreeSet<Edge>(dsg.getDSGEdges());
		DSGWeight = dsg.getDSGWeight();
		
		nNode = mergCAG.getnNode();
		nodeInd = mergCAG.getnodeInd();
		edgeT = mergCAG.getedgeT();
		oriS = mergCAG.getOriS();
		oriT = mergCAG.getOriT();
		w = mergCAG.getw();
		nodeW = mergCAG.getNodeW();
		
		final boolean[] usedDSG = dsg.getUsed();
		used = new boolean[nNode];
		for (int i = 0; i < nNode; i++) { used[i] = usedDSG[i]; }
		
		//DFSBP();
		
		//long st = System.currentTimeMillis();
		preVisit = new boolean[nNode];
		prePath = new int[nNode][maxDepth];
		preGain = new double[nNode];
		visit = new boolean[nNode];
		path = new int[nNode][maxDepth];
		gain = new double[nNode];
		bestGain = new double[nNode];
		bestPath = new int[nNode][maxDepth];
		//System.out.print((System.currentTimeMillis()-st) + "ms\t");
		for (int i = 0; i < maxDepth; i++) {
			//st = System.currentTimeMillis();
			if (bfsBP() == false) break;
			//System.out.print((System.currentTimeMillis()-st) + "ms\t");
		}
		//System.out.println((System.currentTimeMillis()-st) + "ms\t");
		
		clear();
		return DSGWeight;
	}
	
	private static void dfsBP() {
		double[] cost = new double[nNode];
		TreeMap<Integer, MEdge> ohEdges = oneHopEdges(cost);
		
		seq = new MEdge[maxDepth];
		inSeq = new boolean[nNode];
		Arrays.fill(inSeq, false);
		score = new double[maxDepth];
		Iterator<Integer> iter = ohEdges.keySet().iterator();
		while (iter.hasNext()) {
			Integer t = (Integer)iter.next();
			if (used[t]) continue;
			seq[0] = new MEdge(ohEdges.get(t));
			inSeq[t] = true;
			score[0] = nodeW[t] - cost[t];
			dfs(1);
			inSeq[t] = false;
		}
	}
	
	private static boolean bfsBP() {
		// long st = System.currentTimeMillis();
		
		Arrays.fill(bestGain, Double.NEGATIVE_INFINITY);
		Arrays.fill(preVisit, false);
		for (int i = 0; i < nNode; i++) {
			if (used[i]) {
				preVisit[i] = true; preGain[i] = 0;
			}
		}
		
		for (int l = 0; l < maxDepth; l++) {	// l: length of the previous path
			Arrays.fill(gain, Double.NEGATIVE_INFINITY);
			Arrays.fill(visit, false);
			for (int i = 0; i < nNode; i++) {
				if (!preVisit[i]) continue;
				for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
					int t = edgeT[j];
					if (used[t]) continue;
					boolean overlap = false;
					for (int k = 0; k < l; k++) {
						if (edgeT[prePath[i][k]] == t) overlap = true;
					}
					if (overlap) continue;
					if (preGain[i] + nodeW[t] - w[j] > gain[t]) {
						visit[t] = true;
						gain[t] = preGain[i] + nodeW[t] - w[j];
						for (int k = 0; k < l; k++) {
							path[t][k] = prePath[i][k];
						}
						path[t][l] = j;
						if (gain[t] > bestGain[t]) {
							bestGain[t] = gain[t];
							for (int k = 0; k <= l; k++) {
								bestPath[t][k] = path[t][k];
							}
							for (int k = l + 1; k < maxDepth; k++) {
								bestPath[t][k] = -1;
							}
						}
					}
				}
			}
			boolean[] tempVisit = preVisit;
			preVisit = visit;
			visit = tempVisit;
			int[][] tempPath = prePath;
			prePath = path;
			path = tempPath;
			double[] tempGain = preGain;
			preGain = gain;
			gain = tempGain;
		}
		//System.out.print((System.currentTimeMillis()-st) + "ms\t");
		//st = System.currentTimeMillis();
		
		// expand subgraph using greedy strategy, i.e. expanding paths according to gains (large first)
		// only consider nodes with positive gain and positive nodeW
		int[] sortedGainInd = new int[nNode];
		int sortedGainCnt = 0;
		for (int i = 0; i < nNode; i++) {
			if (!used[i] && bestGain[i] > 0 && nodeW[i] > 0) {
				sortedGainInd[sortedGainCnt++] = i;
			}
		}
		quickSort(0, sortedGainCnt - 1, sortedGainInd, bestGain);
		
		boolean change = false;
		for (int i = 0; i < sortedGainCnt; i++) {
			boolean overlap = false;
			int t = sortedGainInd[i];
			for (int k = 0; k < maxDepth; k++) {
				if (bestPath[t][k] == -1) break;
				if (used[edgeT[bestPath[t][k]]]) {
					overlap = true;
					break;
				}
			}
			if (!overlap) {
				for (int k = 0; k < maxDepth; k++) {
					if (bestPath[t][k] == -1) break;
					used[edgeT[bestPath[t][k]]] = true;
					DSGEdges.add(new Edge(oriS[bestPath[t][k]], oriT[bestPath[t][k]]));
					DSGEdges.addAll(mergCAG.getEdgeInNodeByID(edgeT[bestPath[t][k]]));
				}
				DSGWeight += bestGain[t];
				change = true;
			}
		}
		//System.out.println((System.currentTimeMillis()-st) + "ms");
		//st = System.currentTimeMillis();
		
		return change;
	}
	
	private static void quickSort(final int l, final int r, int[] ind, final double[] gain) {
		if (l >= r) return;
		int x = l, y = r;
		int base = ind[l];
		while (x < y) {
			while (x < y && gain[ind[y]] < gain[base]) y--;
			if (x < y) {
				ind[x] = ind[y];
				x++;
			}
			while (x < y && gain[ind[x]] > gain[base]) x++;
			if (x < y) {
				ind[y] = ind[x];
				y--;
			}
		}
		ind[x] = base;
		if (l < x - 1) quickSort(l, x - 1, ind, gain);
		if (x + 1 < r) quickSort(x + 1, r, ind, gain);
	}
	
	private static TreeMap<Integer, MEdge> oneHopEdges(double[] cost) {
		TreeMap<Integer, MEdge> ohEdges = new TreeMap<Integer, MEdge>();
		Arrays.fill(cost, 0);
		for (int i = 0; i < nNode; i++) {
			if (!used[i]) continue;
			int s = i;
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				int t = edgeT[j];
				if (used[t]) continue;
				
				if (!ohEdges.containsKey(t) || cost[t] > w[j]) {
					MEdge me = new MEdge(s, t, oriS[j], oriT[j]);
					ohEdges.put(t, me);
					cost[t] = w[j];
				}
			}
		}
		//System.out.println("One hop edge size: " + ohEdges.size());
		return ohEdges;
	}
	
	private static void dfs(final int depth) {
		if (depth == maxDepth) {
			int best = -1;
			for (int i = 0; i < depth; i++) {
				if (used[seq[i].getT()]) break;
				if (best == -1 || score[i] >= score[best]) best = i;
			}
			if (best != -1 && score[best] > 0) {
				DSGWeight += score[best];
				TreeSet<Edge> newEdges = new TreeSet<Edge>();
				for (int i = 0; i <= best; i++) {
					used[seq[i].getT()] = true;
					DSGEdges.add(new Edge(seq[i].getOriS(), seq[i].getOriT()));
					DSGEdges.addAll(mergCAG.getEdgeInNodeByID(seq[i].getT()));
				}
				//System.out.print(score[best] + " " + ag.GetSumOfEdgeW(newEdges) + "\t");
			}
			return;
		}
		
		int s = seq[depth - 1].getT();
		boolean expanded = false;
		for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
			if (used[seq[0].getT()]) return;
			int t = edgeT[j];
			if (inSeq[t] || used[t]) continue;
			seq[depth] = new MEdge(s, t, oriS[j], oriT[j]);
			inSeq[t] = true;
			score[depth] = score[depth - 1] + nodeW[t] - w[j];
			dfs(depth + 1);
			inSeq[t] = false;
			expanded = true;
		}
		
		if (!expanded) {
			int best = -1;
			for (int i = 0; i < depth; i++) {
				if (used[seq[i].getT()]) break;
				if (best == -1 || score[i] >= score[best]) best = i;
			}
			if (best != -1 && score[best] > 0) {
				DSGWeight += score[best];
				for (int i = 0; i <= best; i++) {
					used[seq[i].getT()] = true;
					DSGEdges.add(new Edge(seq[i].getOriS(), seq[i].getOriT()));
					DSGEdges.addAll(mergCAG.getEdgeInNodeByID(seq[i].getT()));
				}
			}
		}
	}
	
	private static void clear() {
		//used = null;
		seq = null;
		inSeq = null;
		score = null;
		preVisit = visit = null;
		prePath = path = bestPath = null;
		preGain = gain = bestGain = null;
	}
	
	public static TreeSet<Edge> getDSGEdges() {
		return DSGEdges;
		//return null;
	}
	
	public static final boolean[] getUsed() {
		return used;
	}
	
	// get sum of weights of nodes which are not selected
	public static double getRemainNodeW() {
		final double[] nodeW = mergCAG.getNodeW();
		double result = 0;
		for (int i = 0; i < nNode; i++) {
			if (!used[i]) result += nodeW[i];
		}
		return result;
	}
	
	public static double getUsedNodeW() {
		final double[] nodeW = mergCAG.getNodeW();
		double result = 0;
		for (int i = 0; i < nNode; i++) {
			if (used[i]) result += nodeW[i];
		}
		return result;
	}
	
	public static double getAllNodeW() {
		final double[] nodeW = mergCAG.getNodeW();
		double result = 0;
		for (int i = 0; i < nNode; i++) {
			result += nodeW[i];
		}
		return result;
	}
}
