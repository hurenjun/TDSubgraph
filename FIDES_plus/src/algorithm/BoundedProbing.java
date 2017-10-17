package algorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.TreeSet;

import graph.DenseSubGraph;
import graph.Edge;
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
	static boolean[] inSeq;
	static double[] score;
	static boolean[] preVisit, visit;
	static ExpandingPath[] prePath, path, bestPath;
	static final int MAX_DEPTH = 4;
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
		DSGEdges = new TreeSet<Edge>(dsg.listDSGEdges());
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
		
		//long st = System.currentTimeMillis();
		preVisit = new boolean[nNode];
		prePath = new ExpandingPath[nNode];
		visit = new boolean[nNode];
		path = new ExpandingPath[nNode];
		bestPath = new ExpandingPath[nNode];
		//System.out.print((System.currentTimeMillis()-st) + "ms\t");
		for (int i = 0; i < MAX_DEPTH; i++) {
			if (bfsBP() == false) break;
		}
		//System.out.println((System.currentTimeMillis()-st) + "ms\t");
		
		clear();
		return DSGWeight;
	}
	
	private static boolean bfsBP() {
		// long st = System.currentTimeMillis();
		int alSize = mergCAG.getnEdge() * MAX_DEPTH * MAX_DEPTH;
		ArrayList<ExpandingPath> expand = new ArrayList<ExpandingPath>(alSize);
		Arrays.fill(bestPath, null);
		
		// initialize directly reachable nodes by finding the shortest paths
		int[] cheapEToV = new int[nNode];
		Arrays.fill(cheapEToV, -1);
		for (int s = 0; s < nNode; s++) {
			if (used[s]) {
				for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
					int t = edgeT[j];
					if (t != s && !used[t] && (cheapEToV[t] == -1 || w[j] < w[cheapEToV[t]])) 
						cheapEToV[t] = j;
				}
			}
		}
		Arrays.fill(preVisit, false);
		for (int t = 0; t < nNode; t++) {
			if (cheapEToV[t] != -1) {
				double profit = nodeW[t] - w[cheapEToV[t]];
				ExpandingPath ep = new ExpandingPath(cheapEToV[t], profit);
				prePath[t] = ep;
				preVisit[t] = true;
				if (bestPath[t] == null || profit > bestPath[t].profit) {
					bestPath[t] = ep;
				}
			}
		}
		
		// explore farther nodes iteratively
		for (int l = 1; l < MAX_DEPTH; l++) {	// l: length of the previous path
			Arrays.fill(visit, false);
			Arrays.fill(path, null);
			for (int s = 0; s < nNode; s++) {
				if (!preVisit[s]) continue;
				for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
					int t = edgeT[j];
					if (used[t]) continue;
					ExpandingPath ep = prePath[s];
					boolean overlap = false;
					for (int k = 0; k < ep.length; k++) {
						if (edgeT[ep.pathEid[k]] == t) overlap = true;
					}
					if (overlap) continue;
					
					double profit = ep.profit + nodeW[t] - w[j];
					int[] pathEid = new int[l + 1];
					for (int k = 0; k < l; k++) {
						pathEid[k] = ep.pathEid[k];
					}
					pathEid[l] = j;
					ExpandingPath epNow = new ExpandingPath(pathEid, l + 1, profit);
					
					if (path[t] == null || profit > path[t].profit) {
						visit[t] = true;
						path[t] = epNow;
					}
					if (bestPath[t] == null || profit > bestPath[t].profit) {
						bestPath[t] = epNow;
					}
//					if (profit >= 0 && nodeW[t] >= 0) {
//						expand.add(epNow);
//					}
				}
			}
			boolean[] tempVisit = preVisit;
			preVisit = visit;
			visit = tempVisit;
			ExpandingPath[] tempPath = prePath;
			prePath = path;
			path = tempPath;
		}
		//System.out.print((System.currentTimeMillis()-st) + "ms\t");
		//st = System.currentTimeMillis();
		
		// expand subgraph using greedy strategy, i.e. expanding paths according to gains (large first)
		// only consider nodes with positive gain and positive nodeW
		boolean change = false;
		
		for (int s = 0; s < nNode; s++) {
			if (bestPath[s] != null && bestPath[s].profit > 0 && nodeW[s] > 0) {
				expand.add(bestPath[s]);
			}
		}
		
		Collections.sort(expand);
		for (ExpandingPath ep: expand) {
			boolean overlap = false;
			for (int k = 0; k < ep.length; k++) {
				if (used[edgeT[ep.pathEid[k]]])
					overlap = true;
			}
			if (!overlap) {
				for (int k = 0; k < ep.length; k++) {
					used[edgeT[ep.pathEid[k]]] = true;
					DSGEdges.add(new Edge(oriS[ep.pathEid[k]], oriT[ep.pathEid[k]]));
					DSGEdges.addAll(mergCAG.getEdgeInNodeByID(edgeT[ep.pathEid[k]]));
					DSGWeight += ep.profit;
					change = true;
				}
			}
		}
		//System.out.println((System.currentTimeMillis()-st) + "ms");
		//st = System.currentTimeMillis();
		
		return change;
	}
	
	private static void clear() {
		//used = null;
		inSeq = null;
		score = null;
		preVisit = visit = null;
		prePath = path = bestPath = null;
	}
	
	public static TreeSet<Edge> getDSGEdges() {
		return DSGEdges;
		//return null;
	}
	
	public static final boolean[] getUsed() {
		return used;
	}
}

class ExpandingPath implements Comparable<ExpandingPath> {
	int length; 
	double profit;
	int[] pathEid;
	public ExpandingPath(final int[] pathEid, final int length, final double profit) {
		this.profit = profit;
		this.length = length;
		this.pathEid = new int[pathEid.length];
		for (int i = 0; i < pathEid.length; i++) {
			this.pathEid[i] = pathEid[i];
		}
	}
	public ExpandingPath(final int eid, final double profit) {
		this.profit = profit;
		this.length = 1;
		this.pathEid = new int[1];
		this.pathEid[0] = eid;
	}
	
	@Override
	public int compareTo(ExpandingPath o) {
		if (profit > o.profit) return -1;	// place higher profit ahead
		if (profit < o.profit) return 1;
		if (length > o.length) return -1;	// place longer length ahead
		if (length < o.length) return 1;
		for (int i = 0; i < length; i++) {
			if (pathEid[i] < o.pathEid[i]) return -1;
			if (pathEid[i] > o.pathEid[i]) return 1;
		}
		return 0;
	}
	
}
