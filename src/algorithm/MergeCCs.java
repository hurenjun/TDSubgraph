package algorithm;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import graph.Edge;
import graph.LinkCutTree;
import graph.MGraph;
import graph.WST;

public class MergeCCs {
//	public static class CompEdge implements Comparable<CompEdge>{
//		int cid, s, t;
//		float w;
//		public CompEdge(final int _cid, final int _s, final int _t, final float _w) {
//			cid = _cid; s = _s; t = _t;
//			w = _w;
//		}
//		@Override
//		public int compareTo(CompEdge o) {
//			if (cid < o.cid || (cid == o.cid && s < o.s) || (cid == o.cid && s == o.s && t < o.t))
//				return -1;
//			else if (cid == o.cid && s == o.s && t == o.t) 
//				return 0;
//			else return 1;
//		}
//		
//		public boolean equals(CompEdge o) {
//			if (cid == o.cid && s == o.s && t == o.t) return true;
//			else return false;
//		}
//	}
	
	static int nNode, nEdge;
	static boolean validNode[];
	static int cid[];
	static int edgeCid[];
	static double nodeW[];
	static boolean compChange[];
	static TreeMap<Integer, WST> edgeCC[];	
	static LinkCutTree lct;
	static TreeSet<Edge> selectEdge;
	static int finalNodeID[];
	/* 
	 * member variables descriptions:
	 * 1) nNode & nEdge: number of nodes and edges (remain fixed although nodes will be merged);
	 * 2) validNode[i]: whether node i is valid (not merge or not remove after merging);
	 * 3) cid[i] / edgeCid[j]: component id of node i / edge j;
	 * 4) nodeW[i]: node weight of node i, update after merging;
	 * 5) compChange[i]: whether component i is updated;
	 * 6) edgeCC[i]: edges of connected component i, initially each component contain single node,
	 * 	  edgeCC[i] need update after merging
	 */
	
	public static MGraph doMergeCCs(MGraph convAG) {
		nNode = convAG.getnNode();
		nEdge = convAG.getnEdge();
		final int[] nodeInd = convAG.getnodeInd();
		final int[] edgeT = convAG.getedgeT();
		final int[] oriS = convAG.getOriS();
		final int[] oriT = convAG.getOriT();
		final double w[] = convAG.getw();
		final double[] oriNodeW = convAG.getNodeW();
		edgeCC = new TreeMap[nNode];
		for (int i = 0; i < nNode; i++) {
			edgeCC[i] = new TreeMap<Integer, WST>();
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				edgeCC[i].put(edgeT[j], new WST(w[j], oriS[j], oriT[j]));	
			}
		}
		// initially, no nodes(components) are further merged, so the node weight remains same;
		nodeW = new double[nNode]; convAG.copyNodeW(nodeW);
		validNode = new boolean[nNode];
		Arrays.fill(validNode, true);
		cid = new int[nNode];
		for (int i = 0; i < nNode; i++) {cid[i] = i;}
		edgeCid = new int[nEdge];
		Arrays.fill(edgeCid, -1);
		compChange = new boolean[nNode];
		while (true) {
			Arrays.fill(compChange, false);
			boolean change = merge();			
			if (!change) break;
			updateNodeW(nodeInd, edgeT, w, oriNodeW);
		}
		
		MGraph mergCAG = convAG.merge(cid, edgeCid, nodeW, edgeCC);
		
		clear();
		
		return mergCAG;
	}
	
	// mergeCCs by dynamic trees
	public static MGraph doMergeCCsDTrees(MGraph convAG) {
		nNode = convAG.getnNode();
		nEdge = convAG.getnEdge();
		final int[] nodeInd = convAG.getnodeInd();
		final int[] edgeT = convAG.getedgeT();
		final int[] oriS = convAG.getOriS();
		final int[] oriT = convAG.getOriT();
		final double w[] = convAG.getw();
		nodeW = new double[nNode]; convAG.copyNodeW(nodeW);
		TreeSet<Integer>[] nodesInComp = new TreeSet[nNode];
		edgeCC = new TreeMap[nNode];
		cid = new int[nNode];
		for (int i = 0; i < nNode; i++) {
			cid[i] = i;
			nodesInComp[i] = new TreeSet<Integer>();
			nodesInComp[i].add(i);
			edgeCC[i] = new TreeMap<Integer, WST>();
		}
		lct = new LinkCutTree(nNode);
		selectEdge = new TreeSet<Edge>();
		
		while (true) {
			boolean change = false;
			for (int i = 0; i < nNode; i++) {
				int p = findRoot(i);
				if (nodeW[p] <= 0) continue;
				for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
					int q = findRoot(edgeT[j]);
					if (p >= q) continue;
					if (nodeW[p] >= w[j] && nodeW[q] >= w[j]) {
						mergeComp(p, q, nodesInComp, nodeInd, edgeT, w);
						change = true;
					}
				}
			}
			if (!change) break;
		}
		
		edgeCid = new int[nEdge];
		for (int i = 0; i < nNode; i++) {
			int p = findRoot(i);
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				int q = findRoot(edgeT[j]);
				if (selectEdge.contains(new Edge(i, edgeT[j]))) edgeCid[j] = p;
				else edgeCid[j] = -1;
				if (p < q) {
					if (!edgeCC[p].containsKey(q) || edgeCC[p].get(q).getW() > w[j]) {
						edgeCC[p].put(q, new WST(w[j], oriS[j], oriT[j]));
					}
				}
			}
		}
		MGraph mergCAG = convAG.merge(cid, edgeCid, nodeW, edgeCC);
		clear();
		return mergCAG;
	}
	
	// merge nodes if there exists node pair s and t such that nodeW[s] >= w(s, t)
	// and nodeW[t] >= w(s, t)
	private static boolean merge() {
		boolean change = false;
		for (int i = 0; i < nNode; i++) {
			if (!validNode[i]) continue;
			int s = i;
			Iterator<Entry<Integer, WST>> iterS = edgeCC[s].entrySet().iterator();
			while (iterS.hasNext()) {
				Entry<Integer, WST> entry = (Entry<Integer, WST>)iterS.next();
				int t = entry.getKey();
				if (s >= t || !validNode[t]) continue;
				// here s < t && validNode[t] = true
				double w = entry.getValue().getW();
				if (nodeW[s] >= w && nodeW[t] >= w) {	// merge t to s (always large to small)
					validNode[t] = false;
					cid[t] = s;
					compChange[s] = true;
					change = true;
					nodeW[s] = nodeW[s] + nodeW[t] - w;
					Iterator<Entry<Integer, WST>> iterT = edgeCC[t].entrySet().iterator();
					while (iterT.hasNext()) {
						Entry<Integer, WST> entryT = (Entry<Integer, WST>)iterT.next();
						int nbT = entryT.getKey();
						if (edgeCC[nbT] != null) edgeCC[nbT].remove(t);
						if (!validNode[nbT] || nbT == s) continue;
						double nbW = entryT.getValue().getW();
						if (!edgeCC[s].containsKey(nbT) || nbW < edgeCC[s].get(nbT).getW()) {
							edgeCC[s].put(nbT, new WST(entryT.getValue()));
							edgeCC[nbT].put(s, new WST(entryT.getValue()));
						}
					}
					edgeCC[t] = null;
					iterS = edgeCC[s].entrySet().iterator();
				}
			}
		}
		return change;
	}
	
	private static void updateNodeW(final int[] nodeInd, final int edgeT[], final double w[], 
			final double oriNodeW[]) {
		for (int i = 0; i < nNode; i++) {
			findRoot(i);
			if (compChange[i]) nodeW[i] = 0;	// nodeW of component i need to be updated
		}
		// nodeW of component = sum of original node weight of nodes - sum of edge weights of MST
		for (int i = 0; i < nNode; i++) {		// add up node weight of nodes in each components
			if (compChange[cid[i]]) nodeW[cid[i]] += oriNodeW[i];
		}
		
		// select edges of changed components
		int selectS[] = new int[nEdge / 2];
		int selectT[] = new int[nEdge / 2];
		double selectW[] = new double[nEdge / 2];
		int selectEid[] = new int[nEdge / 2];
		int cnt = 0;
		for (int i = 0; i < nNode; i++) {
			int s = i;
			for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
				int t = edgeT[j];
				if (cid[s] != cid[t] || s > t) continue;	// only remains edges within components
				if (!compChange[cid[s]]) continue;	// remove edges within unchanged components
				// all edges are initially set to not belong to any component
				// select edges of components by MST
				edgeCid[j] = -1; 
				selectS[cnt] = s; selectT[cnt] = t;
				selectW[cnt] = w[j];
				selectEid[cnt] = j;
				cnt++;
			}
		}
		
		// MST using selected edges
		SpanningTree.doMinimumST(nNode, cnt, selectS, selectT, selectW, selectEid);
		for (int i = 0; i < cnt; i++) {
			if (selectEid[i] != -1) {	// this edge is further selected by MST
				edgeCid[selectEid[i]] = cid[selectS[i]];
				nodeW[cid[selectS[i]]] -= selectW[i];	// update nodeW of component
			}
		}
	}
	
	private static void mergeComp(final int p, final int q, TreeSet<Integer> nodesInComp[], 
			final int[] nodeInd, final int edgeT[], final double w[]) {
		Iterator<Integer> iter = nodesInComp[p].size() < nodesInComp[q].size() ? 
				nodesInComp[p].iterator() : nodesInComp[q].iterator();
		int target = nodesInComp[p].size() < nodesInComp[q].size() ? q : p;
		boolean first = true;
		while (iter.hasNext()) {
			int s = (Integer)iter.next();
			for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
				int t = edgeT[j];
				if (findRoot(t) == target) {
					if (first) {	// directly add the edge
						//int root = lct.findRoot(s);
						lct.evert(s);
						lct.link(s, t, w[j]);
						nodeW[p] -= w[j];
						int minST = Math.min(s, t);
						selectEdge.add(new Edge(minST, s + t - minST));
						first = false;
					}
					else {			// verify and update if need
						int child = lct.maximumEdgeVW(s, t);
						if (lct.getMaxEdgeWeight() > w[j]) {
							int parent = lct.getRealParent(child);
							lct.cut(child);
							//int root = lct.findRoot(s);
							lct.evert(s);
							lct.link(s, t, w[j]);
							nodeW[p] = nodeW[p] - w[j] + lct.getMaxEdgeWeight();
							int minCP = Math.min(child, parent);
							if (!selectEdge.remove(new Edge(minCP, child + parent - minCP)))
								System.out.println("not delete");
							int minST = Math.min(s, t);
							selectEdge.add(new Edge(minST, s + t - minST));
						}
					}
				}
			}
		}
		cid[q] = p;
		nodeW[p] += nodeW[q];
		nodeW[q] = 0;
		nodesInComp[p].addAll(nodesInComp[q]);
		nodesInComp[q] = null;
	}
	
	private static int findRoot(int k) {
		if (k == cid[k]) return k;
		cid[k] = findRoot(cid[k]);
		return cid[k];
	}
	
	// delete space used in static methods
	private static void clear() {
		finalNodeID = new int[nNode];
		HashMap<Integer, Integer> cidMap = new HashMap<Integer, Integer>();
		int compCnt = 0;
		for (int i = 0; i < nNode; i++) {
			if (!cidMap.containsKey(cid[i])) cidMap.put(cid[i], compCnt++);
		}
		for (int i = 0; i < nNode; i++) {
			finalNodeID[i] = cidMap.get(cid[i]);
		}
		
		cid = null;
		validNode = null;
		nodeW = null;
		compChange = null;
		edgeCC = null;
		edgeCid = null;
		lct.clear();
		selectEdge = null;
	}
	
	public static final int[] getFinalNodeID() {
		return finalNodeID;
	}
}

