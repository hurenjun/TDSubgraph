package algorithm;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import graph.Edge;
import graph.LinkCutTree;
import graph.MEdge;
import graph.MGraph;
import graph.WST;

public class MergeCCs {	
	static int nNode, nEdge;
	static int[] cid, bestCid;
	static double[] nodeW, bestNodeW;
	static MGraph thisConvAG;
	static HashMap<Integer, Integer> id2name;
	
	static int[] degree;
	static TreeSet<WST>[] acrossEdge;
	static TreeMap<Integer, MEdge>[] minCCedge, bestMinCCedge;	
	static ArrayList<Edge> mergNodePair;
	static LinkCutTree lct;
	static TreeSet<Edge> selectEdge, bestSelectEdge;	// insertion & deletion 
	static int[] finalNodeID;
	static double globalOptima;
	/*
	 * member variables descriptions:
	 * 1) nNode & nEdge: number of nodes and edges (remain fixed although nodes will be merged);
	 * 3) cid[i]: component id of node i;
	 * 4) nodeW[i]: node weight of node i, update after merging;
	 * 5) compChange[i]: whether component i is updated;
	 * 6) edgeCC[i]: edges of connected component i, initially each component contain single node,
	 * 	  edgeCC[i] need update after merging
	 */
	
	public static MGraph doMergeCCsDTrees(MGraph convAG, HashMap<Integer, Integer> map) {
		id2name = map;
		return doMergeCCsDTrees(convAG);
	}
	
	/**
	 * strong merging with dynamic trees
	 * @param convAG converted aggregate graph
	 * @return merged converted graph
	 */
	public static MGraph doMergeCCsDTrees(MGraph convAG) {
		thisConvAG = convAG;
		
		nNode = convAG.getnNode();
		nEdge = convAG.getnEdge();
		final int[] nodeInd = convAG.getnodeInd();
		final int[] edgeT = convAG.getedgeT();
		final int[] oriS = convAG.getOriS();
		final int[] oriT = convAG.getOriT();
		final double[] w = convAG.getw();
		
		nodeW = new double[nNode]; convAG.copyNodeW(nodeW);		
		cid = new int[nNode];
		degree = new int[nNode];
		minCCedge = new TreeMap[nNode];
		mergNodePair = new ArrayList<Edge>();
		mergNodePair.ensureCapacity(nEdge);
		acrossEdge = new TreeSet[nNode];
		for (int i = 0; i < nNode; i++) {
			cid[i] = i;
			degree[i] = nodeInd[i + 1] - nodeInd[i];
			minCCedge[i] = new TreeMap<Integer, MEdge>();
			acrossEdge[i] = new TreeSet<WST>();
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				int t = edgeT[j];
				acrossEdge[i].add(new WST(w[j], i, t));
				minCCedge[i].put(t, new MEdge(i, t, oriS[j], oriT[j], w[j]));
				//minCCedge[i].put(t, new WST(w[j], oriS[j], oriT[j]));
				if (nodeW[i] >= w[j] && nodeW[t] >= w[j] && i < t) 
					mergNodePair.add(new Edge(i, t));
			}
		}
		lct = new LinkCutTree(nNode);
		selectEdge = new TreeSet<Edge>();
		
		int left = 0;
		while (left < mergNodePair.size()) {
			int p = findRoot(mergNodePair.get(left).getS());
			int q = findRoot(mergNodePair.get(left).getT());
			if (p != q)	{
				double wp = nodeW[p], wq = nodeW[q];
				mergeNodeWeight(p, q);
				mergeEdges(p, q, wp, wq);
			}
			
			left++;
		}
		
		// TODO: node optimization
		globalOptima = Double.NEGATIVE_INFINITY;
		nodeOptimization();
		
		//cid = bestCid; selectEdge = bestSelectEdge;
		//nodeW = bestNodeW; minCCedge = bestMinCCedge;
		
		int[] edgeCid = new int[nEdge];
		for (int i = 0; i < nNode; i++) {
			int p = findRoot(i);
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				if (selectEdge.contains(new Edge(i, edgeT[j]))) edgeCid[j] = p;
				else edgeCid[j] = -1;
			}
		}
		MGraph mergCAG = convAG.merge(cid, edgeCid, nodeW, minCCedge);
		//printSelectEdge();
		clear();
		return mergCAG;
	}
	
	
	/**
	 * Node optimization: directly merging two nodes using two links (one for linking and 
	 * the other replacing an existing link within nodes) if it increases node weights of 
	 * both nodes.
	 */
	private static void nodeOptimization() {
		// construct the potential node_opt_links
		ArrayList<NodeOptLink> noptLink = new ArrayList<NodeOptLink>(nEdge);
		for (int i = 0; i < nNode; i++) {
			if (findRoot(i) == i) {
				double nw1 = nodeW[i];
				for (WST wst: acrossEdge[i]) {
					int q = findRoot(wst.getT());
					if (q != i) {
						double nw2 = nodeW[q];
						noptLink.add(new NodeOptLink(wst, nw1, nw2));
					}
				}
			}
		}
		
		// node optimization using three heuristics
		//TODO
		//nodeOptGivenComp(noptLink, NodeOptLink.getLEWComparator());
		nodeOptGivenComp(noptLink, NodeOptLink.getHNWComparator());
		//nodeOptGivenComp(noptLink, NodeOptLink.getHNWLEWComparator());
	}
	
	/**
	 * Three heuristics will be exploited. And we choose the one giving the best results.
	 * (1) Heuristics Least Edge Weight: always merging using the link of the least weight. 
	 * Note that the other link is the one of least weight between the two nodes to be
	 * merged.   
	 * (2) Heuristics Highest Node Weight: always merging the two nodes of highest node weights, 
	 * i.e., finding the highest one first, and then finding the highest one connecting to 
	 * the first node.
	 * (3) Heuristics Highest Node Weight Least Edge Weight: always merging the node of the highest 
	 * node weight with another using the link of the least weight.  
	 */
	private static void nodeOptGivenComp(ArrayList<NodeOptLink> noptLink, 
			Comparator<NodeOptLink> comp) {
		// data backup
//		int[] cidCopy = new int[nNode]; 
//		System.arraycopy(cid, 0, cidCopy,0, nNode);
//		double[] nodeWCopy = new double[nNode];
//		System.arraycopy(nodeW, 0, nodeWCopy,0, nNode);
//		int[] degreeCopy = new int[nNode];
//		System.arraycopy(degree, 0, degreeCopy,0, nNode);
//		TreeSet<WST>[] acrossEdgeCopy = new TreeSet[nNode];
//		TreeMap<Integer, MEdge>[] minCCedgeCopy = new TreeMap[nNode];	
//		for (int i = 0; i < nNode; i++) {
//			acrossEdgeCopy[i] = new TreeSet<WST>(acrossEdge[i]);
//			minCCedgeCopy[i] = new TreeMap<Integer, MEdge>(minCCedge[i]);
//		}
//		LinkCutTree lctCopy = new LinkCutTree(lct);
//		TreeSet<Edge> selectEdgeCopy = new TreeSet<Edge>(selectEdge); 
		
		// sort node optimization links, and execute one by one
		noptLink.sort(comp);
//		int noptCount = 0;
//		int oriSize = mergNodePair.size();
		int left = mergNodePair.size();
		for (int index = 0; index < noptLink.size(); index++) {
			int s = noptLink.get(index).s;
			int t = noptLink.get(index).t;
			double w = noptLink.get(index).w;
			int p = findRoot(s);
			int q = findRoot(t);
			if (p != q) {
				MEdge me = minCCedge[p].get(q);
				if (!me.equalST(s, t)) {
					double wRep = -1;
					if (p == findRoot(me.getS())) {	
						// (s & wst.s)/(t & wst.t) belong to the same trees
						lct.maximumEdgeVW(s, me.getS());
						wRep = lct.getMaxEdgeWeight();
						lct.maximumEdgeVW(t, me.getT());
						wRep = Math.max(wRep, lct.getMaxEdgeWeight());
					} 
					else {
						// (s & wst.t)/(t & wst.s) belong to the same trees
						lct.maximumEdgeVW(s, me.getT());
						wRep = lct.getMaxEdgeWeight();
						lct.maximumEdgeVW(t, me.getS());
						wRep = Math.max(wRep, lct.getMaxEdgeWeight());
					}
					double cost = w + me.getW() - wRep; 
					if (nodeW[p] >= cost && nodeW[q] >= cost) {
						//logNodeOpt(noptLink.get(index), me, wRep);
						double wp = nodeW[p], wq = nodeW[q];
						mergeNodeWeight(p, q);
						mergeEdges(p, q, wp, wq);
						//noptCount++;
					}
				}
			}
			
			// if new mergeable node pairs exist
			while (left < mergNodePair.size()) {
				p = findRoot(mergNodePair.get(left).getS());
				q = findRoot(mergNodePair.get(left).getT());
				if (p != q)	{
					double wp = nodeW[p], wq = nodeW[q];
					mergeNodeWeight(p, q);
					mergeEdges(p, q, wp, wq);
				}
				left++;
			}
		}
		
//		while (mergNodePair.size() > oriSize) {
//			mergNodePair.remove(mergNodePair.size() - 1);
//		}
		
		// update best solution according to max_node_weight
//		double maxNW = 0;
//		for (int i = 0; i < nNode; i++) {
//			if (nodeW[i] > maxNW) 
//				maxNW = nodeW[i];
//		}
//		if (noptCount > globalOptima) {
//			bestCid = cid;
//			bestSelectEdge = selectEdge;
//			bestNodeW = nodeW;
//			bestMinCCedge = minCCedge; 
//			globalOptima = noptCount;
//		}
		
		// data recover
//		cid = cidCopy;
//		nodeW = nodeWCopy;
//		degree = degreeCopy;
//		acrossEdge = acrossEdgeCopy;
//		minCCedge = minCCedgeCopy;
//		lct = lctCopy;
//		selectEdge = selectEdgeCopy;
	}
	
	/**
	 * Merge nodes p and q which are guaranteed to be merge-able. It uses all
	 * edges of the lower-degree node. Since the final merged node maintains
	 * a MST, it uses all across edges between p and q to optimize the MST, by
	 * leveraging Link-Cut-Tree in O(log n) time for each across edge. Edges 
	 * will be merged by another procedure (e.g., mergeEdges()).
	 * @param p node to be merged
	 * @param q another node to be merged
	 */
	private static void mergeNodeWeight(final int p, final int q) {
		int highDeg = degree[p] <= degree[q] ? q : p;
		boolean first = true;
		for (WST wst: acrossEdge[p + q - highDeg]) {
			int s = wst.getS();	
			int t = wst.getT();
			if (findRoot(t) == highDeg) {
				if (first) {	// directly add the edge
					//int root = lct.findRoot(s);
					lct.evert(s);
					lct.link(s, t, wst.getW());
					nodeW[highDeg] -= wst.getW();
					int minST = Math.min(s, t);
					selectEdge.add(new Edge(minST, s + t - minST));
					first = false;
				}
				else {			// verify and update if need
					int child = lct.maximumEdgeVW(s, t);
					if (lct.getMaxEdgeWeight() > wst.getW()) {
						int parent = lct.getRealParent(child);
						lct.cut(child);
						//int root = lct.findRoot(s);
						lct.evert(s);
						lct.link(s, t, wst.getW());
						nodeW[highDeg] = nodeW[highDeg] - wst.getW() + lct.getMaxEdgeWeight();
						int minCP = Math.min(child, parent);
						if (!selectEdge.remove(new Edge(minCP, child + parent - minCP)))
							System.out.println("not delete");
						int minST = Math.min(s, t);
						selectEdge.add(new Edge(minST, s + t - minST));
					}
				}
			}
		}
		
		// merging lower degree to higher degree
		cid[p + q - highDeg] = highDeg;
		nodeW[highDeg] += nodeW[p + q - highDeg];
		degree[highDeg] += degree[p + q - highDeg];
		nodeW[p + q - highDeg] = 0;
		degree[p + q - highDeg] = 0;
	}
	
	/**
	 * Probe newly mergeable node pairs and relink edges of lower-degree node to
	 * the higher-degree node.
	 * @param p root of one node
	 * @param q root of another node
	 * @param wp previous weight of node p
	 * @param wq previous weight of node q
	 * @param mergNodePair list of mergeable node pairs
	 */
	private static void mergeEdges(final int p, final int q, final double wp,
			final double wq) {		
		double wpq = (degree[p] == 0) ? nodeW[q] : nodeW[p];
		// subsetP: > wp && <= wpq
		TreeSet<WST> subsetP = (TreeSet<WST>)acrossEdge[p].subSet( 
				new WST(wp, nNode, nNode), new WST(wpq, nNode, nNode));	
		for (WST wst: subsetP) {
			int r1 = findRoot(wst.getS()), r2 = findRoot(wst.getT());
			if (r1 != r2 && nodeW[r1] >= wst.getW() && nodeW[r2] >= wst.getW())
				mergNodePair.add(new Edge(wst.getS(), wst.getT()));
		}
		// subsetQ: > wq && <= wpq
		TreeSet<WST> subsetQ = (TreeSet<WST>)acrossEdge[q].subSet(
				new WST(wq, nNode, nNode), new WST(wpq, nNode, nNode));
		for (WST wst: subsetQ) {
			int r1 = findRoot(wst.getS()), r2 = findRoot(wst.getT());
			if (r1 != r2 && nodeW[r1] >= wst.getW() && nodeW[r2] >= wst.getW())
				mergNodePair.add(new Edge(wst.getS(), wst.getT()));
		}
		
		int del = (degree[p] == 0) ? p : q;
		for (WST wst: acrossEdge[del]) {
			if (findRoot(wst.getS()) != findRoot(wst.getT())) 
				acrossEdge[p + q - del].add(wst);
		}
		TreeMap<Integer, MEdge> mapTarget = minCCedge[p + q - del];
		for (Entry<Integer, MEdge> entry: minCCedge[del].entrySet()) {
			Integer key = entry.getKey();		
			MEdge value = entry.getValue();
			if (!mapTarget.containsKey(key) || mapTarget.get(key).getW() > value.getW()) {
				mapTarget.put(key, new MEdge(value));
				minCCedge[key].put(p + q - del, new MEdge(value));
			}
			minCCedge[key].remove(del);
		}
		mapTarget.remove(p); mapTarget.remove(q);
		acrossEdge[del].clear();
		minCCedge[del].clear();
	}
	
	private static int findRoot(int k) {
		if (k == cid[k]) return k;
		cid[k] = findRoot(cid[k]);
		return cid[k];
	}
	
	private static void logNodeOpt(final NodeOptLink nol, final MEdge me, final double wRep) {
		if (id2name == null) 
			return;
		
		final int[] nodeInd = thisConvAG.getnodeInd();
		final int[] edgeT = thisConvAG.getedgeT();
		final int[] oriS = thisConvAG.getOriS();
		final int[] oriT = thisConvAG.getOriT();
		System.out.print("node opt: ");
		for (int j = nodeInd[nol.s]; j < nodeInd[nol.s + 1]; j++) {
			if (edgeT[j] == nol.t) {
				System.out.print(id2name.get(oriS[j]) + " " + id2name.get(oriT[j]) 
						+ " " + nol.w + "\t");
			}
		}
		for (int j = nodeInd[me.getS()]; j < nodeInd[me.getS() + 1]; j++) {
			if (edgeT[j] == me.getT()) {
				System.out.print(id2name.get(oriS[j]) + " " + id2name.get(oriT[j])
						+ " " + me.getW() + "\t");
			}
		}
		System.out.println(wRep);
	}
	
	private static void printSelectEdge() {
		final int[] nodeInd = thisConvAG.getnodeInd();
		final int[] edgeT = thisConvAG.getedgeT();
		final int[] oriS = thisConvAG.getOriS();
		final int[] oriT = thisConvAG.getOriT();
		for (Edge e: selectEdge) {
			int s = e.getS(), t = e.getT();
			for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
				if (edgeT[j] == t) 
					System.out.print(id2name.get(oriS[j]) + "," + id2name.get(oriT[j]) + "\t");
			}
		}
		System.out.println();
	}
	
	/**
	 *  delete space used in static methods
	 */
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
		degree = null;
		nodeW = null;
		lct.clear();
		selectEdge = null;
		mergNodePair = null;
		thisConvAG = null;
		id2name = null;
		
		for (int i = 0; i < nNode; i++) {
			minCCedge[i].clear();
			acrossEdge[i].clear();
		}
		minCCedge = null;
		acrossEdge = null;
	}
	
	public static final int[] getFinalNodeID() {
		return finalNodeID;
	}
}

class NodeOptLink{
	int s, t;
	double w, nw1, nw2; // guarantee nw1 >= nw2;
	
	public NodeOptLink(final int s, final int t, final double w, 
			final double nw1, final double nw2) {
		this.s = s; this.t = t;
		this.w = w;
		this.nw1 = (nw1 >= nw2) ? nw1 : nw2;
		this.nw2 = (nw1 < nw2) ? nw1 : nw2;
	}
	
	public NodeOptLink(final WST wst, final double nw1, final double nw2) {
		this.s = wst.getS(); this.t = wst.getT();
		this.w = wst.getW();
		this.nw1 = (nw1 >= nw2) ? nw1 : nw2;
		this.nw2 = (nw1 < nw2) ? nw1 : nw2;
	}
	
	/**
	 * edges of lower edge weight come earlier.
	 */
	static Comparator<NodeOptLink> getLEWComparator() {
		return new Comparator<NodeOptLink>() {
			@Override
			public int compare(NodeOptLink arg0, NodeOptLink arg1) {
				if (arg0.w < arg1.w) return -1;
				else if (arg0.w > arg1.w) return 1;
				return 0;
			}
		};
	}
	
	/**
	 * edges between nodes of higher node weights come earlier.
	 */
	static Comparator<NodeOptLink> getHNWComparator() {
		return new Comparator<NodeOptLink>() {
			@Override
			public int compare(NodeOptLink arg0, NodeOptLink arg1) {
				if (arg0.nw1 > arg1.nw1 || arg0.nw1 == arg1.nw1 && arg0.nw2 > arg1.nw2)
					return -1;
				else if (arg0.nw1 == arg1.nw1 && arg0.nw2 == arg1.nw2)
					return 0;
				else return 1;
			}
		};
	}
	
	/**
	 * edges connecting to nodes of higher node weights and of lower edge weights come earlier.
	 */
	static Comparator<NodeOptLink> getHNWLEWComparator() {
		return new Comparator<NodeOptLink>() {
			@Override
			public int compare(NodeOptLink arg0, NodeOptLink arg1) {
				if (arg0.nw1 > arg1.nw1 || arg0.nw1 == arg1.nw1 && arg0.w < arg1.w) 
					return -1;
				else if (arg0.nw1 == arg1.nw1 && arg0.w == arg1.w) return 0;
				return 1;
			}
		};
	}
}
