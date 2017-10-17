package algorithm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import graph.Edge;
import graph.Graph;
import graph.MEdge;
import graph.MGraph;

public class ConvertAG {
	static int cid[];
	static double p[];
	static int edgeCid[];
	static ArrayList<MEdge> convAgW;
	static int nNode, nEdge;
	/*
	 * member variable description：
	 * 1) cid[i]: component id of node i;
	 * 2) p[i]: total positive edge weight of component i;
	 * 3) edgeCid[i]: component id of edge i, -1 if the edge does not belong to any component;
	 * 4) convAgW: edges of converted Aggregate graph, i.e. edges connecting different components;
	 *    In additional to two end points of an edge, MEdge also store the original node id of the 
	 *    edge.
	 */
	
	/**
	 * doConvertAG: using only positive edges, merge each connected component into a 
	 * 	single point, whose node weight is the sum of all positive edge weights within
	 * 	the connected component. And the weights of edges between these merged nodes, i.e. 
	 * 	connected components, are minimum absolute value of original edges of ag from node
	 *  of one component to node of another component.
	 * 
	 * Example: input ag has 4 nodes a, b, c and d, with edges: 
	 * 	(a, b) = 2, (c, d) = 2, (a, c) = -2, (b, d) = -4;
	 *  Thus the converted aggregate graph convAG has two nodes A and B, where A is connected
	 *  component of nodes a and b, B is the connected component of nodes c, and d.
	 *  The edge between A and B is selected from all edges from original nodes of A to original
	 *  nodes of B, selecting the one has minimum absolute value, e.g. 2 in this example.
	 * 
	 * The output of doConvertAg is a MGraph (merge graph), which extends from Graph, i.e. having
	 * 	same member variables of Graph. It additional store the node weights of nodes, the original
	 * 	nodes and edges within each new merged nodes, and the original original end-points of new 
	 * 	edges. 
	 * In the previous example, the node weight of A is 4, original nodes and edges of A is nodes a 
	 * 	and b, and edge (a, b). New converted aggregate graph only has one edge, whose original 
	 * 	end-points is a and c, respectively.
	 */
	public static MGraph doConvertAg(Graph ag) {
		nNode = ag.getnNode();
		nEdge = ag.getnEdge();
		cid = new int[nNode];		// component id of nodes
		final int[] nodeInd = ag.getnodeInd();
		final int[] edgeT = ag.getedgeT();
		final double[] w = ag.getw();
		
		int nComp = components(nNode, nodeInd, edgeT, w);
		
		p = new double[nComp];
		edgeCid = new int[nEdge];
		convAgW = new ArrayList<MEdge>();
		mergeNodes(nNode, nodeInd, edgeT, w);
		
		MGraph convAG = new MGraph(nComp);
		convAG.setNode(p, nNode, cid, edgeCid, nodeInd, edgeT);
		convAG.setEdge(convAgW);
		
		clear();
		return convAG;
	}
	
	// generate components using only positive edges;
	// component of each node is indicated by cid[];
	private static int components(final int n, final int nodeInd[], final int edgeT[], 
			final double w[]) {
		for (int i = 0; i < n; i++) {	// initially, each node belongs to a component with only one node
			cid[i] = i;
		}
		int s, t; 	// edge (s, t)
		int p, q; 	// root cid
		for (int i = 0; i < n; i++) {
			s = i;
			p = findRoot(s);
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				t = edgeT[j];
				if (w[j] < 0 || s > t) continue;
				q = findRoot(t);
				cid[q] = p;
			}
		}
		for (int i = 0; i < n; i++) {
			findRoot(i);
		}
		int[] cidMap = new int[n];
		Arrays.fill(cidMap, -1);
		int ind = 0;
		for (int i = 0; i < n; i++) {
			if (cidMap[cid[i]] == -1) cidMap[cid[i]] = ind++;
		}
		for (int i = 0; i < n; i++) {
			cid[i] = cidMap[cid[i]];
		}
		return ind;
	}
	
	/**
	 *  1) merge each connected component into a single node, selecting nodes and edges for 
	 *  each merged nodes (connected components). 
	 *  2) select edges between each merged nodes.
	 */
	
	private static void mergeNodes(final int n, final int nodeInd[], final int edgeT[], 
			final double w[]) {
		Arrays.fill(p, 0);
		Arrays.fill(edgeCid, -1);
		int[] edgeS = new int[nEdge]; 
		TreeMap<Edge, Integer> exist = new TreeMap<Edge, Integer>();
		for (int i = 0; i < n; i++) {
			int s = i;
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				edgeS[j] = s;
				int t = edgeT[j];
				if (s > t) continue;
				if (cid[s] == cid[t]) { 	// s and t belong to same component
					if (!(w[j] < 0)) {			// edge (s, t) is selected by component
						p[cid[s]] += w[j];
						edgeCid[j] = cid[s];
					}
				}
				else {	// s and t belong to different component, here w[j] must < 0
					Edge e = (cid[s] < cid[t]) ? 
							new Edge(cid[s], cid[t]) : new Edge(cid[t], cid[s]);
					Integer eid = exist.get(e);
					if (eid == null || -w[j] < -w[eid]) {
						exist.put(e, j);
					}
				}
			}
		}
		convAgW.ensureCapacity(exist.size() * 2);
		for (Entry<Edge, Integer> entry: exist.entrySet()) {
			int s = entry.getKey().getS();
			int t = entry.getKey().getT();
			int j = entry.getValue();
			MEdge e1 = new MEdge(s, t, edgeS[j], edgeT[j], -w[j]);
			MEdge e2 = new MEdge(t, s, edgeS[j], edgeT[j], -w[j]);
			convAgW.add(e1);
			convAgW.add(e2);
		}		
		Collections.sort(convAgW);
	}
	
	private static int findRoot(int k) {
		if (k == cid[k]) return k;
		cid[k] = findRoot(cid[k]);
		return cid[k];
	}
	
	private static void clear() {
		cid = null;
		p = null;
		edgeCid = null;
		convAgW = null;
	}
}
