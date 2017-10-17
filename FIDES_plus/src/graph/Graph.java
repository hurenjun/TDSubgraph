package graph;

import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

public class Graph {
	int[] indFrom;		// index of edges starting at i
	int[] to;		// t of edges, while edge is pair (s, t) and t is the right end node
	double[] w;		// Temporal Graph Weight.
	int nNode, nEdge;	// number of nodes and edges
	
	public Graph(final int _nNode) {
		nNode = _nNode; 
	}
	
	public void setEdge(final int _nEdge, final int[] _nodeInd, final int[] _edgeT) {
		nEdge = _nEdge;
		indFrom = null;
		indFrom = new int[nNode + 1];
		to = null;
		to = new int[nEdge];
		for (int i = 0; i < nNode + 1; i++) { 
			indFrom[i] = _nodeInd[i]; 
		}
		for (int i = 0; i < nEdge; i++) { 
			to[i] = _edgeT[i]; 
		}
	}
	
	// using aggregate weight to set edge weight
	public void setW(final float[][] aggreW, final int s, final int t) {
		w = new double[nEdge];
		//boolean zeroW = false;
		if (s == 0)
			for (int i = 0; i < nEdge; i++) {
				w[i] = aggreW[t][i];
			}
		else
			for (int i = 0; i < nEdge; i++) {
				w[i] = aggreW[t][i] - aggreW[s - 1][i];
				//if (w[i] == 0) zeroW = true;
			}
		//if (zeroW) System.out.println("zeor edge weight");
	}
	
	public double computeSumOfEdgeW(final TreeSet<Edge> selectEdge) {
		double sum = 0;
		//int hit = 0;
		for (int i = 0; i < nNode; i++) {
			for (int j = indFrom[i]; j < indFrom[i + 1]; j++) {
				if (selectEdge.contains(new Edge(i, to[j]))) {
					sum += w[j];
					//hit++;
				}
			}
		}
		//System.out.println("size: " + selectEdge.size() + "\thit: " + hit);
		
//		boolean duplicated = false;
//		Iterator<Edge> iter = selectEdge.iterator();
//		while (iter.hasNext()) {
//			Edge e = (Edge)iter.next();
//			if (selectEdge.contains(new Edge(e.getT(), e.getS()))) duplicated = true;
//		}
//		if (duplicated) System.out.print("duplicated\t");
		
		return sum;
	}
	
	public double[] computeSumOfEdgeW(final TreeSet<Edge>[] edgeInNode) {
		double[] soew = new double[edgeInNode.length];
		Arrays.fill(soew, 0);
		TreeMap<Edge, Integer> edge2NID = new TreeMap<Edge, Integer>();
		for (int i = 0; i < edgeInNode.length; i++) {
			Iterator<Edge> iter = edgeInNode[i].iterator();
			while (iter.hasNext()) {
				edge2NID.put(new Edge((Edge)iter.next()), i);
			}
		}
		for (int i = 0; i < nNode; i++) {
			for (int j = indFrom[i]; j < indFrom[i + 1]; j++) {
				if (edge2NID.containsKey(new Edge(i, to[j]))) {
					soew[edge2NID.get(new Edge(i, to[j]))] += w[j];
				}
			}
		}
		return soew;
	}
	
	public int getnNode() {
		return nNode;
	}

	public int getnEdge() {
		return nEdge;
	}
	
	public final int[] getnodeInd() {
		return indFrom;
	}
	
	public final int[] getedgeT() {
		return to;
	}
	
	public final double[] getw() {
		return w;
	}
	
	public void copyW(double[] _w) {
		for (int i = 0; i < nEdge; i++) {
			_w[i] = w[i];
		}
	}
	
	public void copyEdge(int[] s, int[] t, double[] _w, int[] eid) {
		int ind = 0;
		for (int i = 0; i < nNode; i++) {
			for (int j = indFrom[i]; j < indFrom[i + 1]; j++) {
				if (i > to[j]) continue;
				s[ind] = i; t[ind] = to[j]; 
				_w[ind] = w[j]; eid[ind] = j;
				ind++;
			}
		}
	}
	
	public double getEdgeWeight(final int s, final int t) {
		for (int j = indFrom[s]; j < indFrom[s + 1]; j++) {
			if (to[j] == t) return w[j];
		}
		return 0;
	}
	
	public int getPosEdgeCount() {
		int pec = 0;
		for (int i = 0; i < nNode; i++) {
			for (int j = indFrom[i]; j < indFrom[i + 1]; j++) {
				if (w[j] > 0) pec++;
			}
		}
		return pec / 2;
	}
	
	public double getPosCDen() {
		double pcden = 0;
//		int[] hit = new int[nNode];
//		Arrays.fill(hit, 0);
//		for (int i = 0; i < nEdge; i++) {
//			double ew = w[i];
//			if (ew > 0) { 
//				hit[to[i]]++;
//			}
//		}
		for (int i = 0; i < nEdge; i++) {
			double ew = w[i];
			if (ew > 0) { 
				//pcden += ew * hit[to[i]];
				pcden += ew;
			}
		}
		return pcden;
	}

	public void testGraph() {
		System.out.println("node: " + nNode + "\tedge: " + nEdge);
//		for (int i = 0; i < nNode; i++) {
//			System.out.print(i + ": ");
//			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
//				System.out.print(edgeT[j] + " " + w[j] + "\t");
//			}
//			System.out.println();
//		}
	}
}
