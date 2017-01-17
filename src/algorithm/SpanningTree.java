package algorithm;

import java.util.TreeMap;
import java.util.TreeSet;

import graph.Edge;
import graph.MEdge;
import graph.MGraph;

public class SpanningTree {
	private static int cid[];
	//private static int ss, tt, eeid;
	//private static float ww;
	private static int s[], t[], eid[];
	private static double w[];
	
	public static void doMinimumST(final int nNode, final int nEdge, int _s[], int _t[], 
			double _w[], int _eid[]) {
		s = _s; t = _t;
		w = _w;
		eid = _eid;
		cid = new int[nNode];
		for (int i = 0; i < nNode; i++) {cid[i] = i;}
		quickSort(0, nEdge - 1);
		
		for (int i = 0; i < nEdge; i++) {
			int p = findRoot(s[i]);
			int q = findRoot(t[i]);
			if (p == q) eid[i] = -1;
			else cid[q] = p;
		}
		
		clear(false);
	}
	
	// compute minimal spanning tree on g
	// directly modify the edge of g
	public static MGraph doMinimumST(MGraph g) {
		int nEdge = g.getnEdge();
		s = new int[nEdge / 2];
		t = new int[nEdge / 2];
		w = new double[nEdge / 2];
		eid = new int[nEdge / 2];
		g.copyEdge(s, t, w, eid);
		quickSort(0, nEdge / 2 - 1);
		
		int nNode = g.getnNode();
		cid = new int[nNode];
		for (int i = 0; i < nNode; i++) {cid[i] = i;}
		for (int i = 0; i < nEdge / 2; i++) {
			int p = findRoot(s[i]);
			int q = findRoot(t[i]);
			if (p == q) eid[i] = -1;
			else cid[q] = p;
		}
		
		TreeMap<MEdge, Double> selectEdge = new TreeMap<MEdge, Double>();
		final int[] oriS = g.getOriS();
		final int[] oriT = g.getOriT();
		for (int i = 0; i < nEdge / 2; i++) {
			if (eid[i] == -1) continue;
			MEdge e1 = new MEdge(s[i], t[i], oriS[eid[i]], oriT[eid[i]]);
			MEdge e2 = new MEdge(t[i], s[i], oriS[eid[i]], oriT[eid[i]]);
			selectEdge.put(e1, w[i]);
			selectEdge.put(e2, w[i]);
		}
		
		clear(true);
		
		MGraph mst = new MGraph(nNode);
		mst.setNode(g);
		mst.setEdge(selectEdge);
		return mst;
	}
	
	// do minimum spanning tree based on convAG, and only use edges whose end-points u and v 
	// satisfy that used[finalNodeID[u]] = true and used[finalNodeID[v]] = true;
	public static double doMinimumST(final MGraph convAG, final boolean used[], 
			final int finalNodeID[], TreeSet<Edge> dsgEdges) {
		int nEdge = convAG.getnEdge();
		s = new int[nEdge / 2];
		t = new int[nEdge / 2];
		w = new double[nEdge / 2];
		eid = new int[nEdge / 2];
		convAG.copyEdge(s, t, w, eid);
		quickSort(0, nEdge / 2 - 1);
		
		final int[] oriS = convAG.getOriS();
		final int[] oriT = convAG.getOriT();
		final double[] nodeW = convAG.getNodeW();
		
		double density = 0;
		int nNode = convAG.getnNode();
		cid = new int[nNode];
		for (int i = 0; i < nNode; i++) {
			cid[i] = i;
			if (used[finalNodeID[i]]) {
				dsgEdges.addAll(convAG.getEdgeInNodeByID(i));
				density += nodeW[i];
			}
		}
		
		for (int i = 0; i < nEdge / 2; i++) {
			if (!used[finalNodeID[s[i]]] || !used[finalNodeID[t[i]]]) continue;
			int p = findRoot(s[i]);
			int q = findRoot(t[i]);
			if (p != q) {
				cid[q] = p;
				dsgEdges.add(new Edge(oriS[eid[i]], oriT[eid[i]]));
				density -= w[i];
			}
		}
		
		clear(true);
		//System.out.print(density + "\t");
		//return dsgEdges;
		return density;
	}
	
	public static double doMinimumST(final MGraph convAG, final boolean used[], TreeSet<Edge> dsgEdges) {
		int nEdge = convAG.getnEdge();
		s = new int[nEdge / 2];
		t = new int[nEdge / 2];
		w = new double[nEdge / 2];
		eid = new int[nEdge / 2];
		convAG.copyEdge(s, t, w, eid);
		quickSort(0, nEdge / 2 - 1);
		
		final int[] oriS = convAG.getOriS();
		final int[] oriT = convAG.getOriT();
		final double[] nodeW = convAG.getNodeW();
		
		double density = 0;
		int nNode = convAG.getnNode();
		cid = new int[nNode];
		for (int i = 0; i < nNode; i++) {
			cid[i] = i;
			if (used[i]) {
				dsgEdges.addAll(convAG.getEdgeInNodeByID(i));
				density += nodeW[i];
			}
		}
		
		for (int i = 0; i < nEdge / 2; i++) {
			if (!used[s[i]] || !used[t[i]]) continue;
			int p = findRoot(s[i]);
			int q = findRoot(t[i]);
			if (p != q) {
				cid[q] = p;
				dsgEdges.add(new Edge(oriS[eid[i]], oriT[eid[i]]));
				density -= w[i];
			}
		}
		
		clear(true);
		return density;
	}
	
	private static void quickSort(final int l, final int r) {
		if (l >= r) return;
		int x = l, y = r;
		int ss = s[l]; int tt = t[l]; double ww = w[l]; int eeid = eid[l];
		while (x < y) {
			while (x < y && w[y] > ww) y--;
			if (x < y) {
				s[x] = s[y]; t[x] = t[y]; w[x] = w[y]; eid[x] = eid[y];
				x++;
			}
			while (x < y && w[x] < ww) x++;
			if (x < y) {
				s[y] = s[x]; t[y] = t[x]; w[y] = w[x]; eid[y] = eid[x];
				y--;
			}
		}
		s[x] = ss; t[x] = tt; w[x] = ww; eid[x] = eeid;
		if (l < x - 1) quickSort(l, x - 1);
		if (x + 1 < r) quickSort(x + 1, r);
	}
	
	private static int findRoot(int k) {
		if (k == cid[k]) return k;
		cid[k] = findRoot(cid[k]);
		return cid[k];
	}
	
	
	private static void clear(boolean deep) {
		cid = null;
		if (deep) {
			s = null; t = null; w = null; eid = null;
		}
	}
}
