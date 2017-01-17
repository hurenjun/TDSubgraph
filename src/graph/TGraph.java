package graph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

// Temporal Graph: the nodes and edges are fixed, and the weights vary with time;
public class TGraph {
	int[] nodeInd;		// index of edges starting at i
	int[] edgeT;		// t of edges, while edge is pair (s, t) and t is the right end node
	int nTimestamp;		// number of timestamps
	int startT, endT;	// start timestamp(include) and end timestamp(include)
	float[][] tgw;		// Temporal Graph Weight.
	float[][] aggreW;	// aggregate edge weight
	int nNode, nEdge;	// number of nodes and edges
	HashMap<Integer, Integer> name2id, id2name;
	TreeSet<Edge> origEdge;
	TreeMap<Edge, Integer> edge2ind;	// map from an edge to its index in edgeT
	/* 
	 * Variables Description: 
	 * 1) dgw[t][edge_ind] is the weight of edge(indexed by edge_ind) in t-th timestamp;
	 * 2) aggreW[t][edge_ind] is sum of weights of edge(indexed by edge_ind) from timestamp 0 to t;
	 * 3) nodeInd and edgeT are used for storing sparse graph. Further information is available 
	 *    at https://en.wikipedia.org/wiki/Sparse_matrix (Yale);
	 */
	
	// src is the source of temporal graph
	public TGraph(final String src, final int _startT, final int _endT) {		
		startT = _startT; endT = _endT;
		nTimestamp = endT - startT + 1;
		this.constructGraph(src);
		this.loadTGraphWeights(src);
		aggregateWeights();
		//edge2ind.clear();
		System.out.println("original graph, node: " + nNode + "\tedge: " + nEdge);
	}
	
	public TGraph (final TGraph ances, final int _startT, final int _endT) {
		startT = _startT; endT = _endT;
		nTimestamp = endT - startT + 1;
		nodeInd = ances.nodeInd;
		edgeT = ances.edgeT;	
		nNode = ances.nNode; nEdge = ances.nEdge;
		name2id = ances.name2id; id2name = ances.id2name;
		origEdge = ances.origEdge;
		edge2ind = ances.edge2ind;
		
		tgw = new float[nTimestamp][nEdge];
		for (int i = 0; i < nTimestamp; i++) {
			for (int j = 0; j < nEdge; j++) {
				tgw[i][j] = ances.tgw[i + startT][j];
			}
		}
		
		aggregateWeights();
	}
	
	/*
	 * construct the graph using the first timestamp
	 * 1) about the graph file: each line of the file describe an edge in some timestamp.
	 *    each line has four integers: u,v,t,w (separated by ',' with no other space)  
	 *    for weight w of edge (u,v) in time t. Node ids u and v are not guaranteed to be 
	 *    continuous, while the timestamps t are continuous, i.e. all edges in time 0 come
	 *    first, and followed by edges in time 1, 2... 
	 */
	private void constructGraph(final String src) {
		try {
			nNode = nEdge = 0;
			name2id = new HashMap<Integer, Integer>();
			id2name = new HashMap<Integer, Integer>();
			int currNodeInd = 0;
			BufferedReader br = new BufferedReader(new FileReader(src));
			String line = null;
			int firstT = -1;
			int u, v, t;
			char sepChar = ',';
			Set<Edge> edgeSet = new TreeSet<Edge>();
			origEdge = new TreeSet<Edge>();
			while((line = br.readLine()) != null) {
				if (line.length() == 0) continue;
				int sep1 = line.indexOf(sepChar);
				int sep2 = line.indexOf(sepChar, sep1 + 1);
				int sep3 = line.indexOf(sepChar, sep2 + 1);
				u = Integer.valueOf(line.substring(0, sep1));
				v = Integer.valueOf(line.substring(sep1 + 1, sep2));
				t = Integer.valueOf(line.substring(sep2 + 1, sep3));
				if (firstT == -1) firstT = t;
				else if (firstT != t) break;
				if (!name2id.containsKey(u)) {
					name2id.put(u, currNodeInd);
					id2name.put(currNodeInd, u);
					currNodeInd++;
				}
				if (!name2id.containsKey(v)) {
					name2id.put(v, currNodeInd);
					id2name.put(currNodeInd, v);
					currNodeInd++;
				}
				Edge e1 = new Edge(name2id.get(u), name2id.get(v));
				//Edge e1 = new Edge(u, v);
				edgeSet.add(e1);
				origEdge.add(e1);
				Edge e2 = new Edge(name2id.get(v), name2id.get(u));
				//Edge e2 = new Edge(v, u);
				edgeSet.add(e2);
			}
			br.close();
			
			nNode = name2id.size();
			nodeInd = new int[nNode + 1];
			Arrays.fill(nodeInd, -1);
			nodeInd[0] = 0;
			nEdge = edgeSet.size();
			edgeT = new int[nEdge];
			edge2ind = new TreeMap<Edge, Integer>();
			int ind = 0;
			for (Iterator<Edge> iter = edgeSet.iterator(); iter.hasNext(); ) {
				Edge e = new Edge((Edge)iter.next());
				edge2ind.put(e, ind);
				nodeInd[e.getS() + 1] = ind + 1;
				edgeT[ind] = e.getT();
				//iter.next();
				++ind;
			}
			//System.out.println(nEdge + " " + ind);
			for (int i = 0; i < nNode + 1; i++) { 
				if (nodeInd[i] == -1) nodeInd[i] = nodeInd[i - 1];
			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void loadTGraphWeights(final String src) {
		tgw = new float[nTimestamp][nEdge];
		try {
			BufferedReader br = new BufferedReader(new FileReader(src));
			String line = null;
			int u, v, t, edgeInd;
			float w;
			int totalW = 0, posW = 0;
			char sepChar = ',';
			while((line = br.readLine()) != null) {
				if (line.length() == 0) continue;
				int sep1 = line.indexOf(sepChar);
				int sep2 = line.indexOf(sepChar, sep1 + 1);
				int sep3 = line.indexOf(sepChar, sep2 + 1);
				u = Integer.valueOf(line.substring(0, sep1));
				v = Integer.valueOf(line.substring(sep1 + 1, sep2));
				t = Integer.valueOf(line.substring(sep2 + 1, sep3));
				w = Float.valueOf(line.substring(sep3 + 1));
				totalW += 1;
				if (w > 0) posW += 1;
				if (t < startT || t > endT) continue;
				Edge e1 = new Edge(name2id.get(u), name2id.get(v));
				//Edge e1 = new Edge(u, v);
				edgeInd = edge2ind.get(e1);
				tgw[t - startT][edgeInd] = w;
				Edge e2 = new Edge(name2id.get(v), name2id.get(u));
				//Edge e2 = new Edge(v, u);
				edgeInd = edge2ind.get(e2);
				tgw[t - startT][edgeInd] = w;
				//if (u == v) System.out.println(u);
			}
			br.close();
			System.out.println("positive ratio: " + posW * 1.0 / totalW);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// compute aggreW
	private void aggregateWeights() {
		aggreW = new float[nTimestamp][nEdge];
		//double sum = 0;
		for (int i = 0; i < nEdge; i++) {
			aggreW[0][i] = tgw[0][i];
			//sum += tgw[0][i];
		}
		//System.out.print(sum / 1e5 + "\t");
		
		//double maxRatio = 0, minRatio = 0, eqRatio = 0;
		for (int t = 1; t < nTimestamp; t++) {
			//int inc = 0, dec = 0, eq = 0;
			//sum = 0;
			for (int i = 0; i < nEdge; i++) {
				aggreW[t][i] = aggreW[t - 1][i] + tgw[t][i];
				//sum += tgw[t][i];
				//if (tgw[t][i] > tgw[t - 1][i]) inc++;
				//else if (tgw[t][i] < tgw[t - 1][i]) dec++;
				//else eq++;
			}
			//double r = (double)(Math.max(inc, dec) + eq) / nEdge;
			//maxRatio += r;
			//r = (double)(Math.min(inc, dec) + eq) / nEdge;
			//minRatio += r;
			//r = (double) eq / nEdge;
			//eqRatio += r;
			
			//r = (double) inc / nEdge;
			//System.out.print(r + "\t");
			//System.out.print(sum / 1e5 + "\t");
			
			//if (t % 20 == 0) System.out.println();
		}
		//System.out.println();
		//if (nTimestamp > 1) {
		//	maxRatio /= (nTimestamp - 1);
		//	minRatio /= (nTimestamp - 1);
		//	eqRatio /= (nTimestamp - 1);
		//}
		//System.out.println("Avg maxRate: " + maxRatio);
		//System.out.println("Avg minRate: " + minRatio);
		//System.out.println("Avg eqRate: " + eqRatio);
	}
	
	// Input s: start timestamp (include)
	// Input t: end timestamp (include)
	// Output ag: aggregate graph of time interval [s, t]
	public Graph aggregateGraph(final int s, final int t) {
		if (s < startT || s > endT || t < startT || t > endT) {
			//System.out.println("Time Interval Error!");
			return null;
		}
		Graph ag = new Graph(nNode);
		ag.setEdge(nEdge, nodeInd, edgeT);
		ag.setW(aggreW, s - startT, t - startT);
		return ag;
	}
	
	public double[] getCDensity(){
		double[] cden = new double[nTimestamp];	// cohesive density
		// cden[i] = sum of edge weights in timestamp i
		for (int i = 0; i < nTimestamp; i++) {
			cden[i] = 0;
			for (int j = 0; j < nEdge; j++) { cden[i] += tgw[i][j]; }
			cden[i] /= 2;
		}
		
		return cden;
	}
	
	public double getPosCDen(int s, int t) {
		double pcden = 0;
		if (s > t || s >= nTimestamp || t < 0) return 0;
		if (s < 0) s = 0;
		if (t >= nTimestamp) t = nTimestamp - 1;
		for (int i = 0; i < nEdge; i++) {
			double ew = (s == 0) ? aggreW[t][i] : aggreW[t][i] - aggreW[s - 1][i];
			if (ew >= 0) pcden += ew;
		}
		return pcden / 2;
	}
	
	public int getNTimestamp() {
		return nTimestamp;
	}
	
	public int getNNode() {
		return nNode;
	}
	
	public void printSelectEdge(TreeSet<Edge> selectEdge, int s, int t) {
		double den = 0, ew;
		Iterator<Edge> iter = selectEdge.iterator();
		//System.out.println("Selected Edges");
		while (iter.hasNext()) {
			Edge e = (Edge)iter.next();
			int edgeInd = edge2ind.get(e);
			if (s != 0) ew = (aggreW[t][edgeInd] - aggreW[s - 1][edgeInd]);
			else ew = (aggreW[t][edgeInd]);
			//System.out.println(id2name.get(e.getS()) + "\t" + id2name.get(e.getT()) 
			//		 +  "\t" + ew);
			den += ew;
		}
		//System.out.println("End of Selected Edges");
		System.out.println(den);
	}
	
	public void testGraph() {
		System.out.println("node: " + nNode + "\tedge: " + nEdge + "\tnTimestamp: " + nTimestamp);
		int[] testNode = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
		for (int i = 0; i < testNode.length; i++) {
			if (!name2id.containsKey(testNode[i])) continue;
			int s = name2id.get(testNode[i]);
			for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
				System.out.print(id2name.get(s) + " " + id2name.get(edgeT[j]) + ": ");
				for (int k = 0; k < nTimestamp; k++) {
					System.out.print(tgw[k][j] + " ");
				}
				System.out.println();
			}
		}
	}
	
	public void printID2Name() {
		System.out.println(id2name);
	}
	
	public HashMap<Integer, Integer> getId2name() {
		return id2name;
	}
	
	public double getECRate() {
		long ecEdge = 0;
		for (int t = 1; t < nTimestamp; t++) {
			long GE = 0, LE = 0;
			for (int i = 0; i < nEdge; i++) {
				if (tgw[t][i] >= tgw[t - 1][i]) GE++;
				if (tgw[t][i] <= tgw[t - 1][i]) LE++;
			}
			ecEdge += Math.max(GE, LE);
		}
		double ecr = ecEdge;
		if (nEdge > 0) ecr /= nEdge;
		if (nTimestamp > 1) ecr /= (nTimestamp - 1);
		return ecr;
	}
}
