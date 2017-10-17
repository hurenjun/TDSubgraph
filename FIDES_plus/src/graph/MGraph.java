package graph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;


// 1) MGraph extends class Graph, where Graph is the class for general graph.
// 2) MGraph has same feature as Graph. 
// 3) The difference between MGraph and Graph is that nodes of MGraph have node weights and  
// are subgraphs of general graph, for example:
//		Graph g = {{a, b, c, d, e}, {(a, b), (b, c), (a, c), (c, d), (d, e)}}
// 		MGraph mg = {{A, B}, {(A, B)}}, where A = {{a, b, c}, {(a, b), (b, c)}}, 
// 		B = {{d, e}, {(d, e)}}, and nw(A) = 3 / nw(B) = 2;
// 4) MGraph will additionally maintain edges of Graph (e.g., (a,b)) in nodes 
// of MGraph (e.g., A);
public class MGraph extends Graph{
	public TreeSet<Edge>[] edgeInNode;		// (General) edges in (MGraph) nodes
	//TreeSet<Integer>[] nodeInNode;	// (General) nodes in (MGraph) nodes
	double[] nodeW;
	int[] oriS;
	int[] oriT;
	/*
	 * member variables description:
	 * 1) member variables corresponding to node: nNode, nodeInNode, edgeInNode;
	 * 2) member variables corresponding to edge: nEdge, nodeInd, edgeT, w, oriS, oriT
	 */
	
	public MGraph(final int _nNode) {
		super(_nNode);
	}
	
	public void setNode(final double p[], final int n, final int cid[], final int edgeCid[], 
			final int nodeInd[], final int edgeT[]) {
		nodeW = new double[nNode];
		//nodeInNode = new TreeSet[nNode];
		edgeInNode = new TreeSet[nNode];
		for (int i = 0; i < nNode; i++) {
			//nodeInNode[i] = new TreeSet<Integer>();
			edgeInNode[i] = new TreeSet<Edge>();
		}
		for (int i = 0; i < nNode; i++) {
			nodeW[i] = p[i];
		}
		//for (int i = 0; i < n; i++) {
		//	nodeInNode[cid[i]].add(i);
		//}
		
		for (int i = 0; i < n; i++) {
			for (int j = nodeInd[i]; j < nodeInd[i + 1]; j++) {
				if (edgeCid[j] != -1) {
					Edge e = new Edge(i, edgeT[j]);
					edgeInNode[edgeCid[j]].add(e);
				}
			}
		}
	}
	
	public void setNode(final MGraph oriG) {
		nodeW = new double[nNode];
		//nodeInNode = new TreeSet[nNode];
		edgeInNode = new TreeSet[nNode];
		for (int i = 0; i < nNode; i++) {
			nodeW[i] = oriG.nodeW[i];
			//nodeInNode[i] = new TreeSet<Integer>();
			//nodeInNode[i].addAll(oriG.nodeInNode[i]);
			edgeInNode[i] = new TreeSet<Edge>();
			edgeInNode[i].addAll(oriG.edgeInNode[i]);
		}
	}
	
	public void setEdge(final ArrayList<MEdge> edge) {
		indFrom = null; indFrom = new int[nNode + 1];
		Arrays.fill(indFrom, -1);
		indFrom[0] = 0;
		
		nEdge = edge.size();
		to = null; to = new int[nEdge];
		w = null; w = new double[nEdge];
		oriS = null; oriS = new int[nEdge];
		oriT = null; oriT = new int[nEdge];
		
		int ind = 0;
		for (MEdge e: edge) {
			indFrom[e.getS() + 1] = ind + 1;
			to[ind] = e.getT();
			w[ind] = e.getW();
			//if (w[ind] <= 0) System.out.println("Error Edge Weight!");
			oriS[ind] = e.getOriS();
			oriT[ind] = e.getOriT();
			ind++;
		}
		
		for (int i = 0; i < nNode + 1; i++) {
			if (indFrom[i] == -1) indFrom[i] = indFrom[i - 1];
		}
	}
	
	public void copyNodeW(double _nodeW[]) {
		for (int i = 0; i < nNode; i++) {
			_nodeW[i] = nodeW[i];
		}
	}
	
	public final double[] getNodeW() {
		return nodeW;
	}
	
	public final int[] getOriS() {
		return oriS;
	}
	
	public final int[] getOriT() {
		return oriT;
	}
	
	public final TreeSet<Edge> getEdgeInNodeByID(final int nid) {
		return edgeInNode[nid];
	}
	
	// merge nodes into components, and return a new graph where components of 
	// original graph are nodes of new graph;
	public MGraph merge(final int cid[], final int edgeCid[], final double _nodeW[],
			final TreeMap<Integer, MEdge> minCCedge[]) {
		HashMap<Integer, Integer> cidMap = new HashMap<Integer, Integer>();
		int compCnt = 0;
		for (int i = 0; i < nNode; i++) {
			if (!cidMap.containsKey(cid[i])) cidMap.put(cid[i], compCnt++);
		}
//		for (int i = 0; i < nNode; i++) {
//			cid[i] = cidMap.get(cid[i]);
//		}
//		int compNEdge = 0;
//		for (int i = 0; i < nEdge; i++) {
//			if (edgeCid[i] != -1) compNEdge++;
//		}
//		System.out.print("(" + nNode + " " + compCnt + " " + compNEdge + ")\t");
		
		MGraph resultG = new MGraph(compCnt);
		resultG.nodeW = new double[compCnt];
		Arrays.fill(resultG.nodeW, 0);
		//resultG.nodeInNode = new TreeSet[compCnt];
		resultG.edgeInNode = new TreeSet[compCnt];
		for (int i = 0; i < compCnt; i++) {
			//resultG.nodeInNode[i] = new TreeSet<Integer>();
			resultG.edgeInNode[i] = new TreeSet<Edge>();
		}
		
		for (int i = 0; i < nNode; i++) {
			int s = i;
			if (cid[i] == i) resultG.nodeW[cidMap.get(i)] = _nodeW[i];
			//Iterator<Integer> iterNode = nodeInNode[s].iterator();
			//while (iterNode.hasNext()) {
			//	Integer nodeInt = new Integer((Integer)iterNode.next());
			//	resultG.nodeInNode[cidMap.get(cid[s])].add(nodeInt);
			//}
			Iterator<Edge> iterEdge = edgeInNode[s].iterator();
			while (iterEdge.hasNext()) {
				Edge e = new Edge((Edge)iterEdge.next());
				resultG.edgeInNode[cidMap.get(cid[s])].add(e);
			}
			
			for (int j = indFrom[s]; j < indFrom[s + 1]; j++) {
				int t = to[j];
				if (s > t) continue;
				if (edgeCid[j] != -1) {
					resultG.edgeInNode[cidMap.get(cid[edgeCid[j]])]
							.add(new Edge(oriS[j], oriT[j]));
				}
			}
		}
		
		ArrayList<MEdge> resultGEdge = new ArrayList<MEdge>(nEdge);
		for (int i = 0; i < nNode; i++) {
			if (cid[i] != i) continue;
			for (Entry<Integer, MEdge> entry: minCCedge[i].entrySet()) {
				MEdge me = entry.getValue();
				MEdge meAdd = new MEdge(cidMap.get(i), cidMap.get(entry.getKey()), 
						me.getOriS(), me.getOriT(), me.getW());
				resultGEdge.add(meAdd);
			}
		}
		Collections.sort(resultGEdge);
		resultG.setEdge(resultGEdge);
		
		return resultG;
	}
	
	// 1) return a new graph which merges all nodes u where used[u]=true into a single component,
	// while other vertices remains unchanged;
	// 2) if exists multiple edges between nodes and component, select the edge with minimal 
	// weight;
//	public MGraph merge(final boolean used[], final double dsgDens, 
//			final TreeSet<Edge> dsgEdges) 
//	{
//		int nodeIDMap[] = new int[nNode];
//		int ind = 1;
//		for (int i = 0; i < nNode; i++) {
//			if (used[i]) {
//				nodeIDMap[i] = 0;
//			}
//			else nodeIDMap[i] = ind++;
//		}
//		MGraph resultG = new MGraph(ind);
//		
//		// operations correspond to nodes
//		resultG.nodeW = new double[resultG.nNode];
//		resultG.nodeW[0] = dsgDens;
//		//resultG.nodeInNode = new TreeSet[resultG.nNode];
//		resultG.edgeInNode = new TreeSet[resultG.nNode];
//		//resultG.nodeInNode[0] = new TreeSet<Integer>();
//		resultG.edgeInNode[0] = new TreeSet<Edge>(dsgEdges);
//		// the next two lines are for testing
//		//resultG.edgeInNode[0].clear();
//		//System.out.println(dsgEdges.size() + " " + resultG.edgeInNode[0].size());
//		for (int i = 0; i < nNode; i++) {
//			if (used[i]) //resultG.nodeInNode[0].addAll(nodeInNode[i]);
//				continue;
//			else {
//				resultG.nodeW[nodeIDMap[i]] = nodeW[i];
//				//resultG.nodeInNode[nodeIDMap[i]] = new TreeSet<Integer>(nodeInNode[i]);
//				resultG.edgeInNode[nodeIDMap[i]] = new TreeSet<Edge>(edgeInNode[i]);
//			}
//		}
//		
//		// operations correspond to edges
//		TreeMap<MEdge, Double> resultGEdge = new TreeMap<MEdge, Double>();
//		for (int i = 0; i < nNode; i++) {
//			int s = i;
//			for (int j = indFrom[i]; j < indFrom[i + 1]; j++) {
//				int t = to[j];
//				if (s > t) continue; 
//				if (nodeIDMap[s] != nodeIDMap[t]) {
//					MEdge me1 = new MEdge(nodeIDMap[s], nodeIDMap[t], oriS[j], oriT[j]);
//					MEdge me2 = new MEdge(nodeIDMap[t], nodeIDMap[s], oriS[j], oriT[j]);
//					if (!resultGEdge.containsKey(me1)) {
//						resultGEdge.put(me1, w[j]);
//						resultGEdge.put(me2, w[j]);
//					}
//					else if (resultGEdge.get(me1) > w[j]) {
//						resultGEdge.remove(me1);
//						resultGEdge.remove(me2);
//						resultGEdge.put(me1, w[j]);
//						resultGEdge.put(me2, w[j]);
//					}
//				}
//			}
//		}
//		resultG.setEdge(resultGEdge);
//		
//		return resultG;
//	}
	
	public int getValidNode() {
		int vn = 0;
		System.out.print("{");
		for (int i = 0; i < nNode; i++) {
			if (nodeW[i] > 0) {
				vn++;
				System.out.print(nodeW[i] + " ");
			}
		}
		System.out.print("}");
		return vn;
	}
	
	// get sum of edge weights within nodes(components)
	public double[] getNodeSEW(final int cid[], final int edgeCid[]) {
		double[] nodeSEW = new double[nNode];
		Arrays.fill(nodeSEW, 0);
		for (int i = 0; i < nNode; i++) {
			nodeSEW[cid[i]] += nodeW[i];
		}
		for (int i = 0; i < nEdge; i++) {
			if (edgeCid[i] != -1) nodeSEW[cid[edgeCid[i]]] -= w[i];
		}
		return nodeSEW;
	}
	
	// check whether the nodeW of the largest node equals to the sum of weight of edges
	// within the node;
	// here nodes are inherently components;
	public void checkNodeW(Graph ag) {
		System.out.println("node: " + nNode + "\tedge: " + nEdge);
		int rootID = 0;
		for (int i = 1; i < nNode; i++) {
			if (nodeW[i] > nodeW[rootID]) rootID = i;
		}
		
		System.out.println("maximal node weight: " + nodeW[rootID]);
		System.out.println("sum of edge weight: " + ag.computeSumOfEdgeW(edgeInNode[rootID]));
		
//		int mis = 0;
//		for (int i = 0; i < 10000; i++) {
//			if (nodeW[i] != ag.GetSumOfEdgeW(edgeInNode[i])) mis++;
//		}
//		System.out.println("mis:" + mis + "\t");
	}
	
	// print graph information
	public void testGraph(final HashMap<Integer, Integer> id2name) {
		System.out.println("node: " + nNode + "\tedge: " + nEdge);
		for (int i = 0; i < nNode; i++) {
			System.out.print(i + ": ");
			for (int j = indFrom[i]; j < indFrom[i + 1]; j++) {
				System.out.print(to[j] + " " + w[j] + "(" + id2name.get(oriS[j]) + " "
						+ id2name.get(oriT[j]) + ")\t");
			}
			System.out.println();
		}

		for (int i = 0; i < nNode; i++) {
			System.out.print(i + ": ");
			System.out.print("NW{" + nodeW[i] + "}\t");
//			System.out.print("node{");
//			for (Iterator<Integer> iter = nodeInNode[i].iterator(); iter
//					.hasNext();) {
//				Integer node = (Integer) iter.next();
//				System.out.print(node + " ");
//			}
//			System.out.print("}\t");
			System.out.print("edge{");
			for (Iterator<Edge> iter = edgeInNode[i].iterator(); iter.hasNext();) {
				Edge e = (Edge) iter.next();
				System.out.print("(" + id2name.get(e.getS()) + " " + 
						id2name.get(e.getT()) + ") ");
			}
			System.out.println("}");
		}
	}
}
