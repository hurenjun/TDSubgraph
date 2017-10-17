package graph;

import java.util.Arrays;
import java.util.TreeSet;

public class DenseSubGraph {
	MGraph graph;
	double rootedNodeW[];
	TreeSet<Edge> selectEdge;
	int nNode, rootId;
	boolean used[];
	int height[];
	TreeSet<Edge> dsgEdge;	// original s and t
	/*
	 * member variable description:
	 * 1) graph: the underlying minimal spanning tree;
	 * 2) rootedNodeW[i]: maximal node weight of subtrees using i as root, 
	 *    cross ref. 'Strong Pruning';
	 * 3) selectEdge: select edges in MST (not edges of original graph)
	 * 4) nNode & rootID: number of nodes & root id of the final selected sub-tree
	 * 5) used[i]: whether node i is used in the final selected sub-tree
	 * 6) height[i]: height if node i in MST (starting from zero, and the 
	 *    root is randomly selected);
	 * 7) DSGEdge: edges of dense subgraph (edges of original graph that are selected
	 *    as dense sub-graph)
	 */
	
	public DenseSubGraph(MGraph _graph) {
		graph = _graph;
		nNode = graph.getnNode();
		rootedNodeW = new double[nNode];
		height = new int[nNode];
		selectEdge = new TreeSet<Edge>();
		rootId = -1;
	}

	public double[] getRootedNodeW() {
		return rootedNodeW;
	}
	
	public int[] getHeight() {
		return height;
	}

	public TreeSet<Edge> getSelectEdge() {
		return selectEdge;
	}

	public int getnNode() {
		return nNode;
	}
	
	// select dense subgraph
	public void selectDSG() {
		rootId = 0;
		for (int i = 1; i < nNode; i++) {
			if (rootedNodeW[i] > rootedNodeW[rootId]) rootId = i;
		}
		
		final int[] nodeInd = graph.getnodeInd();
		final int[] edgeT = graph.getedgeT();
		final int[] oriS = graph.getOriS();
		final int[] oriT = graph.getOriT();
		
		dsgEdge = new TreeSet<Edge>();
		used = new boolean[nNode];
		Arrays.fill(used, false);
		int queue[] = new int[nNode];
		int l = 0, r = 1;
		queue[l] = rootId;
		while (l < r) {
			int s = queue[l];
			used[s] = true;
			//final TreeSet<Edge> edgeInNode = graph.getEdgeInNodeByID(s);
			//Iterator<Edge> iter = edgeInNode.iterator();
			//while (iter.hasNext()) {
			//	Edge e = (Edge)iter.next();
			//DSGEdge.add(new Edge(e));
			//}
			dsgEdge.addAll(graph.getEdgeInNodeByID(s));
			for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
				int t = edgeT[j];
				if (height[t] == height[s] + 1 && selectEdge.contains(new Edge(s, t))) {
					queue[r++] = t;
					dsgEdge.add(new Edge(oriS[j], oriT[j]));
				}
			}
			l++;
		}
		
		//System.out.print(r + "/" + nNode + "\t");
	}
	
	public double getDSGWeight() {
		return rootedNodeW[rootId];
	}
	
	// get sum of weights of nodes which are not selected
	public double compputeRemainNodeW() {
		final double[] nodeW = graph.getNodeW();
		double result = 0;
		for (int i = 0; i < nNode; i++) {
			if (!used[i]) result += nodeW[i];
		}
		return result;
	}
	
	public double computeUsedNodeW() {
		final double[] nodeW = graph.getNodeW();
		double result = 0;
		for (int i = 0; i < nNode; i++) {
			if (used[i]) result += nodeW[i];
		}
		return result;
	}
	
	public double computeAllNodeW() {
		final double[] nodeW = graph.getNodeW();
		double result = 0;
		for (int i = 0; i < nNode; i++) {
			result += nodeW[i];
		}
		return result;
	}
	
	public TreeSet<Edge> listDSGEdges() {		
		return dsgEdge;		
		//return null;
	}
	
	public boolean[] getUsed() {
		return used;
	}
}
