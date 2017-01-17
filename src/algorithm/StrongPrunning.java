package algorithm;

import java.util.Arrays;
import java.util.TreeSet;

import graph.DenseSubGraph;
import graph.Edge;
import graph.MGraph;

public class StrongPrunning {
	static int nNode, nEdge;
	static double rootedNodeW[];
	static boolean validNode[];
	static int queue[];
	static int height[];
	static TreeSet<Edge> selectEdge;
	/* 
	 * member variables descriptions:
	 * 1) rootedNodeW[i]: maximum node weight of subtree whose root is node i
	 * 2) validNode[i]: whether node i is still valid (not visited)
	 * 3) queue: the visiting sequence, nodes in this queue is sorted according to 
	 *    their height in the tree for specific root;
	 */
	
	
	/*
	 *  1) Select a subtree of given MST that maximizes the NW-Score of subtree;
	 *  2) NW-Score of a subtree is the sum of node weights minus the sum of edge weight;
	 *  3) It has been shown that the NW-Score maximization problem for general graphs is
	 *     NP-hard, without any constant-factor approximation algorithm, but there exists 
	 *     exact algorithm for NW-Score maximization problem for trees, which takes linear 
	 *     time of the number of nodes in trees;
	 *  4) The algorithm consists of tree steps: 
	 *     a) randomly select a node, and construct tree T using the selected node as root.
	 *        by now each nodes in tree has one parent (exclude the root node) and 0 or 1 or
	 *        several children, these are reflected in height of nodes;
	 *     b) for each node v, compute the maximum NW-Score among all subtrees which are rooted 
	 *        at v, (this step can be done in a bottom-up manner, which is equivalent to Strong 
	 *        Pruning);
	 *     c) find the node with maximal NW-Score, and the subtree rooted as this node is the 
	 *        optimal subtree;
	 */ 
	public static DenseSubGraph doStrongPrunning(MGraph g) {
		DenseSubGraph dsg = new DenseSubGraph(g);
		nNode = g.getnNode();
		nEdge = g.getnEdge();
		final int[] nodeInd = g.getnodeInd();
		final int[] edgeT = g.getedgeT();
		final double w[] = g.getw();
		rootedNodeW = dsg.getRootedNodeW();
		g.copyNodeW(rootedNodeW);	// initially rootedNodeW is the node weight of the node
		height = dsg.getHeight();
		selectEdge = dsg.getSelectEdge();
		
		int root = 0;
		bfsTravel(root, nodeInd, edgeT);
		computeRootedNodeW(nodeInd, edgeT, w);
		clear();
		dsg.selectDSG();
		return dsg;
	}
	
	private static void bfsTravel(final int root, final int nodeInd[], final int edgeT[]) {
		validNode = new boolean[nNode];
		Arrays.fill(validNode, true);
		queue = new int[nNode];
		
		int l = 0, r = 1;
		queue[l] = root;
		height[root] = 0;
		validNode[root] = false;
		while (l < r) {
			int s = queue[l];
			for (int j = nodeInd[s]; j < nodeInd[s + 1]; j++) {
				int t = edgeT[j];
				if (!validNode[t]) continue;
				queue[r++] = t;
				validNode[t] = false;
				height[t] = height[s] + 1;
			}
			l++;
		}
	}
	
	private static void computeRootedNodeW(final int nodeInd[], final int edgeT[], final double w[]) { 		
		for (int i = nNode - 1; i >= 0; i--) {
			int t = queue[i];
			for (int j = nodeInd[t]; j < nodeInd[t + 1]; j++) {
				int s = edgeT[j];
				if (height[t] != height[s] + 1) continue;
				// directed edge (s, t)  is in the tree
				if (rootedNodeW[t] >= w[j]) {
					rootedNodeW[s] += (rootedNodeW[t] - w[j]);
					selectEdge.add(new Edge(s, t));
					selectEdge.add(new Edge(t, s));
				}
			}
		}
	}
	
	private static void clear() {
		validNode = null;
		queue = null;
	}
}
