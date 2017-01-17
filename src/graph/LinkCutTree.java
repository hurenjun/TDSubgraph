package graph;

import java.util.Arrays;
import java.util.Random;

public class LinkCutTree {
	private int[] left;
	private int[] right;
	private int[] parent;
	private int[] pathParent;
	private double[] costToParent;
	private int[] realParent;
	private int maxEdgeChild;
	private double maxEdgeWeight;
	private int pathNodeCnt;
	private int[] pathNode;
	private double[] pathNodeWeight;
	/*
	 * member variable description:
	 * 1) left, right and parent: left child, right child, and parent of node (aux. tree);
	 * 2) pathParent: path parent pointer of node (aux. tree, and only root is valid);
	 * 3) costToParent: cost(weight) of edge between node and its parent;
	 * 4) realParent: parent node of each node in the represented tree;
	 */
	
	// build a forest with nNode trees, and initially each tree is a single node
	public LinkCutTree(final int nNode) {
		left = new int[nNode];
		Arrays.fill(left, -1);
		right = new int[nNode];
		Arrays.fill(right, -1);
		parent = new int[nNode];
		Arrays.fill(parent, -1);
		pathParent = new int[nNode];
		Arrays.fill(pathParent, -1);
		costToParent = new double[nNode];
		Arrays.fill(costToParent, 0);
		realParent = new int[nNode];
		Arrays.fill(realParent, -1);
		pathNode = new int[nNode];
		pathNodeWeight = new double[nNode];
	}
	
	private boolean isRoot(final int v) {
		return parent[v] == -1;
		//return parent[v] == -1 || (left[parent[v]] != v && right[parent[v]] != v);
	}

	private void connect(final int ch, final int p, Boolean isLeftChild) {
		if (ch != -1)
			parent[ch] = p;
		if (isLeftChild != null) {
			if (isLeftChild)
				left[p] = ch;
			else
				right[p] = ch;
		}
	}

	private void rotate(final int x) {
		int p = parent[x];
		int g = (p == -1) ? -1 : parent[p];
		boolean isRootP = isRoot(p);
		boolean leftChildX = (x == left[p]);
		
		// after rotating with p, now x becomes root (of aux. tree)
		if (g == -1) {
			pathParent[x] = pathParent[p];
			pathParent[p] = -1;
		}

		// create 3 edges: (x.r(l),p), (p,x), (x,g)
		connect(leftChildX ? right[x] : left[x], p, leftChildX);
		connect(p, x, !leftChildX);
		connect(x, g, !isRootP ? p == left[g] : null);
	}

	private void splay(final int x) {
		while (!isRoot(x)) {
			int p = parent[x];
			int g = (p == -1) ? -1 : parent[p];
			if (!isRoot(p)) {
				// rotate p: left-left or right-right, i.e. zig-zig
				// rotate x: left-right or right-left, i.e. zig-zag
				rotate((x == left[p]) == (p == left[g]) ? p : x);
			}
			rotate(x);
		}
	}

	// access node v will form a preferred path from the root down to v; 
	// the new formed preferred path may change other preferred paths;
	private int access(int v) {
		int backupV = v;
		splay(v);
		if (right[v] != -1) {
			pathParent[right[v]] = v;
			parent[right[v]] = -1;
			right[v] = -1;
		}
		
		while (pathParent[v] != -1) {
			int w = pathParent[v];
			splay(w);
			if (right[w] != -1) {
				pathParent[right[w]] = w;
				parent[right[w]] = -1;
			}
			right[w] = v;
			parent[v] = w;
			pathParent[v] = -1;
			v = w;
		}
		
		splay(backupV);
		return v;
	}
	
	public int findRoot(int v) {
		access(v);
		while (left[v] != -1)
			v = left[v];
		access(v);
		return v;
	}
	
	// v is guarantee to be the root of a tree, and w can be arbitrary node
	// the link operator will let v be a child of w
	public void link(int v, int w, double cost) {
		access(v);
		access(w);
		left[v] = w;
		parent[w] = v;
		pathParent[v] = pathParent[w];
		pathParent[w] = -1;
		costToParent[v] = cost;
		realParent[v] = w;
	}

	public void cut(int v) {
		access(v);
		if (left[v] != -1) {
			parent[left[v]] = -1;
			pathParent[left[v]] = pathParent[v];
			left[v] = -1;
			pathParent[v] = -1;
			costToParent[v] = 0;
			realParent[v] = -1;
		}
	}
	
	// let v be the root of the represented tree containing v;
	// method: access(v) first and change the direction of path from root down to v
	public void evert(int v) {
		if (realParent[v] == -1) return;
		access(v);
		pathNodeCnt = 0;
		extractPath(v);
		for (int i = 1; i < pathNodeCnt; i++) { // i = 0 if the root
			cut(pathNode[i]);
			link(pathNode[i - 1], pathNode[i], pathNodeWeight[i]);
		}
	}

	// return the maximum edge of path from v to w
	public int maximumEdgeVW(int v, int w) {
		//int rootV = findRoot(v);
		//int rootW = findRoot(w);
		maxEdgeChild = -1;
		maxEdgeWeight = -1;
		//if (v != w) {
			access(v);
			int lca = access(w);	// least common ancestor
			access(lca);
			if (v != lca) {
				splay(v); maxEdgeOfSubtree(v);
			}
			if (w != lca) {
				splay(w); maxEdgeOfSubtree(w);
			}
		//}
		return maxEdgeChild;
	}
	
	private void maxEdgeOfSubtree(int v) {
		if (costToParent[v] > maxEdgeWeight) {
			maxEdgeWeight = costToParent[v];
			maxEdgeChild = v;
		}
		if (left[v] != -1) maxEdgeOfSubtree(left[v]);
		if (right[v] != -1) maxEdgeOfSubtree(right[v]);
	}
	
	private void extractPath(int v) {
		if (left[v] != -1) extractPath(left[v]);
		pathNode[pathNodeCnt] = v;
		pathNodeWeight[pathNodeCnt++] = costToParent[v];
		if (right[v] != -1) extractPath(right[v]);
	}
	
	public double getMaxEdgeWeight() {
		return maxEdgeWeight;
	}
	
	public int getRealParent(final int v) {
		return realParent[v];
	}
	
	public void clear() {
		left = right = parent = pathParent = null;
		costToParent = null;
		realParent = null;
		pathNode = null;
		pathNodeWeight = null;
	}
	
	public static void main(String[] args) {
		Random rnd = new Random(System.currentTimeMillis());
//		int nNode = 15;
//		LinkCutTree lct = new LinkCutTree(nNode);
//		lct.link(14, 13, 2);
//		lct.link(13, 11, 3);
//		lct.link(10, 8, 1);  lct.link(11, 8, 4); lct.link(12, 8, 1);
//		lct.link(8, 7, 2); lct.link(9, 7, 2);
//		lct.link(7, 6, 5);
//		lct.link(3, 1, 2); lct.link(4, 1, 3); lct.link(5, 1, 4);
//		lct.link(6, 2, 3);
//		lct.link(1, 0, 6); lct.link(2, 0, 1);
//		
//		for (int i = 0; i < nNode; i++) {
//			System.out.print("root[" + i + "]=" + lct.findRoot(i) + "\t");
//		}
//		System.out.println();
//		
//		lct.cut(8);
//		lct.link(8, 0, 2);
//		for (int i = 0; i < 1; i++) {
//			int v = rnd.nextInt(nNode);
//			int w = rnd.nextInt(nNode);
//			//int v = 13;
//			//int w = 9;
//			System.out.print("path " + v +"\t" + w +": " );
//			int child = lct.maximumEdgeVW(v, w);
//			if (child == -1) System.out.println("v and w are same");
//			else System.out.println("(" + child + ", " + lct.realParent[child] + ")\t" 
//					+ lct.costToParent[child]);
//		}
		
		// stress testing
		long st = System.currentTimeMillis();
		int nNode = 100000;
		LinkCutTree lct = new LinkCutTree(nNode);
		for (int i = 0; i < 100000; i++) {
			int cmd = rnd.nextInt(10);
	        if (cmd == 0) {
	        	int trial = 0;
	        	while (true) {
		        	int v = rnd.nextInt(nNode);
		        	if (lct.findRoot(v) != v) {
		        		lct.cut(v); break;
		        	}
		        	trial++;
		        	if (trial == 1000) break;
	        	}
	        }
	        else if (cmd == 1 || cmd == 2 || cmd == 3) {
	        	int trial = 0;
				while (true) {
					int v = rnd.nextInt(nNode);
					int w = rnd.nextInt(nNode);
					if (lct.findRoot(v) != lct.findRoot(w)) {
						int root = lct.findRoot(v);
						lct.link(root, w, (double) (rnd.nextFloat() + 0.5));
						break;
					}
					trial++;
		        	if (trial == 1000) break;
				}
	        }
	        else {
				int trial = 0;
				while (true) {
					int v = rnd.nextInt(nNode);
					int w = rnd.nextInt(nNode);
					if (lct.findRoot(v) == lct.findRoot(w)) {
						lct.maximumEdgeVW(v, w);
						break;
					}
					trial++;
					if (trial == 1000) break;
				}
	        }
		}
		System.out.println(System.currentTimeMillis() - st + "ms");
	}
}