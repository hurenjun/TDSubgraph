package algorithm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;

import graph.Edge;
import graph.TGraph;

public class IntervalSelector {
	static TGraph tg;
	static int T, nNode;
	static double avgChange, xiSmooth;
	static double[] localAvgChange;
	static double[] cden;
	static TreeMap<Edge, Double> pcden;
	/*
	 * member variable description: 
	 * 1) T: number of timestamps; nNode: number of node;
	 * 2) xi: smooth factor, based on average change of cden, 
	 *    i.e. xi = alpha * \Sigma(cden[i + 1] - cden[i]) / (T - 1), where alpha is fixed to 2;
	 * 3) cden[i]: cohesive density of time i, i.e. cden[i] = sum of edge weights in time i;
	 */
	
	/**
	 * DoIntervalSelect: select k time intervals such that 
	 * 1) half of the k intervals are select by maxTInterval, and another half are 
	 *    selected by minTInterval;
	 * 2) algorithm maxTInterval uses local maximums to select candidate intervals, 
	 *    and minTInterval on the other hand uses local minimums;
	 * 3) both algorithms select the k/2 interval from candidates according the positive
	 *    cohesive density of each interval;
	 */
//	public static Edge[] doIntervalSelect(final TGraph tg, final int k) {
//		T = tg.getNTimestamp();
//		nNode = tg.getNNode();
//		cden = tg.getNormCDensity(); // normalized cohesive density
//		avgChange = 0;	
//		for (int i = 0; i < T - 1; i++) {
//			avgChange += Math.abs(cden[i + 1] - cden[i]);
//		}
//		if (T > 1) avgChange /= (T - 1);
//		
//	}
	
	/**
	 * DoIntervalSelect: select k time intervals such that 
	 * 1) half of the k intervals are select by maxTInterval, and another half are 
	 *    selected by minTInterval;
	 * 2) algorithm maxTInterval uses local maximums to select candidate intervals, 
	 *    and minTInterval on the other hand uses local minimums;
	 * 3) both algorithms select the k/2 interval from candidates according the positive
	 *    cohesive density of each interval;
	 */
	public static Edge[] doIntervalSelect(final TGraph tgraph, int k) {
		tg = tgraph;
		T = tg.getNTimestamp();
		nNode = tg.getNNode();
		cden = tg.getNormCDensity(); // normalized cohesive density
		double xi = 1.0;
		
//		for (int i = 0; i < T; i++) {
//			System.out.print(i + " " + cden[i] + "\t");
//			if (i % 10 == 9) System.out.println();
//		}
//		System.out.println();
		
		avgChange = 0;	
		for (int i = 0; i < T - 1; i++) {
			avgChange += Math.abs(cden[i + 1] - cden[i]);
		}
		if (T > 1) avgChange /= (T - 1);
		xiSmooth = avgChange * xi;
		//System.out.println("xiSmooth:" + xiSmooth);
		if (T >= 5) {
			localAvgChange = new double[T];
			localAvgChange[1] = localAvgChange[0] 
					= (Math.abs(cden[1] - cden[0]) 
					+ Math.abs(cden[2] - cden[1])
					+ Math.abs(cden[3] - cden[2])
					+ Math.abs(cden[4] - cden[3])) / 4;
			for (int i = 2; i < T - 2; i++) {
				localAvgChange[i] = (Math.abs(cden[i] - cden[i - 1]) 
						+ Math.abs(cden[i - 1] - cden[i - 2])
						+ Math.abs(cden[i] - cden[i + 1])
						+ Math.abs(cden[i + 1] - cden[i + 2])) / 4;
				
			}
			localAvgChange[T - 2] = localAvgChange[T - 1] 
					= (Math.abs(cden[T - 1] - cden[T - 2]) 
					+ Math.abs(cden[T -2] - cden[T - 3])
					+ Math.abs(cden[T -3] - cden[T - 4])
					+ Math.abs(cden[T -4] - cden[T - 5])) / 4;
			for (int i = 0; i < T; i++) {
				if (localAvgChange[i] > avgChange) 
					localAvgChange[i] = avgChange;
			}
		}
		
		pcden = new TreeMap<Edge, Double>();
		Edge[] topK = new Edge[k];
		for (int i = 0; i < k; i++) {
			topK[i] = new Edge();
		}
		
		// intervals generating by local maxima
		ArrayList<IntvPCden> tuning = maxTInterval(2 * k);	// top-2k intervals 
		for (IntvPCden item: tuning) {
			//System.out.print(item.s + " " + item.t + "\t-->\t");
			earlyTuning(item);
			//System.out.println(item.s + " " + item.t);
		}
		Collections.sort(tuning);
		int count = 0;
		for (IntvPCden item: tuning) {
			boolean dup = false;
			for (int j = 0; j < count; j++) {
				if (topK[j].getS() == item.s && topK[j].getT() == item.t)
					dup = true;
			}
			if (!dup) {
				topK[count++] = new Edge(item.s, item.t);
				if (count == k / 2) break;
			}
		}
		
		// intervals generating by local minima
		tuning = minTInterval(2 * k);
		for (IntvPCden item: tuning) {
			earlyTuning(item);
		}
		Collections.sort(tuning);
		count = 0;
		for (IntvPCden item: tuning) {
			boolean dup = false;
			for (int j = 0; j < k / 2 + count; j++) {
				if (topK[j].getS() == item.s && topK[j].getT() == item.t)
					dup = true;
			}
			if (!dup) {
				topK[k / 2 + count] = new Edge(item.s, item.t);
				count++;
				if (count == k / 2) break;
			}
		}
		
		clear();
		return topK;
	}
	
	/**
	 * select candidate intervals using local maxima
	 */
	private static ArrayList<IntvPCden> maxTInterval(int capacity) {		
		if (T <= 10) {
			ArrayList<IntvPCden> tuning = new ArrayList<IntvPCden>(T * (T + 1) / 2);
			for (int i = 0; i < T; i++) {
				for (int j = i; j < T; j++) {
					tuning.add(new IntvPCden(i, j, tg.getPosCDen(i, j)));
				}
			}
			return tuning;
		}
		
		int[] l = new int[nNode];	// lower sides of maximums
		int[] u = new int[nNode];	// upper sides of maximums
		int h = 0;
		
		for (int i = 0; i < T; i++) {
			if ((i - 4 < 0 || cden[i] > cden[i - 4])
					&& (i - 3 < 0 || cden[i] > cden[i - 3])
					&& (i - 2 < 0 || cden[i] > cden[i - 2])
					&& (i - 1 < 0 || cden[i] > cden[i - 1])
					&& (i + 1 >= T || cden[i] > cden[i + 1])
					&& (i + 2 >= T || cden[i] > cden[i + 2])
					&& (i + 3 >= T || cden[i] > cden[i + 3])
					&& (i + 4 >= T || cden[i] > cden[i + 4])) {
				smoothLocal(l, u, h, i, true);
				h++;
			}
		}
		
		TreeSet<Edge> intervals = new TreeSet<Edge>();
		for (int i = 0; i < h; i++) {
			intervals.add(new Edge(l[i], u[i]));
		}
		h = 0;
		Edge combIntv = null;
		for (Edge e: intervals) {
			if (combIntv == null) combIntv = new Edge(e);
			else {
				if (e.getS() <= combIntv.getT()) {
					if (e.getT() > combIntv.getT()) combIntv.setT(e.getT());
				}
				else {
					l[h] = combIntv.getS(); u[h] = combIntv.getT(); h++;
					combIntv = new Edge(e);
				}
			}
		}
		l[h] = combIntv.getS(); u[h] = combIntv.getT(); h++;
		//System.out.println("# maxima: " + h);
		
//		int count = 0;
//		for (int i = 0; i < h; i++) {
//			System.out.print("[" + l[i] + "," + u[i] + "]\t");
//			if (count % 10 == 9) System.out.println();
//			count++;
//		}
//		System.out.println("end of maxima");		
		
		TreeSet<IntvPCden> top2k = new TreeSet<IntvPCden>();
		for (int step = 0; step < h; step++) {
			for (int start = 0; start < h - step; start++) {
				int left = l[start], right = u[start + step];
				if (left > right) continue;
				if (top2k.size() < capacity) {
					double pcd = mapSearch(left, right);
					top2k.add(new IntvPCden(left, right, pcd));
				} 
				else {
					double pcd = mapSearch(left, right);
					if (pcd > top2k.last().d) {
						top2k.pollLast();
						top2k.add(new IntvPCden(left, right, pcd));
					}
				}
			}
		}
		
		ArrayList<IntvPCden> tuning = new ArrayList<IntvPCden>(capacity);
		for (IntvPCden item: top2k) {
			tuning.add(new IntvPCden(item));
		}
		return tuning;
	}
	
	/**
	 * select candidate intervals using local minimums
	 */
	private static ArrayList<IntvPCden> minTInterval(int capacity) {
		if (T <= 10) {
			ArrayList<IntvPCden> tuning = new ArrayList<IntvPCden>(T * (T + 1) / 2);
			for (int i = 0; i < T; i++) {
				for (int j = i; j < T; j++) {
					tuning.add(new IntvPCden(i, j, tg.getPosCDen(i, j)));
				}
			}
			return tuning;
		}
		
		int[] l = new int[nNode];	// lower sides of maximums
		int[] u = new int[nNode];	// upper sides of maximums
		int h = 0;
		
		for (int i = 0; i < T; i++) {
			if ((i - 4 < 0 || cden[i] < cden[i - 4])
					&& (i - 3 < 0 || cden[i] < cden[i - 3])
					&& (i - 2 < 0 || cden[i] < cden[i - 2])
					&& (i - 1 < 0 || cden[i] < cden[i - 1])
					&& (i + 1 >= T || cden[i] < cden[i + 1])
					&& (i + 2 >= T || cden[i] < cden[i + 2])
					&& (i + 3 >= T || cden[i] < cden[i + 3])
					&& (i + 4 >= T || cden[i] < cden[i + 4])) {
				smoothLocal(l, u, h, i, true);
				h++;
			}
		}
		
		TreeSet<Edge> intervals = new TreeSet<Edge>();
		for (int i = 0; i < h; i++) {
			intervals.add(new Edge(l[i], u[i]));
		}
		h = 0;
		Edge combIntv = null;
		for (Edge e: intervals) {
			if (combIntv == null) combIntv = new Edge(e);
			else {
				if (e.getS() <= combIntv.getT()) {
					if (e.getT() > combIntv.getT()) combIntv.setT(e.getT());
				}
				else {
					l[h] = combIntv.getS(); u[h] = combIntv.getT(); h++;
					combIntv = new Edge(e);
				}
			}
		}
		l[h] = combIntv.getS(); u[h] = combIntv.getT(); h++;
		//System.out.println("# minima: " + h);
		
//		int count = 0;
//		for (int i = 0; i < h; i++) {
//			System.out.print("[" + l[i] + "," + u[i] + "]\t");
//			if (count % 10 == 9) System.out.println();
//			count++;
//		}
//		System.out.println("end of minima");
		
		for (int i = 0; i < h; i++) {
			for (int j = i + 1; j < h; j++) {
				if (u[i] <= l[j]) {
					intervals.add(new Edge(u[i], l[j]));
				}
			}
		}
		TreeSet<IntvPCden> top2k = new TreeSet<IntvPCden>();
		for (int step = 1; step < h; step++) {
			for (int start = 0; start < h - step; start++) {
				int left = u[start], right = l[start + step];
				if (left > right) continue;
				if (top2k.size() < capacity) {
					double pcd = mapSearch(left, right);
					top2k.add(new IntvPCden(left, right, pcd));
				} 
				else {
					double pcd = mapSearch(left, right);
					if (pcd > top2k.last().d) {
						top2k.pollLast();
						top2k.add(new IntvPCden(left, right, pcd));
					}
				}
			}
		}
		
		ArrayList<IntvPCden> tuning = new ArrayList<IntvPCden>(capacity);
		for (IntvPCden item: top2k) {
			tuning.add(new IntvPCden(item));
		}
		return tuning;
	}
	
//	private static void smooth(int l[], int u[], int h, int t, double xi) {
//		if (xi >= 3 * avgChange) xi = 3 * avgChange;
//		l[h] = t;
//		while (l[h] > 0 && Math.abs(cden[l[h] - 1] - cden[t]) < xi) { 
//			l[h]--; 
//			if (h > 0 && l[h] == u[h - 1]) break;
//		}
//		u[h] = t;
//		while (u[h] < T - 1 && Math.abs(cden[u[h] + 1] - cden[t]) < xi) { 
//			u[h]++; 
//		}
//	}
	
	private static void smoothLocal(int l[], int u[], int h, int t, boolean maxima) {
		l[h] = t;
		//double epn = maxima ? Math.max(cden[l[h]], 0) : Math.max(-cden[l[h]], 0);
		double epn = maxima ? cden[l[h]] : -cden[l[h]];
		double xi = Math.exp(epn) * localAvgChange[l[h]];
		if (xi >= 3 * avgChange) xi = 3 * avgChange;
		while (l[h] > 0 && Math.abs(cden[l[h] - 1] - cden[l[h]]) < xi) { 
			l[h]--; 
			//epn = maxima ? Math.max(cden[l[h]], 0) : Math.max(-cden[l[h]], 0);
			epn = maxima ? cden[l[h]] : -cden[l[h]];
			xi = Math.exp(epn) * localAvgChange[l[h]];
			if (xi >= 3 * avgChange) xi = 3 * avgChange;
			if (h > 0 && l[h] == u[h - 1]) break;
		}
		u[h] = t;
		//epn = maxima ? Math.max(cden[u[h]], 0) : Math.max(-cden[u[h]], 0);
		epn = maxima ? cden[u[h]] : -cden[u[h]];
		xi = Math.exp(epn) * localAvgChange[u[h]];
		if (xi >= 3 * avgChange) xi = 3 * avgChange;
		while (u[h] < T - 1 && Math.abs(cden[u[h] + 1] - cden[u[h]]) < xi) { 
			u[h]++; 
			//epn = maxima ? Math.max(cden[u[h]], 0) : Math.max(-cden[u[h]], 0);
			epn = maxima ? cden[u[h]] : -cden[u[h]];
			xi = Math.exp(epn) * localAvgChange[u[h]];
			if (xi >= 3 * avgChange) xi = 3 * avgChange;
		}
	}
	
	private static void earlyTuning(IntvPCden intv) {
		int base = 1;
		int incTime = 0;
		int[] s = new int[4], t = new int[4];
		double[] pcd = new double[4];
		int maxInd = -1;
		while (true) {
			boolean change = false;
			s[0] = (intv.s - base >= 0) ? intv.s - base : 0;
			t[0] = intv.t;
			pcd[0] = mapSearch(s[0], t[0]);
			s[1] = (intv.s + base <= intv.t) ? intv.s + base : intv.t;
			t[1] = intv.t;
			pcd[1] = mapSearch(s[1], t[1]);
			s[2] = intv.s;
			t[2] = (intv.t + base < tg.getNTimestamp()) ? intv.t + base : tg.getNTimestamp() - 1;
			pcd[2] = mapSearch(s[2], t[2]);
			s[3] = intv.s;
			t[3] = (intv.t - base >= intv.s) ? intv.t - base : intv.s;
			pcd[3] = mapSearch(s[3], t[3]);
			for (int i = 0; i < 4; i++) {
				if (pcd[i] >= pcd[0] && pcd[i] >= pcd[1] && pcd[i] >= pcd[2] && pcd[i] >= pcd[3]) {
					maxInd = i; break;
				}
			}
			if (pcd[maxInd] > intv.d) {
				//System.out.print(intv.s + "," + intv.t + "\t");
				intv.s = s[maxInd]; intv.t = t[maxInd]; intv.d = pcd[maxInd];
				change = true;
				//System.out.println(intv.s + "," + intv.t);
			}
			if (!change) break;
			incTime++;
			if (incTime == 4) {
				base *= 2;
				incTime = 0;
			}
		}
	}
	private static double mapSearch(final int s, final int t) {
		if (s > t) return 0;
		Edge key = new Edge(s, t);
		Double value = pcden.get(key);
		if (value == null) {
			value = tg.getPosCDen(s, t);
			pcden.put(key, value);
		}
		return value;		
	}
	
	/**
	 * Tune the interval intv by enlarging or shrinking the boundary.
	 * The interval is enlarged/shrunk as long as the sum of weights of selected edges is 
	 * positive. 
	 * @param tg temporal graph
	 * @param dsgEdge selectEdges
	 * @param intv the original interval
	 * @param tunedIntv the tuned interval
	 */
	public static double tuning(final TGraph tg, final TreeSet<Edge> dsgEdge, final Edge intv, 
			final double oriWeight, Edge tunedIntv) {
		// select Edge id of dsgEdge
		int[] eid = tg.computeEids(dsgEdge);
		
		HashMap<Edge, Double> density = new HashMap<Edge, Double>();
		
		tunedIntv.copy(intv);
		double finalWeight = oriWeight;
		int[] s = new int[4], t = new int[4];
		double[] pcd = new double[4];
		int maxInd = -1, base = 1, incTime = 0;
		while (true) {
			boolean change = false;
			s[0] = (tunedIntv.getS() - base >= 0) ? tunedIntv.getS() - base : 0;
			t[0] = tunedIntv.getT();
			s[1] = (tunedIntv.getS() + base <= tunedIntv.getT()) ? 
					tunedIntv.getS() + base : tunedIntv.getT();
			t[1] = tunedIntv.getT();
			s[2] = tunedIntv.getS();
			t[2] = (tunedIntv.getT() + base < tg.getNTimestamp()) ? 
					tunedIntv.getT() + base : tg.getNTimestamp() - 1;
			s[3] = tunedIntv.getS();
			t[3] = (tunedIntv.getT() - base >= tunedIntv.getS()) ? 
					tunedIntv.getT() - base : tunedIntv.getS();
					
			for (int i = 0; i < 4; i++) {
				Edge key = new Edge(s[i], t[i]);
				Double value = density.get(key);
				if (value == null) {
					value = tg.computeSumOfEdgeWeight(s[i], t[i], eid);
					density.put(key, value);
				}
				pcd[i] = value;
			}
			for (int i = 0; i < 4; i++) {
				if (pcd[i] >= pcd[0] && pcd[i] >= pcd[1] && pcd[i] >= pcd[2] && pcd[i] >= pcd[3]) {
					maxInd = i; break;
				}
			}
			if (pcd[maxInd] > finalWeight) {
				//System.out.print(tunedIntv.getS() + "," + tunedIntv.getT() + "\t");
				tunedIntv.setS(s[maxInd]); tunedIntv.setT(t[maxInd]); 
				finalWeight = pcd[maxInd];
				change = true;
				//System.out.println(tunedIntv.getS() + "," + tunedIntv.getT());
			}
			if (!change) break;
			incTime++;
			if (incTime == 4) {
				base *= 2;
				incTime = 0;
			}
		}
		
		return finalWeight;
	}
	
	private static void clear() {
		cden = null;
		localAvgChange = null;
//		int cnt = 0;
//		for (Entry<Edge, Double> item: pcden.entrySet()) {
//			Edge key = item.getKey();
//			System.out.print("[" + key.getS() + "," + key.getT() + "]" + item.getValue() + "\t");
//			cnt++;
//			if (cnt % 10 == 0) System.out.println();
//		}
		pcden = null;
	}
}

class IntvPCden implements Comparable<IntvPCden>{
	int s, t; double d;
	public IntvPCden(final int s, final int t, final double pcden) {
		this.s = s; this.t = t; this.d = pcden;
	}
	public IntvPCden(IntvPCden item) {
		this.s = item.s; this.t = item.t; this.d = item.d;
	}
	@Override // placing items of higher d earlier
	public int compareTo(IntvPCden o) {
		if (this.d < o.d) return 1;
		else if (this.d > o.d) return -1;
		return 0;
	}
}
