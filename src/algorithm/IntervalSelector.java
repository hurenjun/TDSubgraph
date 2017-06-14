package algorithm;

import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeSet;

import graph.Edge;
import graph.TGraph;

public class IntervalSelector {
	static int T, nNode;
	static double xiBase, xiSmooth;
	static double[] cden;
	/*
	 * member variable description: 
	 * 1) T: number of timestamps; nNode: number of node;
	 * 2) xi: smooth factor, based on average change of cden, 
	 *    i.e. xi = alpha * \Sigma(cden[i + 1] - cden[i]) / (T - 1), where alpha is fixed to 2;
	 * 3) cden[i]: cohesive density of time i, i.e. cden[i] = sum of edge weights in time i;
	 */
	
	/*
	 * DoIntervalSelect: select k time intervals according to 
	 * 1) half of the k intervals are select by maxTInterval, while another half are 
	 *    selected by minTInterval;
	 * 2) algorithm maxTInterval uses local maximums to select candidate intervals, 
	 *    and minTInterval on the other hand uses local minimums;
	 * 3) both algorithms select the k/2 interval from candidates according the positive
	 *    cohesive density of each interval;
	 */
	public static Edge[] doIntervalSelect(final TGraph tg, final int k, final double xi) {
		T = tg.getNTimestamp();
		nNode = tg.getNNode();
		cden = tg.getCDensity(); // cohesive density
		
		for (int i = 0; i < T; i++) {
			System.out.print(i + ":" + cden[i] + "\t");
			if (i % 10 == 9) System.out.println();
		}
		System.out.println();
		
		xiBase = 0;	
		for (int i = 0; i < T - 1; i++) {
			xiBase += Math.abs(cden[i + 1] - cden[i]);
		}
		if (T > 1) xiBase /= (T - 1);	// average change of cden
		xiSmooth = xiBase * xi;
		
		Edge[] topK = new Edge[k];
		for (int i = 0; i < k; i++) { topK[i] = new Edge();}
		double[] posCDen = new double[k];
		Arrays.fill(posCDen, -1);
		
		TreeSet<Edge> candIntvMax = maxTInterval();
		//System.out.println("candIntvMax: " + candIntvMax.size());
		Iterator<Edge> iter = candIntvMax.iterator();
		while (iter.hasNext()) {
			Edge e = (Edge)iter.next();
			double pcden = tg.getPosCDen(e.getS(), e.getT());
			if (pcden > posCDen[k / 2 - 1]) {
				topK[k / 2 - 1].copy(e);
				posCDen[k / 2 - 1] = pcden;
				int i = k / 2 - 1;
				while (i > 0 && posCDen[i - 1] < posCDen[i]) {
					Edge tempE = new Edge(topK[i - 1]);
					topK[i - 1].copy(topK[i]);
					topK[i].copy(tempE);
					double tempPcden = posCDen[i - 1];
					posCDen[i - 1] = posCDen[i];
					posCDen[i] = tempPcden;
					i--;
				}
			}
		}
		for (int i = 0; i < k / 2; i++) {
			if (topK[i].equals(new Edge())) continue;
			int base = 1;
			int incTime = 0;
			while (true) {
				boolean change = false;
				Edge e = topK[i];
				double pcden = tg.getPosCDen(e.getS() - base, e.getT());
				if (pcden > posCDen[i]) {
					e.setS(e.getS() - base >= 0 ? e.getS() - base : 0);
					posCDen[i] = pcden;
					change = true;
				}
				pcden = tg.getPosCDen(e.getS(), e.getT() + base);
				if (pcden > posCDen[i]) {
					e.setT(e.getT() + base < tg.getNTimestamp() ? e.getT() + base : tg.getNTimestamp() - 1);
					posCDen[i] = pcden;
					change = true;
				}
				if (!change) break;
				incTime++;
				if (incTime == 4) {
					base *= 2;
					incTime = 0;
				}
			}
		}
		
		TreeSet<Edge> candIntvMin = minTInterval();
		//System.out.println("candIntvMin: " + candIntvMin.size());
		iter = candIntvMin.iterator();
		while (iter.hasNext()) {
			Edge e = (Edge)iter.next();
			double pcden = tg.getPosCDen(e.getS(), e.getT());
			if (pcden > posCDen[k - 1]) {
				topK[k - 1].copy(e);
				posCDen[k - 1] = pcden;
				int i = k - 1;
				while (i > k / 2 && posCDen[i - 1] < posCDen[i]) {
					Edge tempE = new Edge(topK[i - 1]);
					topK[i - 1].copy(topK[i]);
					topK[i].copy(tempE);
					double tempPcden = posCDen[i - 1];
					posCDen[i - 1] = posCDen[i];
					posCDen[i] = tempPcden;
					i--;
				}
			}
		}
		for (int i = k / 2; i < k; i++) {
			if (topK[i].equals(new Edge())) continue;
			int base = 1, incTime = 0;
			while (true) {
				boolean change = false;
				Edge e = topK[i];
				double pcden = tg.getPosCDen(e.getS() + base, e.getT());
				if (pcden > posCDen[i]) {
					e.setS(e.getS() + base);
					posCDen[i] = pcden;
					change = true;
				}
				pcden = tg.getPosCDen(e.getS(), e.getT() - base);
				if (pcden > posCDen[i]) {
					e.setT(e.getT() - base);
					posCDen[i] = pcden;
					change = true;
				}
				if (!change) break;
				incTime++;
				if (incTime == 4) {
					base *= 2;
					incTime = 0;
				}
			}
		}
		
		clear();
		return topK;
	}
	
	// select candidate intervals using local maximums
	private static TreeSet<Edge> maxTInterval() {
		if (T <= 10) {
			TreeSet<Edge> intervals = new TreeSet<Edge>();
			for (int i = 0; i < T; i++) {
				for (int j = i; j < T; j++) {
					intervals.add(new Edge(i, j));
				}
			}
			return intervals;
		}
		
		int[] l = new int[nNode];	// lower sides of maximums
		int[] u = new int[nNode];	// upper sides of maximums
		int h = 0;
		//System.out.print("local maxima: ");
		if (cden[0] >= cden[1]) {	// consider the special case, 0
			smooth(l, u, h, 0, xiSmooth);
			h++;
			//System.out.print(0 + " ");
		}
		
		for (int i = 1; i < T - 1; i++) {
			if (cden[i] > cden[i - 1] && cden[i] >= cden[i + 1] || 
					cden[i] >= cden[i - 1] && cden[i] > cden[i + 1]) {
				smooth(l, u, h, i, xiSmooth);
				h++;
				//System.out.print(i + " ");
			}
		}
		
		if (T >= 2 && cden[T - 1] >= cden[T - 2]) {
			smooth(l, u, h, T - 1, xiSmooth);
			h++;
			//System.out.print(T - 1 + " ");
		}
		//System.out.println();
		
		TreeSet<Edge> intervals = new TreeSet<Edge>();
		for (int i = 0; i < h; i++) {
			intervals.add(new Edge(l[i], u[i]));
		}
		h = 0;
		Iterator<Edge> iter = intervals.iterator();
		Edge combIntv = null;
		while (iter.hasNext()) {
			Edge e = (Edge)iter.next();
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
//		System.out.print("maxTInterval: " + h + "\t");
//		for (int i = 0; i < h; i++) {
//			System.out.print("[" + l[i] + "," + u[i] + "]\t");
//		}
//		System.out.println();
//		for (int i = 0; i < h; i++) {
//			double base = 0;
//			for (int j = l[i]; j <= u[i]; j++) {
//				base += cden[j];
//			}
//			base /= (u[i] - l[i] + 1);
//			smooth(l, u, i, base, xiSmooth);
//		}
		
		intervals.clear();
		for (int i = 0; i < h; i++) {
			for (int j = i; j < h; j++) {
				if (l[i] <= u[j]) {
					intervals.add(new Edge(l[i], u[j]));
				}
			}
		}
		return intervals;
	}
	
	// select candidate intervals using local minimums
	private static TreeSet<Edge> minTInterval() {
		if (T <= 10) {
			TreeSet<Edge> intervals = new TreeSet<Edge>();
			for (int i = 0; i < T; i++) {
				for (int j = i; j < T; j++) {
					intervals.add(new Edge(i, j));
				}
			}
			return intervals;
		}
		
		int[] l = new int[nNode];	// lower sides of maximums
		int[] u = new int[nNode];	// upper sides of maximums
		int h = 0;
		//System.out.print("local minima: ");
		if (cden[0] <= cden[1]) {	// consider the special case, 0
			smooth(l, u, h, 0, xiSmooth);
			h++;
			//System.out.print(0 + " ");
		}
		
		for (int i = 1; i < T - 1; i++) {
			if (cden[i] < cden[i - 1] && cden[i] <= cden[i + 1] || 
					cden[i] <= cden[i - 1] && cden[i] < cden[i + 1]) {
				smooth(l, u, h, i, xiSmooth);
				h++;
				//System.out.print(i + " ");
			}
		}
		
		if (T >= 2 && cden[T - 1] <= cden[T - 2]) {
			smooth(l, u, h, T - 1, xiSmooth);
			h++;
			//System.out.print(T - 1 + " ");
		}
		//System.out.println();
		
		TreeSet<Edge> intervals = new TreeSet<Edge>();
		for (int i = 0; i < h; i++) {
			intervals.add(new Edge(l[i], u[i]));
		}
		h = 0;
		Iterator<Edge> iter = intervals.iterator();
		Edge combIntv = null;
		while (iter.hasNext()) {
			Edge e = (Edge)iter.next();
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
		
		intervals.clear();
		for (int i = 0; i < h; i++) {
			for (int j = i + 1; j < h; j++) {
				if (u[i] <= l[j]) {
					intervals.add(new Edge(u[i], l[j]));
				}
			}
		}
		return intervals;
	}
	
	private static void smooth(int l[], int u[], int h, int t, double xi) {
		l[h] = t;
		while (l[h] > 0 && Math.abs(cden[l[h] - 1] - cden[t]) < xi) { l[h]--; }
		u[h] = t;
		while (u[h] < T - 1 && Math.abs(cden[u[h] + 1] - cden[t]) < xi) { u[h]++; }
	}
	
	private static void clear() {
		cden = null;
	}
}
