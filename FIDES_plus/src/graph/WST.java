package graph;

 public class WST implements Comparable<WST> {
	double w;
	int s, t;
	public WST(final double w, final int s, final int t) {
		this.w = w; this.s = s; this.t = t;
	}
	
	public WST(final WST o) {
		w = o.w; s = o.s; t = o.t;
	}

	public double getW() {
		return w;
	}

	public int getS() {
		return s;
	}

	public int getT() {
		return t;
	}
	
	public boolean equalST(final int s, final int t) {
		if (this.s == s && this.t == t || this.s == t && this.t == s)
			return true;
		else return false;
	}

	@Override
	public int compareTo(WST o) {
		if (w < o.w) return -1;
		else if (w > o.w) return 1;
		if (s < o.s || s == o.s && t < o.t) return -1;
		else if (s == o.s && t == o.t) return 0;
		else return 1;
	}
	
	
}