package graph;

 public class WST {
	double w;
	int s, t;
	public WST(final double _w, final int _s, final int _t) {
		w = _w; s = _s; t = _t;
	}
	
	public WST(final WST other) {
		w = other.w; s = other.s; t = other.t;
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
	
	
}