package graph;

// maintain an edge with two endpoints s and t;
public class Edge implements Comparable<Edge> {
	private int s, t;
	
	public Edge() {
		s = t = -1;
	}
	
	public Edge(final int _s, final int _t) {
		s = _s;
		t = _t;
	}
	
	public Edge(final Edge e) {
		s = e.s;
		t = e.t;
	}
	
	public int getS() {
		return s;
	}

	public int getT() {
		return t;
	}
	
	public void setS(final int _s) {
		s = _s;
	}
	
	public void setT(final int _t) {
		t = _t;
	}
	
	public void copy(final Edge e) {
		s = e.s;
		t = e.t;
	}

	@Override
	public int compareTo(final Edge edgeO) {
		if (s < edgeO.s || s == edgeO.s && t < edgeO.t) return -1;
		else if (s == edgeO.s && t == edgeO.t) return 0;
		else return 1;
	}
	
	public boolean equals(final Edge edgeO) 
    { 
        if (s == edgeO.s && t == edgeO.t) return true;
        return false;
    } 
}
