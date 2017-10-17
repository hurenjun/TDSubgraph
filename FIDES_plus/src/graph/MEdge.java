package graph;

public class MEdge implements Comparable<MEdge> {
	int oriS, oriT, s, t;
	double w;
	
	public MEdge(final int s, final int t, final int oriS, final int oriT, final double w) {
		this.s = s; this.t = t;
		this.oriS = oriS; this.oriT = oriT;
		this.w = w;
	}

	public MEdge(final MEdge me) {
		s = me.s; t = me.t;
		oriS = me.oriS; oriT = me.oriT;
		w = me.w;
	}
	
	public int getS() {
		return s;
	}

	public int getT() {
		return t;
	}
	
	public int getOriS() {
		return oriS;
	}
	
	public int getOriT() {
		return oriT;
	}
	
	public double getW() {
		return w;
	}
	
	@Override
	public int compareTo(final MEdge edgeO) {
		if (s < edgeO.s || s == edgeO.s && t < edgeO.t) return -1;
		else if (s == edgeO.s && t == edgeO.t) return 0;
		else return 1;
	}
	
	public boolean equals(final MEdge edgeO) 
    { 
        if (s == edgeO.s && t == edgeO.t) return true;
        return false;
    } 
	
	public boolean equalST(final int s, final int t) {
		if (this.s == s && this.t == t || this.s == t && this.t == s)
			return true;
		else return false;
	}
}
