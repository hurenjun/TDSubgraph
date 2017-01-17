package graph;

public class MEdge implements Comparable<MEdge> {
	int oriS, oriT, s, t;
	
	public MEdge(final int _s, final int _t, final int _oriS, final int _oriT) {
		s = _s;
		t = _t;
		oriS = _oriS;
		oriT = _oriT;
	}

	public MEdge(final MEdge me) {
		s = me.s;
		t = me.t;
		oriS = me.oriS;
		oriT = me.oriT;
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
	
}
