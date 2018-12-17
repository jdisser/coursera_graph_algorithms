import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashSet;
import java.lang.Math;

class Point {
	public double x;
	public double y;
	
	public Point(int x, int y) {
		this.x = (double) x;
		this.y = (double) y;
	}
	
	public double distance(Point t) {
		return Math.sqrt(Math.pow((x - t.x), 2.0) + Math.pow((y - t.y), 2.0));
	}
}


class Node {
	private static final long INFINITY = Long.MAX_VALUE/4;
	public int index;
	public double dist;				//this is the potential distance from the source in the graph
	public double distR;			//this is the potential distance from the target in the reverse graph
	public boolean visited;
	public boolean visitedR;
	public int pindex;
	public Point coord;
	
	public Node(int i, Point coord) {
        this.index = i;
        this.dist = INFINITY;
        this.distR = INFINITY;
        this.visited = false;
        this.visitedR = false;
        this.pindex = -1;
        this.coord = coord;
	}
	
	public void resetNode() {
		this.dist = INFINITY;
        this.distR = INFINITY;
        this.visited = false;
        this.visitedR = false;
	}		
}

class Edge {
	public int target;
	public int source;
	public double length;
	public double lengthP;
	public Node u;
	public Node v;
	
	
	public Edge(Node u, Node v, int l) {
		this.target = u.index;
		this.source = v.index;
		this.u = u;
		this.v = v;
		this.length = (double)l;
		this.lengthP = this.length;	//initially this is the same
	}

	
	public void setLengthP(Node start, Node finish) { 		//parameters reference the current search directions (start = target in reverse)
		double pfu = u.coord.distance(finish.coord);
		double pfv = v.coord.distance(finish.coord);
		double pru = u.coord.distance(start.coord);
		double prv = v.coord.distance(start.coord);
		
		double pafu = (pfu - pru)/2;
		double pafv = (pfv - prv)/2;
		
		lengthP = length - pafu + pafv;
		
	}
	
	
}

class BiGraph {
	
	public ArrayList<ArrayList<Edge>> graph;		//adjacency list of edges
	public ArrayList<ArrayList<Edge>> graphR;		//adjacency list of reversed edges
	public ArrayList<Node> heap;					//priority queue with method to update key values
	public ArrayList<Node> heapR;					//priority queue for reversed graph
	public ArrayList<Node> map;						//returns Nodes by vertex (index)
	public HashSet<Node> working;					//all nodes processed by query that must be reset
	
	public final long INFINITY = Long.MAX_VALUE / 4;


	
	public int root;
	public int n;

	
	public BiGraph(int n) {
		this.n = n;

        this.graph = new ArrayList<ArrayList<Edge>>();
        this.graphR = new ArrayList<ArrayList<Edge>>();
        
        this.heap = new ArrayList<Node>();
        this.heapR = new ArrayList<Node>();
        
        this.map = new ArrayList<Node>();
        
        this.working = new HashSet<Node>();
        
        for (int i = 0; i < n; i++) {					//graph indexes and nodes are 0 indexed
        	
            this.graph.add(new ArrayList<Edge>());
            this.graphR.add(new ArrayList<Edge>());          
  
        }
        
        this.root = -1;

	}
	
	public void addNode(int i, Point p) {
		Node n = new Node(i,p);
		map.add(i, n);
	}
	
	public void initializeQueues() {
		
		heap.clear();
		heapR.clear();
		
	}
	
	public void resetWorkingNodes() {
		if(!working.isEmpty()) {
			for(Node n: working) {
				n.resetNode();
			}
			working.clear();
		}
	}
	
	public void enQueue(Node n, ArrayList<Node> h) {
		
		int end = h.size();
		h.add(end, n);
		minSiftUp(end, h);
//		working.add(n);
		
	}
	
	
	public void addEdges(int s, int t, int c) {		//(source, target, length in 1 based indexing)

		s -= 1;
		t -= 1;								//raw vertexes are 1 based, map indices are 0 based
		
		Node source = map.get(s);
		Node target = map.get(t);
		
		
		Edge e = new Edge(source, target, c);
        Edge er = new Edge(target, source, c);

        graph.get(s).add(e);
		graphR.get(t).add(er);

	}
	
	
	private void minSiftDown(int i, ArrayList<Node> h) {
		
//		System.out.println("minSiftDown 1: " + i);
//		printHeap();
		
    	int mini = i;
    	int li = h.size() - 1;	//0 based last index
    	int left = 2*i + 1;
    	int rght = 2*i + 2;
    	
    	if(h == heap) {
    		if(left <= li && h.get(left).dist < h.get(mini).dist) {
        		mini = left;
        	}
        	if(rght <= li && h.get(rght).dist < h.get(mini).dist) {
        		mini = rght;
        	}
        	if(mini != i) {
    	   		swap(i, mini, h);
    	   		minSiftDown(mini, h);
        	}
    	} else {
    		if(left <= li && h.get(left).distR < h.get(mini).distR) {
        		mini = left;
        	}
        	if(rght <= li && h.get(rght).distR < h.get(mini).distR) {
        		mini = rght;
        	}
        	if(mini != i) {
    	   		swap(i, mini, h);
    	   		minSiftDown(mini, h);
        	}
    	}
 
    }
	
	private void printHeap(ArrayList<Node> h) {
		System.out.print("heap: ");
		for(Node n : h) {
			System.out.print(n.index + ":" + n.dist + " ");
		}
		System.out.println(" ");
	}
	
	
	
	private void minSiftUp(int i, ArrayList<Node> h) {
		
		int pi;
//		System.out.println(" minSiftUp i: " + i);
//		printHeap(h);
		if(i > 0) {
			pi = (i - 1)/2;
			
//			System.out.println(" minSiftUp pi: " + pi);
			if(h == heap) {
				if(h.get(pi).dist > h.get(i).dist) {
					swap(pi, i, h);
					minSiftUp(pi, h);
				}
			} else {
				if(h.get(pi).distR > h.get(i).distR) {
					swap(pi, i, h);
					minSiftUp(pi, h);
				}
			}
			
			
		}

	}
	
	private Node getMin(ArrayList<Node> h) {
		
		Node nm = h.get(0);
		swap(0, h.size() - 1, h);
		h.remove(h.size() - 1);
		minSiftDown(0, h);
		
		return nm;
	}
	
	private void decreaseKey(int i, long d, ArrayList<Node> h) {
		
//		System.out.println(" decreaseKey i: " + i + " d: " + d);
		Node dn = map.get(i);
		if(h == heapR) {
			dn.distR = d;
		} else {
			dn.dist = d;
		}
		int dni = h.indexOf(dn);
		minSiftUp(dni, h);
	}

    
    private void swap(int i, int n, ArrayList<Node> h) {
    	
    	Node tn = h.get(i);
    	h.set(i, h.get(n));
    	h.set(n, tn);
    	
    }
    
    private double shortestPath(double cDist) {
    	double d = INFINITY;
    	double dn = 0;
    	for(Node n: working) {
    		dn = n.dist + n.distR;
    		d = Math.min(d, dn);
    	}
    	return d + cDist;					//convert the potential distance into the real distance
    }
	
	public double biDijkstra(int s, int t) {	//s & t are 0 indexed
		
		//implements meet in the middle bidirectional Dijkstra's algorithm

//		long start = System.nanoTime();
		
//		int processed = 0;
		
		initializeQueues();
		
		resetWorkingNodes();
		
		Node sn = map.get(s);
		Node tn = map.get(t);
		
		enQueue(sn, heap);
		enQueue(tn, heapR);
		
		decreaseKey(s, 0, heap);				//start the algorithm on this node in the forward direction
		decreaseKey(t, 0, heapR);				//and from this one in the reverse direction (Reverse Graph => heapR)
		
		
		working.add(sn);				//add the initial nodes to the working set
		working.add(tn);
		
		
		double result = -1;						//if after processing all nodes in the graph this is unchanged the target is unreachable
		
		double cDist = sn.coord.distance(tn.coord);		//shortest distance adjustment
		
		if(s == t) {
			return 0;							//I found myself!!
		}
		
		while(!heap.isEmpty() || !heapR.isEmpty()) {
	
			if(!heap.isEmpty()) {						//process the next node in the forward graph
				
				Node r = getMin(heap);
//				System.out.println("processing Node: " + r.index);
				for(Edge e : graph.get(r.index)) {
					
					e.setLengthP(sn, tn);				//Set the potential edge length in the forward direction
					
					Node tt = e.v;
				
					if(tt.dist > r.dist + e.lengthP) {

							working.add(tt);
							tt.dist = r.dist + e.lengthP;
							enQueue(tt, heap);

						
//						tt.pindex = r.index;				//min path is my daddy
					}
//					++processed;
				}
		
				r.visited = true;
				if(r.visitedR == true) {				//stop when a node has been processed from both ends
				
					result = shortestPath(cDist);
					break;					
				}
					
			}
			
			if(!heapR.isEmpty()) {						//process the next node in the reverse graph
				
				Node rr = getMin(heapR);	
//				System.out.println("processing Node: " + rr.index);
				
				for(Edge er : graphR.get(rr.index)) {
					
					er.setLengthP(tn, sn);				//start is target in Gr
					
					Node ttr = er.v;
				
					if(ttr.distR > rr.distR + er.lengthP) {
															

							working.add(ttr);
							ttr.distR = rr.distR + er.lengthP;
							enQueue(ttr, heapR);

						
//						ttr.pindex = rr.index;
					}
//					++processed;
				}
				rr.visitedR = true;
				if(rr.visited == true) {
				
					result = shortestPath(cDist);
					break;
				}
			}
			
		}
//		long finish = System.nanoTime();
//		long elapsed = (finish - start) / 1000000;
//		System.out.println("BiDijkstra processed edges: " + processed + " ms: " + elapsed);
		return result;
    }

}



public class DistWithCoords {
    

    public static void main(String args[]) {
        Scanner in = new Scanner(System.in);
        int n = in.nextInt();
        int m = in.nextInt();
        
        BiGraph g = new BiGraph(n);

        for (int i = 0; i < n; i++) { 
            int x, y;
            x = in.nextInt();
            y = in.nextInt();
            
            Point p = new Point(x , y);
            g.addNode(i, p);

        }

        for (int i = 0; i < m; i++) {
            int x, y, c;
            x = in.nextInt();
            y = in.nextInt();
            c = in.nextInt();
            
            g.addEdges(x, y, c);
            
 
        }

        int t = in.nextInt();

        for (int i = 0; i < t; i++) {
            int u, v;
            u = in.nextInt();
            v = in.nextInt();
            System.out.println((long)g.biDijkstra(u-1, v-1));
        }
        in.close();
    }
}
