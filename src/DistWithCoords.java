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
	public int index;				//this is the map index for back reference
	public long dist;				//this is the distance from the source in the graph
	public long distR;				//this is the distance from the target in the reverse graph
	public double k;				//k(v) = d(s,v) + pf(v)   actual distance from s plus potential function to t
	public double kr;				//kr(v) = dr(s,v) + pf(v)   actual distance from s plus potential function to t
	public boolean processed;
	public boolean processedR;
	public boolean queued;
	public boolean queuedR;
	public int pindex;
	public Point coord;
	
	public Node(int i, Point coord) {
        this.index = i;
        this.dist = INFINITY;
        this.distR = INFINITY;
        this.processed = false;
        this.processedR = false;
        this.queued = false;
        this.queuedR = false;
        this.pindex = -1;
        this.coord = coord;
	}
	
	public void resetNode() {
		this.dist = INFINITY;
        this.distR = INFINITY;
        this.processed = false;
        this.processedR = false;
        this.queued = false;
        this.queuedR = false;
	}
	
	private double pf(Node start, Node finish) {
		double pf = (coord.distance(finish.coord) - coord.distance(start.coord))/2 + finish.coord.distance(start.coord)/2;
		return pf;
	}

	public double setK(Node start, Node finish) {
		k = dist + pf(start, finish);
		return k;
	}
	
	public double setKr(Node start, Node finish) {
		kr = distR + pf(finish, start);					//start = finish for reverse graph potential
		return kr;
	}
}

class Edge {
	public int target;
	public int source;
	public long length;
	public double lengthP;
	public Node u;
	public Node v;
	
	
	public Edge(Node u, Node v, int l) {
		this.target = u.index;
		this.source = v.index;
		this.u = u;
		this.v = v;
		this.length = l;
		this.lengthP = this.length;	//initially this is the same
	}

	//this method is for BiDijkstra with modified edges, conflicts with grader (double vs long)
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
	public long mu = INFINITY;

	
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
    		if(left <= li && h.get(left).k < h.get(mini).k) {
        		mini = left;
        	}
        	if(rght <= li && h.get(rght).k < h.get(mini).k) {
        		mini = rght;
        	}
        	if(mini != i) {
    	   		swap(i, mini, h);
    	   		minSiftDown(mini, h);
        	}
    	} else {
    		if(left <= li && h.get(left).kr < h.get(mini).kr) {
        		mini = left;
        	}
        	if(rght <= li && h.get(rght).kr < h.get(mini).kr) {
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
				if(h.get(pi).k > h.get(i).k) {
					swap(pi, i, h);
					minSiftUp(pi, h);
				}
			} else {
				if(h.get(pi).kr > h.get(i).kr) {
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
	
	private void decreaseKey(Node dn, long d, Node start, Node finish, ArrayList<Node> h) {
		
//		System.out.println(" decreaseKey i: " + i + " d: " + d);
		
		if(h == heapR) {
			dn.distR = d;
			dn.setKr(start, finish);
		} else {
			dn.dist = d;
			dn.setK(start, finish);
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
	
	public long biAStar(int s, int t) {	//s & t are 0 indexed inegers from the graph data
		
		//implements  bidirectional A* algorithm
//		long start = System.nanoTime();	
//		int processed = 0;
		
		double prt;
		long mu = INFINITY;
		
		initializeQueues();
		resetWorkingNodes();
		Node sn = map.get(s);
		Node tn = map.get(t);
		
		sn.dist = 0;
		sn.setK(sn, tn);
		
		tn.distR = 0;
		tn.setKr(sn, tn);
		
		prt = tn.kr;									//reverse potential of target since tn.distr = 0;
		
		enQueue(sn, heap);
		sn.queued = true;
		enQueue(tn, heapR);
		tn.queuedR = true;	
		working.add(sn);								//add the initial nodes to the working set
		working.add(tn);	
		long result = -1;								//if after processing all nodes in the graph this is unchanged the target is unreachable
		
		
		if(s == t) {
			return 0;									//I found myself!!
		}
		
		while(!heap.isEmpty() || !heapR.isEmpty()) {
	
			if(!heap.isEmpty()) {						//process the next node in the forward graph
				
				Node processing = getMin(heap);
				processing.queued = false;
				
//				System.out.println("processing Node: " + processing.index);
				
				for(Edge e : graph.get(processing.index)) {
					
					Node tt = e.v;
					
					long td = processing.dist + e.length;
				
					if(tt.dist > td) {
							
							working.add(tt);
							
							if(tt.queued == false) {
								tt.dist = td;
								tt.setK(sn, tn);
								enQueue(tt, heap);
								tt.queued = true;
							} else {
								decreaseKey(tt, td, sn, tn, heap);
							}

//						tt.pindex = processing.index;				//min path is my daddy
					}
//					++processed;
				}
				
				processing.processed = true;
				
				if(processing.processedR == true) {				
					long tp = processing.dist + processing.distR;
					if( tp < mu) {
						mu = tp; 
						result = mu;
					}
					if(!heap.isEmpty()) {
						if(heap.get(0).k + heapR.get(0).kr > mu + prt) {	//stop when the shortest queued nodes have estimated paths longer than the shortest one found
							break;
						}
					}

										
				}
					
			}
			
			if(!heapR.isEmpty()) {						//process the next node in the reverse graph
				
				Node processingR = getMin(heapR);
				processingR.queuedR = false;
				
//				System.out.println("processing Node: " + rr.index);
				
				for(Edge er : graphR.get(processingR.index)) {
			
					Node ttr = er.v;
					
					long tdr = processingR.distR + er.length;
				
					if(ttr.distR > tdr) {

							working.add(ttr);
							if(ttr.queuedR == false) {
								ttr.distR = tdr;
								ttr.setKr(sn, tn);
								enQueue(ttr, heapR);
								ttr.queuedR = true;
							} else {
								decreaseKey(ttr, tdr, sn, tn, heapR);
							}
					
//						ttr.pindex = rr.index;
					}
//					++processed;
				}
				processingR.processedR = true;
				if(processingR.processed == true) {
					long tp = processingR.dist + processingR.distR;
					if( tp < mu) {
						mu = tp;
						result = mu;
					}
					if(!heapR.isEmpty()) {
						if(heap.get(0).k + heapR.get(0).kr > mu + prt) {
							break;
						}
					}

					
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
            System.out.println(g.biAStar(u-1, v-1));
        }
        in.close();
    }
}
