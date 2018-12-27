import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.io.BufferedReader;
import java.io.IOException;
import java.lang.Math;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

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

class Barrier {
	public Point ll;
	public Point ur;
	
	public Barrier(int llx, int lly, int urx, int ury, int L) {
		this.ll = new Point((llx - 1) * L,(lly - 1) * L);
		this.ur = new Point((urx - 1) * L, (ury - 1) * L);
	}
	
	public boolean contains(Point px) {
//		System.out.println("Contains: p(" + px.x + "," + px.y + ") ll(" + ll.x + "," + ll.y + ") ur(" + ur.x + "," + ur.y + ")");
		if(px.x <= ur.x && px.x >= ll.x && px.y >= ll.y && px.y <= ur.y)
			return true;
		else
			return false;
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
	public boolean active;
	public int pindex;
	public int pindexR;
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
        this.pindexR  = -1;
        this.coord = coord;
        this.active = true;
	}
	
	public void resetNode() {
		this.dist = INFINITY;
        this.distR = INFINITY;
        this.processed = false;
        this.processedR = false;
        this.queued = false;
        this.queuedR = false;
        this.pindex = -1;
        this.pindexR = -1;
	}
	
	private double pf(Node start, Node finish) {
		double pf = (coord.distance(finish.coord) - coord.distance(start.coord))/2 + finish.coord.distance(start.coord)/2;
		return pf;
	}

	public double setK(Node start, Node finish) {
		k = dist + pf(start, finish);
		return k;
	}
	
	public double setK(Node start, Node finish, boolean astar) {
		if(astar)
			k = dist + pf(start, finish);
		else
			k = dist;
		return k;
	}
	
	public double setKr(Node start, Node finish) {
		kr = distR + pf(finish, start);					//start = finish for reverse graph potential
		return kr;
	}
	
	public double setKr(Node start, Node finish, boolean astar) {
		
		if(astar)
			kr = distR + pf(finish, start);					//start = finish for reverse graph potential
		else
			kr = distR;
		
		return kr;
	}
	
	public static int nodeNumber(int i, int j, int w) {	//1 index position of node in test graph w x h
		return w*(j - 1) + i;
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
	public HashMap<Integer,Node> map;				//returns Nodes by vertex (index)
	public HashSet<Node> working;					//all nodes processed by query that must be reset
	
	public ArrayList<Barrier> barriers;				//for generating sets of inactive nodes for test graphs
	
	public final long INFINITY = Long.MAX_VALUE / 4;


	
	public int root;
	public int cNode;
	public int n;
	public long mu = INFINITY;
	public int processed;

	
	public BiGraph(int n) {
		this.n = n;

        this.graph = new ArrayList<ArrayList<Edge>>();
        this.graphR = new ArrayList<ArrayList<Edge>>();
        
        this.heap = new ArrayList<Node>();
        this.heapR = new ArrayList<Node>();
        
        this.map = new HashMap<Integer,Node>(n+1);
        
        this.working = new HashSet<Node>();
        
        for (int i = 0; i <= n; i++) {					//graph indexes and nodes are 1 indexed, index 0 not used
        	
            this.graph.add(new ArrayList<Edge>());
            this.graphR.add(new ArrayList<Edge>());          
  
        }
        
        this.root = -1;
        this.cNode = -1;
        
        this.processed = 0;
        
        this.barriers = new ArrayList<Barrier>();

	}
	
	
	
	public void genTestGraph(int W, int H, int L, int C) {

		int nodes = 0;
		int edges = 0;
		int inactive = 0;
		
		//generate a complete grid of active and inactive nodes
		for(int i = 1; i <= W; ++i) {
			for(int j = 1; j <= H; ++j) {
				Point p = new Point((i - 1) * L, (j - 1) * L);
				Node n = addNode(Node.nodeNumber(i,j,W),p);
				for(Barrier b: barriers) {
					if(b.contains(p)) {
						n.active = false;
						++inactive;
					}
				}
				++nodes;
			}
		}
		
		//generate the set of edges connecting the adjacent active nodes
		//note graph will be directed but bidirectional with different distances in each direction
		for(int i = 1; i <= W; ++ i) {
			for(int j = 1; j <= H; ++j) {
				Node n = map.get(Node.nodeNumber(i, j, W));
				if(n.active) {
					for(int ii = -1; ii <= 1; ++ii) {
						for(int jj = -1; jj <= 1; ++jj) {
							if(jj == 0 && ii == 0)	//no self loops
								continue;
							if(i + ii > 0 && j + jj > 0 && i + ii <= W && j + jj <= H) {
								Node t = map.get(Node.nodeNumber(i + ii, j + jj, W));
								if(t.active) {
									if(ii * jj == 0) {
										addEdges(n.index,t.index,L + (int)Math.ceil(Math.random() * L / C));
										++edges;
									} else {
										addEdges(n.index,t.index,(L * 3) / 2 + (int)Math.ceil(Math.random() * L / C));
										++edges;
									}
									
								}
							}
						}
					}
				}
			}
		}
		System.out.println("Generated edges: " + edges + " nodes: " + nodes + " inactive: " + inactive);	
	}
	
	public void addBarrier(int llx, int lly, int urx, int ury, int L) {
		Barrier b = new Barrier(llx,lly,urx,ury,L);
		barriers.add(b);
	}
	
	public Node addNode(int i, Point p) {
		Node n = new Node(i,p);
		map.put(i, n);
		return n;
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

// map changed to hashmap now uses 1 based indexing
//		s -= 1;
//		t -= 1;								//raw vertexes are 1 based, map indices are 0 based
		
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
	
	private void decreaseKey(Node dn, long d, Node start, Node finish, ArrayList<Node> h, boolean astar) {
		
//		System.out.println(" decreaseKey i: " + i + " d: " + d);
		
		if(h == heapR) {
			dn.distR = d;
			if(astar)
				dn.setKr(start, finish);
			else
				dn.setKr(start, finish, false);
		} else {
			dn.dist = d;
			if(astar)
				dn.setK(start, finish);
			else
				dn.setK(start, finish, false);
		}
		int dni = h.indexOf(dn);
		minSiftUp(dni, h);
	}
	
	private void decreaseKey(Node dn, long d, Node start, Node finish, ArrayList<Node> h, boolean astar, boolean biDirectional) {
		
//		System.out.println(" decreaseKey i: " + i + " d: " + d);
		if(biDirectional) {
			if(h == heapR) {
				dn.distR = d;
				if(astar)
					dn.setKr(start, finish);
				else
					dn.setKr(start, finish, false);
			} else {
				dn.dist = d;
				if(astar)
					dn.setK(start, finish);
				else
					dn.setK(start, finish, false);
			}
			int dni = h.indexOf(dn);
			minSiftUp(dni, h);
		} else {
			h = heap;							//if not bidirectional default to the forward direction
			dn.dist = d;
			if(astar)
				dn.setK(start, finish);
			else
				dn.setK(start, finish, false);
			int dni = h.indexOf(dn);
			minSiftUp(dni, h);
		}
		
	}

    
    private void swap(int i, int n, ArrayList<Node> h) {
    	
    	Node tn = h.get(i);
    	h.set(i, h.get(n));
    	h.set(n, tn);
    	
    }
    
    private long shortestPath() {
    	long d = INFINITY;
    	long dn = 0;
    	for(Node n: working) {
    		dn = n.dist + n.distR;
    		d = Math.min(d, dn);
    	}
    	return d;					
    }
	
	public long biAStar(int s, int t) {	//s & t are 0 indexed integers from the graph data
		
		//implements  bidirectional A* algorithm
//		long start = System.nanoTime();	
		processed = 0;
		int dblProcessed = 0;
		
		double prt;
		long mu = INFINITY;
		
		initializeQueues();
		resetWorkingNodes();
		
		cNode = -1;
		
		Node sn = map.get(s);
		Node tn = map.get(t);
		
		sn.dist = 0;
		sn.setK(sn, tn);
		sn.pindex = sn.index;
		
		tn.distR = 0;
		tn.setKr(sn, tn);
		tn.pindexR = tn.index;
		
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

						tt.pindex = processing.index;				//min path is my daddy
					}
					
				}
				
				processing.processed = true;
				++processed;
				
				if(processing.processedR == true) {	
					/*
					long tp = processing.dist + processing.distR;
					++dblProcessed;
					if( tp < mu) {
						mu = tp; 
						result = mu;
						cNode = processing.index;
					}
					if(!heap.isEmpty()) {
						if(heap.get(0).k > mu && heapR.get(0).kr > mu) {
							System.out.println("Exiting on break");
							break;
						}
					}
					 */
					result = shortestPath();
					break;
										
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
					
						ttr.pindexR = processingR.index;
					}
					
				}
				processingR.processedR = true;
				++processed;
				if(processingR.processed == true) {
					/*
					long tp = processingR.dist + processingR.distR;
					++dblProcessed;
					if( tp < mu) {
						mu = tp;
						result = mu;
					}
					if(!heapR.isEmpty()) {
						if(heap.get(0).k > mu && heapR.get(0).kr > mu) {
							break;
						}
					}
					*/
					result = shortestPath();
					break;
				}
			}
			
		}
//		long finish = System.nanoTime();
//		long elapsed = (finish - start) / 1000000;
		System.out.println("BiAStar double processed nodes: " + dblProcessed);
		return result;
    }
	
	public long dijkstra(int s, int t) {
		
		//this is a checking method using the simple Dijkstra algorithm for testing
		
		
		double prt = 0;									//in simple Dijkstra there is no potential
		long mu = INFINITY;
		processed = 0;
		
		initializeQueues();
		resetWorkingNodes();
		
		cNode = -1;
		
		Node sn = map.get(s);
		Node tn = map.get(t);
		
		sn.dist = 0;
		sn.setK(sn, tn, false);							//false => don't use the A* potential in the K value
	
		enQueue(sn, heap);
		sn.queued = true;
		working.add(sn);								//add the initial nodes to the working set	
		
		if(s == t) {
			return 0;									//I found myself!!
		}
		
		
	
		while(!heap.isEmpty()) {						//process the next node in the forward graph
			
			Node processing = getMin(heap);
			processing.queued = false;
			
//				System.out.println("processing Node: " + processing.index);
			if(processing != tn) {
				for(Edge e : graph.get(processing.index)) {
					
					Node tt = e.v;	
					
					long td = processing.dist + e.length;
				
					if(tt.dist > td) {
							
							working.add(tt);
							
							if(tt.queued == false) {	
								tt.dist = td;
								tt.setK(sn, tn, false);
								enQueue(tt, heap);
								tt.queued = true;
							} else {
								decreaseKey(tt, td, sn, tn, heap, false, false);		//use the unidirectional non potential version of the method
							}

//							tt.pindex = processing.index;				//min path is my daddy
					}
						
				}
				
				processing.processed = true;
				++processed;
				
			} else {							//processing the target node here
				
				if(processing.dist < mu)
					mu = processing.dist;
				processing.processed = true;
			}
		
			
			
			if(tn.processed == true) {				
				if(!heap.isEmpty()) {
					if(heap.get(0).k > mu) {	//stop when the shortest queued node distance is longer than the shortest path found
						break;
					}
				}
			}				
		}
//		System.out.println("Dijkstra processed edges: " + processed);
		return mu == INFINITY? -1: mu;

	}
}



public class DistWithCoords {
    

    public static void main(String args[]) {
    	
		
    	
    	if(args.length != 0) {

    		if(args[0].equals("test")) {
    			int W = 19;						//node grid W x H
        		int H = 15;
        		int L = 32;						//node spacing
        		int C = 4;						//random edge extension factor (L / C)
        		
        		BiGraph g = new BiGraph(W * H);
        		
        		System.out.println("Graph parameters: W: " + W + " H: " + H + " L: " + L + " C: " + C);
        		
        		g.addBarrier(4,4,4,10,L);			//set of inactive nodes
        		g.addBarrier(1,10,9,10,L);
        		g.addBarrier(9,7,9,13,L);
        		g.addBarrier(7,7,17,7,L);
        		g.addBarrier(5,13,15,13,L);			//use (1,13,15,13,L) for unreachable
        		g.addBarrier(10,3,13,6,L);
        		
        		for(Barrier b : g.barriers) {
        			System.out.println("barrier: (" + b.ll.x + ", " + b.ll.y + ", " + b.ur.x + ", " + b.ur.y + ")" );
        		}
        		
        		long start = System.nanoTime();
        		g.genTestGraph(W, H, L, C);
        		long finish = System.nanoTime();
        		long elapsed = (finish - start)/1000000;
        		
        		System.out.println("Graph generation time: " + elapsed + " ms" );
        		
        		start = System.nanoTime();
        		
        		int strt = Node.nodeNumber(3, 2, W);
        		int trgt = Node.nodeNumber(7, 12, W);
        		
        		System.out.println("biAStar start: " + strt + " target: " + trgt + " Distance: " + g.biAStar(strt, trgt));
        		
        		finish = System.nanoTime();
        		elapsed = (finish - start) / 1000000;
        		System.out.println("biAStar time: " + elapsed + " ms" );
        		
        		
        		start = System.nanoTime();
        		
        		
        		System.out.println("dijkstra start: " + strt + " target: " + trgt + " Distance: " + g.dijkstra(strt, trgt));
        		
        		finish = System.nanoTime();
        		elapsed = (finish - start) / 1000000;
        		System.out.println("dijkstra time: " + elapsed + " ms" );    		
        		
        		

    		} else {
    			
    			//to run a test file invoke the class with a path parameter ie: java DistWithCoords distTests/01
    			
    			String p = args[0];
    			String pp = p + ".a";
    			
    			Path fp = Paths.get(p);
    			Path ep = Paths.get(pp);
    			Charset cset = Charset.forName("US-ASCII");
    			String s;
    			
    			
    			try(BufferedReader br = Files.newBufferedReader(fp, cset)){
    				s = br.readLine();
    				String[] params = s.split(" ");
    				int n = Integer.valueOf(params[0]);
    				int m = Integer.valueOf(params[1]);
    				
    				BiGraph g = new BiGraph(n + 1);
    				
    				System.out.println("Running file: " + fp.toString());
    				System.out.println("Points: " + n + " Edges: " + m);
    				
    				int points = 0;
    				
    				
    				for(int i = 1; i <= n; ++i){
    					s = br.readLine();
    					params = s.split(" ");
    					
    					Point pnt = new Point(Integer.valueOf(params[0]) , Integer.valueOf(params[1]));
    		            g.addNode(i, pnt);
    		            ++points;
    					
    				}
    				
    				System.out.println("Generated points: " + points);
    				
    				int edges = 0;
    				
    				for(int j = 1; j <= m; ++j) {
    					s = br.readLine();
    					params = s.split(" ");
    					
    					g.addEdges(Integer.valueOf(params[0]), Integer.valueOf(params[1]), Integer.valueOf(params[2]));
    					++edges;
    				}
    				
    				System.out.println("Generated edges: " + edges);
    				
    				s = br.readLine();
    				
    				int tests = Integer.valueOf(s);
    				
    				System.out.println("Running tests: " + tests);
    				
    				String[] expected = new String[tests];
    				
    				
    				try(BufferedReader bre = Files.newBufferedReader(ep, cset)){
    					for(int jj = 0; jj < tests; ++jj) {
    						expected[jj] = bre.readLine();
    					}
    				} catch (IOException e) {
    					e.printStackTrace();
    				}
    				
    				int biAStarSum = 0;
    				int dijkstraSum = 0;
    				
    				for(int k = 0; k < tests; ++k) {
    					s = br.readLine();
    					params = s.split(" ");
    					
    					int start = Integer.valueOf(params[0]);
    					int target = Integer.valueOf(params[1]);
    					
    					long bistarDist = g.biAStar(start, target);
    					int bistarNodes = g.processed;
    					biAStarSum += bistarNodes;
    					long dijkstraDist = g.dijkstra(start, target);
    					int dijkstraNodes = g.processed;
    					dijkstraSum += dijkstraNodes;
    					long expectedDist = Long.valueOf(expected[k]);
    					
    					System.out.println("Test: " + k + " biAStar: " + bistarDist + " nodes: "+ bistarNodes + " Dijkstra: " + dijkstraDist + " nodes: "+ dijkstraNodes +" Expected: " + expectedDist );
    					if(bistarDist != expectedDist)
    						System.out.println("NOT EXPECTED!!");
    					if(bistarDist != dijkstraDist)
    						System.out.println("NOT VERIFIED!!");
    					
    				}
    				
    				System.out.println("BiAStar Sum: " + biAStarSum + " Dijkstra Sum: " + dijkstraSum);
    				
    				
    				
    			} catch (IOException e) {
    				e.printStackTrace();
    			}
    			

    		}
   		
    	} else {		//no initial parameter passed
    	
    	
    	
	        Scanner in = new Scanner(System.in);
	        int n = in.nextInt();
	        int m = in.nextInt();
	        
	        BiGraph g = new BiGraph(n + 1);
	
	        for (int i = 1; i <= n; i++) { 
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
	            System.out.println(g.biAStar(u, v));
	        }
	        in.close();
	    }
    }
}
