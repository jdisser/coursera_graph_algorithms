import java.util.Scanner;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

class Node {
	private static final long INFINITY = Long.MAX_VALUE/4;
	public int index;				//this is the map index for back reference
	public long dist;				//this is the distance from the source in the graph
	public long distR;				//this is the distance from the target in the reverse graph
	public boolean processed;
	public boolean processedR;
	public boolean queued;
	public boolean queuedR;
	public boolean active;			//graph state false => not available in graph
	public boolean contracted;		//contracted state true => contracted
	public Node parent;
	public Node parentR;
//	public HashSet<Edge> bList;		//bucket list for single hop backward Edges search during contraction
	public int level;				//contraction hueristic
	public int priority;			//priority for contraction
	public int neighbors;			//number of contracted neighbors
	public int rank;				//order of node in contracted graph
	public int edgeDiff;			//shortcuts - inDegree - outDegree
	public int hops;				//number of edges from source used to terminate contraction
	public long key;				//The key used in the priority queue
	public long keyR;				//the reverse search key
	public long shortcutDist;		//the minimum distance u->v->w for a shortcut to be created
	public int shortcuts;
	
	
	
	public Node(int i) {
        this.index = i;
        this.dist = INFINITY;
        this.distR = INFINITY;
        this.processed = false;
        this.processedR = false;
        this.queued = false;
        this.queuedR = false;
        this.contracted = false;
        this.parent = null;
        this.parentR  = null;
//        this.bList = new HashSet<Edge>();
        this.active = true;
        this.level = 0;
        this.priority = 0;
        this.neighbors = 0;
        this.rank = 0;
        this.edgeDiff = 0;
        this.hops = Integer.MAX_VALUE;
        this.key = 0;
        this.keyR = 0;
        this.shortcutDist = INFINITY;
        this.shortcuts = 0;
	}
	
	
	//Note: this method does NOT reset the shortcuts, contracted, level
	public void resetNode(boolean reactivate) {
		this.dist = INFINITY;
        this.distR = INFINITY;
        this.processed = false;
        this.processedR = false;
        this.queued = false;
        this.queuedR = false;
        this.parent = null;
        this.parentR = null;
        this.key = 0;
        this.keyR = 0;
        //can blist be eliminated and the graphR list be used instead??? No but blist not implemented
//        this.bList.clear();
        if(reactivate)
        	this.active = true;
//        this.level = 0;
        //TODO: determine how to handle the heuristic properties for resets during prioritizing, contracting and searching
        this.priority = 0;
        this.neighbors = 0;
        this.rank = 0;
        this.edgeDiff = 0;
        this.hops = Integer.MAX_VALUE;
        this.shortcutDist = INFINITY;
	}

	
	public static int nodeNumber(int i, int j, int w) {	//1 index position of node in test graph w x h
		return w*(j - 1) + i;
	}
	
	
	public long setK() {
		key = dist ;
		return key;
	}
	
	public long setK(boolean prioritize) {
		if(prioritize)
			key = priority;
		else
			key = dist;
		return key;
	}
	
	public long setKr() {
		keyR = distR;					
		return keyR;
	}
	
	public long setKr(boolean prioritize) {
		
		if(prioritize)
			keyR = priority;					
		else
			keyR = distR;
		
		return keyR;
	}
	
	public void setPriority() {
		this.priority = this.edgeDiff + this.neighbors + this.shortcuts + this.level;
	}
	
}



class Edge {
	public int target;
	public int source;
	public long length;
	public Node u;
	public Node v;
	public Node cn;				//shortcut u--cn-->v
	public boolean shortcut;
	
	
	public Edge(Node u, Node v, int l, boolean shortcut) {
		this.target = u.index;
		this.source = v.index;
		this.u = u;
		this.v = v;
		this.length = l;
		this.shortcut = shortcut;
	}
	
	//shortcut constructor u-->v-->w => u --> v(was w) & cn <- v
	public Edge(Node u, Node v, Node w, long l) {
		
		this.u = u;
		this.v = w;
		this.cn = v;
		this.length = l;
	}
	
	//Edge constructor u --> v && l(u,v) = l
	public Edge(Node u, Node v, long l) {
		this.u = u;
		this.v = v;
		this.cn = null;
		this.length = l;
	}


	
	
	
}

class BiGraph {
	
	public ArrayList<ArrayList<Edge>> graph;		//adjacency list of edges
	public ArrayList<ArrayList<Edge>> graphR;		//adjacency list of reversed edges
	public ArrayList<Node> heap;					//priority queue with method to update key values
	public ArrayList<Node> heapR;					//priority queue for reversed graph
	public ArrayList<Node> preProc;					//priority queue for graph preprocessing
	public HashMap<Integer,Node> map;				//returns Nodes by vertex (index)
	public HashSet<Node> working;					//all nodes processed by query that must be reset
	
	public final long INFINITY = Long.MAX_VALUE / 4;

	public int n;
	public long mu = INFINITY;
	public int processed;
	
	TableHash hTable;

	
	public BiGraph(int n, TableHash ht) {
		this.n = n;
		this.hTable = ht;			//inject a hash table object

        this.graph = new ArrayList<ArrayList<Edge>>();
        this.graphR = new ArrayList<ArrayList<Edge>>();
        
        this.heap = new ArrayList<Node>();
        this.heapR = new ArrayList<Node>();
        
        this.preProc = new ArrayList<Node>();
        
        
        
        this.map = new HashMap<Integer,Node>(n+1);
        
        this.working = new HashSet<Node>();
        
        for (int i = 0; i <= n; i++) {					//graph indexes and nodes are 1 indexed, index 0 not used
        	
            this.graph.add(new ArrayList<Edge>());
            this.graphR.add(new ArrayList<Edge>());          
  
        }

        this.processed = 0;
  
	}
	
	
	
	
	
	
	
	public Node addNode(int i) {
		Node n = new Node(i);
		map.put(i, n);
		return n;
	}
	
	public void initializeQueues() {
		
		heap.clear();
		heapR.clear();
		
	}
	
	public void resetWorkingNodes(boolean reactivate) {
		if(!working.isEmpty()) {
			for(Node n: working) {
				n.resetNode(reactivate);
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
    	int lhi;
    	int rhi;
    	int minhi;
    	
    	if(h == heap) {
    		if(left <= li && h.get(left).key < h.get(mini).key) {
        		mini = left;
        	}
        	if(rght <= li && h.get(rght).key < h.get(mini).key) {
        		mini = rght;
        	}
        	if(mini != i) {
    	   		swap(i, mini, h);
    	   		minSiftDown(mini, h);
        	}
    	}
    	if(h == heapR){
    		if(left <= li && h.get(left).keyR < h.get(mini).keyR) {
        		mini = left;
        	}
        	if(rght <= li && h.get(rght).keyR < h.get(mini).keyR) {
        		mini = rght;
        	}
        	if(mini != i) {
    	   		swap(i, mini, h);
    	   		minSiftDown(mini, h);
        	}
    	}
    	if(h == preProc) {
    		
    		lhi = hTable.hash(h.get(left).index);					//use the hash of the index to break ties
    		rhi = hTable.hash(h.get(rght).index);
    		minhi = hTable.hash(h.get(mini).index);
    		
    		if(left <= li && h.get(left).key < h.get(mini).key) {
        		mini = left;
        	} else {
        		if(left <= li && h.get(left).key == h.get(mini).key) {
        			if(lhi < minhi)
        				mini = left;
        		}
        	}
        	if(rght <= li && h.get(rght).key < h.get(mini).key) {
        		mini = rght;
        	} else {
        		if(rght <= li && h.get(rght).key == h.get(mini).key) {
        			if(rhi < mini)
        				mini = rght;
        		}
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
		int hi;
		int hpi;
//		System.out.println(" minSiftUp i: " + i);
//		printHeap(h);
		if(i > 0) {
			pi = (i - 1)/2;
			
//			System.out.println(" minSiftUp pi: " + pi);
			if(h == heap) {
				if(h.get(pi).key > h.get(i).key) {
					swap(pi, i, h);
					minSiftUp(pi, h);
				}
			} 
			if(h == heapR) {
				if(h.get(pi).keyR > h.get(i).keyR) {
					swap(pi, i, h);
					minSiftUp(pi, h);
				}
			}
			if(h == preProc) {
				hi = hTable.hash(h.get(i).index);					//use the hash of the index to break ties
	    		hpi = hTable.hash(h.get(pi).index);
	    		
	    		if(h.get(pi).key > h.get(i).key) {
					swap(pi, i, h);
					minSiftUp(pi, h);
				} else {
					if(h.get(pi).key == h.get(i).key) {
						if(hpi > hi) {
							swap(pi, i, h);
							minSiftUp(pi, h);
						}
					}
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

	private void decreaseKey(Node dn, long d, ArrayList<Node> h) {
		
//		System.out.println(" decreaseKey i: " + i + " d: " + d);
		
		if(h == heapR) {
			dn.distR = d;
			dn.setKr();	//no reverse search in priority queue
		} 
		if(h == heap) {
			dn.dist = d;
			dn.setK();
		}
		if(h == preProc) {
			dn.setPriority();
			dn.setK(true);
		}		
		
		int dni = h.indexOf(dn);
		minSiftUp(dni, h);
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
		resetWorkingNodes(true);

		
		Node sn = map.get(s);
		Node tn = map.get(t);
		
		sn.dist = 0;
		sn.setK();
		sn.parent = sn;
		
		tn.distR = 0;
		tn.setKr();
		tn.parentR = tn;
		
		prt = tn.keyR;									//reverse potential of target since tn.distr = 0;
		
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
								tt.setK();
								enQueue(tt, heap);
								tt.queued = true;
							} else {
								decreaseKey(tt, td, heap);
							}

						tt.parent = processing;				//min path is my daddy
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
								ttr.setKr();
								enQueue(ttr, heapR);
								ttr.queuedR = true;
							} else {
								decreaseKey(ttr, tdr, heapR);
							}
					
						ttr.parentR = processingR;
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
//		System.out.println("BiAStar double processed nodes: " + dblProcessed);
		return result;
    }
	
	public long dijkstra(int s, int t) {
		
		//this is a checking method using the simple Dijkstra algorithm for testing
		
		

		long mu = INFINITY;
		processed = 0;
		
		initializeQueues();
		resetWorkingNodes(true);

		
		Node sn = map.get(s);
		Node tn = map.get(t);
		
		sn.dist = 0;
		sn.setK();							
	
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
								tt.setK();
								enQueue(tt, heap);
								tt.queued = true;
							} else {
								decreaseKey(tt, td, heap);		//use the unidirectional non potential version of the method
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
					if(heap.get(0).key > mu) {	//stop when the shortest queued node distance is longer than the shortest path found
						break;
					}
				}
			}				
		}
//		System.out.println("Dijkstra processed edges: " + processed);
		return mu == INFINITY? -1: mu;

	}
	
	public int shortcut(Node v, boolean contract, int hops) {
		
		//if create, create shortcuts else return the edge difference shortcuts - ins - outs
		
		int shortcuts = 0;
		
		ArrayList<Edge> inEdges = graphR.get(v.index);
		ArrayList<Edge> outEdges = graph.get(v.index);
		
		int ins = inEdges.size();
		int outs = outEdges.size();
		
		
		ArrayList<Node> us = new ArrayList<Node>();
		ArrayList<Node> ws = new ArrayList<Node>();

		
		for(Edge e : inEdges) {			//get source nodes
			us.add(e.u);
		}
		
		for(Edge e : outEdges) {		//get target nodes
			ws.add(e.v);
		}
		
		
		
		v.active = false;				//remove v from the active graph
		
		for(Node u : us) {				//check for a witness path to each target
			
			resetWorkingNodes(true);	//TODO: verify that this method can be called for each source
			
			long mu = 0;
			processed = 0;
			initializeQueues();
			
			long uvDist = 0;
			u.dist = 0;
			u.hops = 0;
			
			for(Edge e : inEdges) {						//get d(u,v) for this source node
				if(e.u == u) {
					uvDist = e.length;
				}
			}
			
			long maxShortcut = 0;
			
			for(Edge e : outEdges) {					//for each target set the maximum shortcut distances d(u,v) + d(v,w) from this source
				e.v.shortcutDist = uvDist + e.length;
				maxShortcut = Math.max(maxShortcut, e.v.shortcutDist);
				e.v.dist = INFINITY;
			}

			long minRevDist = INFINITY;
			
			for(Node w : ws) {								//use the graphR edge list as a blist!!
				for(Edge e : graphR.get(w.index)) {
					minRevDist = Math.min(minRevDist, e.length);
					if(e.v == u) {
						w.dist = e.length;
						w.parent = u;
						w.hops = 1;
					}
				}
			}
			
			long dijkstraStop = maxShortcut - minRevDist;

			u.setK();							
		
			enQueue(u, heap);
			u.queued = true;
			working.add(u);									//add the initial node to the working set
			
			
		
			while(!heap.isEmpty()) {						//process the next node in the forward graph
				
				Node x = getMin(heap);
				x.queued = false;
				mu = Math.max(mu, x.dist);
				if(mu >= dijkstraStop || x.hops > hops)		//stopping conditions are max d(u,w) > max (d(u,v)+d(v,w)) - min(x,w) || hops > hops parameter 
					break;
				
//					System.out.println("processing Node: " + processing.index);

				for(Edge e : graph.get(x.index)) {
					
					Node tt = e.v;	
					
					if(tt.active) {										//ignore contracted nodes
						long td = x.dist + e.length;
						
						if(tt.dist > td) {
								if(ws.contains(tt)) {
									//tt is a target node and is never queued
									working.add(tt);
									tt.dist = td;
									tt.parent = x;
									tt.hops = x.hops + 1;
								} else {
									if(tt.queued == false) {
										tt.hops = x.hops + 1;
										tt.dist = td;
										tt.setK();
										enQueue(tt, heap);
										tt.queued = true;
										working.add(tt);
									} else {
										decreaseKey(tt, td, heap);		//use the unidirectional non potential version of the method
										tt.hops = x.hops + 1;
									}
								}
//									tt.pindex = processing.index;				//min path is my daddy
						}
					}
						
				}
					
				++processed;
				x.processed = true;

			}
//			System.out.println("Dijkstra processed edges: " + processed);

			//Count and optionally generate shortcuts
			for(Node w : ws) {
				if(w.dist > w.shortcutDist) {
					//The witness path is longer or non existent
					++shortcuts;
					if(contract) {
						Edge sc = new Edge(u , v , w , w.shortcutDist);
						graph.get(u.index).add(sc);
					}
					
				}
			}
		}

		working.add(v);
		resetWorkingNodes(true);	//reactivate the contracted node!! It remains in the graph!!
		v.shortcuts = shortcuts;
		if(contract) {
			v.contracted = true;
			for(Node w : ws) {		//increase contracted neighbors
				++w.neighbors;
			}
			for(Node u : us) {
				++u.neighbors;
				u.level = Math.max(u.level, v.level + 1);
			}
		}
		v.edgeDiff = shortcuts - ins - outs;
		return v.edgeDiff;			//return the edge difference
	}
}

class TableHash{
	
	ArrayList<ArrayList<Integer>> hashTable;
	ArrayList<Integer> masks;
	
	int bits;
	int r;
	int k;
	
	
	public TableHash(int bits, int r) {
		
		this.bits = bits;
		this.r = r;
		
		this.hashTable = new ArrayList<ArrayList<Integer>>();
		for(int i = 0; i < r; ++i) {
			this.hashTable.add(new ArrayList<Integer>());
		}
		
		int m = 2 ^ (bits + 1) - 1;
		this.k = bits/r;
		
		for(int i = 0; i < r; ++i) {
			
			ArrayList<Integer> l = hashTable.get(i);
			
			for(int j = 0; j < k; ++j) {
				l.add(Long.valueOf((Math.round(Math.random() * m))).intValue());
			}
		}
		
		masks = new ArrayList<Integer>();
		
		for(int i = 0; i < r; ++i) {
			int msk = (2^k - 1) << k * i;
			masks.add(i, msk);
		}
		
		
	}
	
	public int hash(int x) {
		
		int result = 0;
		
		for(int i = 0; i < r; ++i) {
			result ^= hashTable.get(i).get((x & masks.get(i)) >> k * i);
		}
		
		return result;
		
		
	}
	
}

public class DistPreprocessSmall {
    
    public static void main(String args[]) {
    	
		
    	
    	if(args.length != 0) {

    		
    			//TODO: rewrite this section to process the test graphs without the euclidian coordinates
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
    					
    					
    		            g.addNode(i);
    		            ++points;
    					
    				}
    				
    				System.out.println("Generated Nodes: " + points);
    				
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
    				
    				//TODO: rewrite this section to use the biDijkstraCh() method instead of the biStar
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
    			

    		
   		
    	} else {		//no initial parameter passed
    	
    	
    	
	        Scanner in = new Scanner(System.in);
	        int n = in.nextInt();
	        int m = in.nextInt();
	        
	        BiGraph g = new BiGraph(n + 1);
	
	        for (int i = 1; i <= n; i++) { 
	            g.addNode(i);
	        }
	
	        for (int i = 0; i < m; i++) {
	            int x, y, c;
	            x = in.nextInt();
	            y = in.nextInt();
	            c = in.nextInt();
	            
	            g.addEdges(x, y, c);
	            
	 
	        }
	
	        
	        g.preprocess();
	        System.out.println("Ready");
	        
	        int t = in.nextInt();
	
	        for (int i = 0; i < t; i++) {
	            int u, v;
	            u = in.nextInt();
	            v = in.nextInt();
	            System.out.println(g.biDijkstraCh(u, v));
	        }
	        in.close();
	    }
    }
}
