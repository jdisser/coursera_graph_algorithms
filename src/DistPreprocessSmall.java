import java.util.Scanner;
import java.util.function.BiPredicate;
import java.util.function.Consumer;
import java.util.function.Function;
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
	public boolean processedP;
	public boolean queued;			//in heap
	public boolean queuedR;			//in heapR
	public boolean queuedP;			//in preProc
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
	public long shortcutDist;		//the minimum distance u->v->w for a shortcut to be created
	public int shortcuts;
	
	
	
	public Node(int i) {
        this.index = i;
        this.dist = INFINITY;
        this.distR = INFINITY;
        this.processed = false;
        this.processedR = false;
        this.processedP = false;
        this.queued = false;
        this.queuedR = false;
        this.queuedP = false;
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
        this.shortcutDist = INFINITY;
        this.shortcuts = 0;
	}
	
	
	public static int nodeNumber(int i, int j, int w) {	//1 index position of node in test graph w x h
		return w*(j - 1) + i;
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
	public PriorityNodeQ heap;						//priority queue with method to update key values
	public PriorityNodeQ heapR;						//priority queue for reversed graph
	public PriorityNodeQ preProc;					//priority queue for graph preprocessing
	public HashMap<Integer,Node> map;				//returns Nodes by vertex (index)
	
	
	public final long INFINITY = Long.MAX_VALUE / 4;

	public int n;
	public long mu = INFINITY;
	public int processed;
	
	public int nRank;
	public int maxHop;
	
	TableHash hTable;

	
	public BiGraph(int n, TableHash ht) {
		this.n = n;
		this.hTable = ht;			//inject a hash table object
		this.nRank = 0;
		this.maxHop = 5;			//adjust this for best performance

        this.graph = new ArrayList<ArrayList<Edge>>();
        this.graphR = new ArrayList<ArrayList<Edge>>();
        
        Consumer<Node> resetHeap = (nd) -> {
        	nd.dist = INFINITY;
            nd.processed = false;
            nd.queued = false;
            nd.parent = null;
            nd.hops = Integer.MAX_VALUE;		//this field is used in the shortcut method to limit the dijkstra search
            nd.shortcutDist = INFINITY;			//used in shortcut to determine witness paths
            nd.active = true;					//used to remove the node during the shortcut search
        };
        
        Function<Node, Long> getHeapKey = (nd) -> { return nd.dist; };
        
        Consumer<Node> resetHeapR = (nd) -> {
        	nd.distR = INFINITY;
            nd.processedR = false;
            nd.queuedR = false;
            nd.parentR = null;
        };
        
        Function<Node, Long> getHeapRKey = (nd) -> { return nd.distR; };
        
        Consumer<Node> resetPreProc = (nd) -> {

            nd.queuedP = false;
            nd.processedP = false;
            nd.priority = 0;
            //these Node properties are persisted after buildPriorityQueue and shortcut(true) contraction
//            this.level = 0;
//            this.neighbors = 0;
//            this.rank = 0;
//            this.edgeDiff = 0;
        };
        
        Function<Node, Integer> getPreProcKey = (nd) -> { return nd.priority; };
        
        this.heap = new PriorityNodeQ(ht, (a , b) -> a.dist < b.dist, (a, b) -> a.dist == b.dist, resetHeap, getHeapKey, "heap", false);
        this.heapR = new PriorityNodeQ(ht, (a , b) -> a.distR < b.distR, (a, b) -> a.distR == b.distR, resetHeapR, getHeapRKey, "heapR", false);
        this.preProc = new PriorityNodeQ(ht, (a , b) -> a.priority < b.priority, (a, b) -> a.priority == b.priority, resetPreProc, getPreProcKey, "preProc", true);
        
        
        
        this.map = new HashMap<Integer,Node>(n+1);
        
        
        
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
	
	/*
	public void initializeQueues() {
		
		heap.initializeQueue();
		heapR.initializeQueue();
		preProc.initializeQueue();
		
	}
	*/
	
	public void addEdges(int s, int t, int c) {		//(source, target, length in 1 based indexing)

		
		Node source = map.get(s);
		Node target = map.get(t);
		
		
		Edge e = new Edge(source, target, c);
        Edge er = new Edge(target, source, c);

        graph.get(s).add(e);
		graphR.get(t).add(er);

	}
	
	//TODO: determine if this should be in the PriorityNodeQ Class
	public Node minPriority(Node a, Node b) {
		//returns minimum priority node with hash tiebreaker
		
		Node min = null;
		
		if(a == null || b==null)
			return min;
		
		if(a.priority < b.priority)
			min = a;
		else if(b.priority < a.priority)
			min = b;
		else if(hTable.hash(a.index) < hTable.hash(b.index))
			min = a;
		else
			min = b;
			
		
		return min;
		
	}

	//TODO: need a combined working set!!
	/*
    private long shortestPath() {
    	long d = INFINITY;
    	long dn = 0;
    	for(Node n: working) {
    		dn = n.dist + n.distR;
    		d = Math.min(d, dn);
    	}
    	return d;					
    }
	*/
	
    //TODO  This method is to be rewritten as a bidirectional dijkstra with rank
    /*
     *
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
					
//					long tp = processing.dist + processing.distR;
//					++dblProcessed;
//					if( tp < mu) {
//						mu = tp; 
//						result = mu;
//						cNode = processing.index;
//					}
//					if(!heap.isEmpty()) {
//						if(heap.get(0).k > mu && heapR.get(0).kr > mu) {
//							System.out.println("Exiting on break");
//							break;
//						}
//					}
					
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
					
//					long tp = processingR.dist + processingR.distR;
//					++dblProcessed;
//					if( tp < mu) {
//						mu = tp;
//						result = mu;
//					}
//					if(!heapR.isEmpty()) {
//						if(heap.get(0).k > mu && heapR.get(0).kr > mu) {
//							break;
//						}
//					}
					
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
	*/
    
    
	public long dijkstra(int s, int t) {
		
		long mu = INFINITY;
		processed = 0;
		
		heap.initializeQueue();
		heap.resetWorkingNodes();

		
		Node sn = map.get(s);
		Node tn = map.get(t);
		
		sn.dist = 0;
								
	
		heap.enQueue(sn);
		sn.queued = true;
		
		
		if(s == t) {
			return 0;									//I found myself!!
		}
		
		
	
		while(!heap.isEmpty()) {						//process the next node in the forward graph
			
			Node processing = heap.getMin();
			processing.queued = false;
			
//				System.out.println("processing Node: " + processing.index);
			if(processing != tn) {
				for(Edge e : graph.get(processing.index)) {
					
					Node tt = e.v;	
					
					long td = processing.dist + e.length;
				
					if(tt.dist > td) {
							
							if(tt.queued == false) {	
								tt.dist = td;
								heap.enQueue(tt);
								tt.queued = true;
							} else {
								tt.dist = td;
								heap.decreaseKey(tt);
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
					if(heap.peek().dist > mu) {	//stop when the shortest queued node distance is longer than the shortest path found
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
		
		System.out.println("Shortcut contract: " + contract);
		
		int shortcuts = 0;
		
		heap.initializeQueue();		
		
		ArrayList<Edge> inEdges = graphR.get(v.index);
		ArrayList<Edge> outEdges = graph.get(v.index);
		
		int ins = inEdges.size();
		int outs = outEdges.size();
		
		
		ArrayList<Node> us = new ArrayList<Node>();
		ArrayList<Node> ws = new ArrayList<Node>();

		
		for(Edge e : inEdges) {			//get source nodes
			if(!contract)
				us.add(e.u);
			else
				if(!e.u.contracted)		//while contracting the graph ignore previously contracted nodes
					us.add(e.u);
		}
		
		for(Edge e : outEdges) {		//get target nodes
			if(!contract)
				ws.add(e.v);
			else
				if(!e.v.contracted)
					ws.add(e.v);
		}
		
		
		
		v.active = false;				//remove v from the active graph
		
		for(Node u : us) {				//check for a witness path to each target
			
			heap.resetWorkingNodes();	
			
			long mu = 0;
			processed = 0;
			heap.initializeQueue();
			
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

			heap.enQueue(u);
			u.queued = true;
			
			while(!heap.isEmpty()) {						//process the next node in the forward graph
				
				Node x = heap.getMin();
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
									heap.working.add(tt);
									tt.dist = td;
									tt.parent = x;
									tt.hops = x.hops + 1;
								} else {
									if(tt.queued == false) {
										tt.hops = x.hops + 1;
										tt.dist = td;
										heap.enQueue(tt);
										tt.queued = true;
									} else {
										tt.dist = td;
										heap.decreaseKey(tt);		
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
						Edge scr = new Edge(w, v, u, w.shortcutDist);
						graphR.get(w.index).add(scr);
					}
					
				}
			}
		}

		heap.working.add(v);
		heap.resetWorkingNodes();							
		v.shortcuts = shortcuts;
		v.edgeDiff = shortcuts - ins - outs;				//set the edge difference NOTE: call setPriority on a node to update the priority property
		if(contract) {
			v.contracted = true;
			v.rank = ++nRank;
			for(Node w : ws) {								//increase contracted neighbors
				++w.neighbors;
			}
			for(Node u : us) {
				++u.neighbors;
				u.level = Math.max(u.level, v.level + 1); 	//update level of neighbors
			}
		}
		
		return v.edgeDiff;									//return the edge difference
	}
	
	public void buildPriorityQueue() {
		
		int priority = 0;
		
		for(int i = 1; i < n; ++i) {
			System.out.println();
			System.out.print("Processing Node: " + i);
			Node n = map.get(i);
			priority = shortcut(n, false, maxHop);
			System.out.println(" Priority: " + priority);
			n.setPriority();
			preProc.enQueue(n);
			n.queuedP = true;
		}
	}
	
	public void contractGraph() {
		int loops = 0;
		while(!preProc.isEmpty()) {
			++loops;
			if(loops > 100)
				break;
			Node n = preProc.getMin();
			System.out.println("Pop: " + n.index + " priority: " + n.priority);
			shortcut(n, false, maxHop);
			n.setPriority();
			System.out.println("Reprioritized: " + n.index + " priority: " + n.priority);
			Node min = preProc.peek();
			System.out.println("preProc min: " + min.index + " priority: " + min.priority);
			if(n != minPriority(n, min)) {
				System.out.println("enQueue: " + n.index + " priority: " + n.priority);
				preProc.enQueue(n);
				continue;
			} else {
				shortcut(n, true, maxHop);
				n.setPriority();
				n.processedP = true;
				System.out.println("contracted: " + n.index + " priority: " + n.priority);
			}
			++loops;
			if(loops > 100)
				break;
				
		}
		preProc.resetWorkingNodes();
	}

}

class TableHash{
	
	ArrayList<ArrayList<Integer>> hashTable;
	ArrayList<Integer> masks;
	
	int bits;
	int r;
	int k;
	int n;
	
	
	public TableHash(int bits, int r) {
		
		this.bits = bits;
		this.r = r;
		
		this.hashTable = new ArrayList<ArrayList<Integer>>();
		for(int i = 0; i < r; ++i) {
			this.hashTable.add(new ArrayList<Integer>());
		}
		
		int m = (1 << (bits + 1)) - 1;		//same as 2 ^ (bits + 1) - 1 in ansi c
		this.k = bits/r;
		this.n = (1 << (k + 1)) - 1;
		
//		System.out.println("Bits: " + bits + " r: " + r + " m: " + m + " k: "+ k + " n: " + n);
		
		
		for(int i = 0; i < r; ++i) {
			
			ArrayList<Integer> l = hashTable.get(i);
			
			for(int j = 0; j < n; ++j) {
				l.add(Long.valueOf((Math.round(Math.random() * m))).intValue());
			}
		}
		
		/*
		System.out.println("HashTable:");
		for(int i = 0; i < r; ++i) {
			for(int j = 0; j < n; ++j) {
				System.out.println("i: " + i + " j: " + j + " htij: " + hashTable.get(i).get(j));
			}
		}
		*/
		
		masks = new ArrayList<Integer>();
		
		for(int i = 0; i < r; ++i) {
			int msk = ((1 << k) - 1) << k * i;
			masks.add(i, msk);
		}
		
		
	}
	
	public int hash(int x) {
		
		int result = 0;
		int hti = 0;
		
		for(int i = 0; i < r; ++i) {
			hti = (x & masks.get(i)) >> k * i;
//			System.out.println("hti: " + hti);
			result ^= hashTable.get(i).get(hti);
		}
		
		return result;
		
		
	}
	
}

class PriorityNodeQ {
	
	/*
	 * This is a general priority queue for graph nodes that can use different node properties as keys
	 * Two lambda expressions are passed in the constructor for the predicate functions a < b and a == b
	 * where a & b are Nodes. The Queue is named for diagnostic printing. If tieBreak is passed as true in the
	 * constructor a hash function is used to order nodes with equal keys.
	 */
	
	
	private ArrayList<Node> heap;
	private TableHash hTable;
	public HashSet<Node> working;					//all nodes processed in this queue that must be reset
	private String name;
	
	//these are passed in as a Lambda functions to allow for different Node properties as keys
	private BiPredicate<Node, Node> minNode;		//test(a,b) -> true if a < b	
	private BiPredicate<Node, Node> equNode;		//test(a,b) -> true if a = b

	
	private Consumer<Node> resetNode;				//resets the properties of nodes used in this queue
	private Function<Node, ? extends Number> getKey;
	
	private boolean tieBreak;						//if tieBreak = true use the hash of the node index as a tiebreaker in sifts
	
	public PriorityNodeQ(
			TableHash hTable,
			BiPredicate<Node,Node> minNode,
			BiPredicate<Node,Node> equNode,
			Consumer<Node> resetNode,
			Function<Node, ? extends Number> getKey,
			String name,
			boolean tieBreak ) {
		
		this.hTable = hTable;
		this.minNode = minNode;
		this.equNode = equNode;
		this.resetNode = resetNode;
		this.getKey = getKey;
		this.name = name;
		this.tieBreak = tieBreak;
	}
	
	public void decreaseKey(Node dn) {
	
		//this method assumes the property used as a key of the Node has been altered prior to the method call
		int dni = heap.indexOf(dn);
		minSiftUp(dni);
	}
	
	public void enQueue(Node n) {
		
		working.add(n);				
		int end = heap.size();
		heap.add(end, n);
		minSiftUp(end);

	}
	
	public Node getMin() {
		
		Node nm = heap.get(0);
		swap(0, heap.size() - 1);
		heap.remove(heap.size() - 1);
		minSiftDown(0);
		
		return nm;
	}
	
	public void initializeQueue() {
		
		heap.clear();

	}
	
	private void minSiftDown(int i) {
		
//		System.out.println("minSiftDown 1: " + i);
//		printHeap();
		
    	int mini = i;
    	int li = heap.size() - 1;	//0 based last index
    	int left = 2*i + 1;
    	int rght = 2*i + 2;
    	int lhi = 0;
    	int rhi = 0;
    	int minhi;
    	
    	if(!tieBreak) {
    		if(left <= li && minNode.test(heap.get(left), heap.get(mini))) {	//use the minNode lambda for the test
        		mini = left;
        	}
        	if(rght <= li && minNode.test(heap.get(rght), heap.get(mini))) {
        		mini = rght;
        	}
        	if(mini != i) {
    	   		swap(i, mini);
    	   		minSiftDown(mini);
        	}
    	}  else {
    		
    		if(left <= li)
    			lhi = hTable.hash(heap.get(left).index);					//use the hash of the index to break ties
    		if(rght <= li)
    			rhi = hTable.hash(heap.get(rght).index);
    		minhi = hTable.hash(heap.get(mini).index);
    		
    		//TODO: there is a small probability that the hash will not be unique ~ 4/3000 handle this case here
    		
    		if(left <= li && minNode.test(heap.get(left), heap.get(mini))) {
        		mini = left;
        	} else {
        		if(left <= li && equNode.test(heap.get(left), heap.get(mini))) {
        			if(lhi < minhi)
        				mini = left;
        		}
        	}
        	if(rght <= li && minNode.test(heap.get(rght), heap.get(mini))) {
        		mini = rght;
        	} else {
        		if(rght <= li && equNode.test(heap.get(rght), heap.get(mini))) {
        			if(rhi < mini)
        				mini = rght;
        		}
        	}
        	if(mini != i) {
    	   		swap(i, mini);
    	   		minSiftDown(mini);
        	}
    	}
 
    }
	
	private void minSiftUp(int i) {
		
		int pi;
		int hi;
		int hpi;
		
		System.out.print(" minSiftUp i: " + i);
//		printHeap();
		System.out.println();
		
		if(i > 0) {
			pi = (i - 1)/2;
			
			
			if(!tieBreak) {
				
				if(minNode.test(heap.get(pi), heap.get(i))) {
					swap(pi, i);
					minSiftUp(pi);
				}
			} else {
				
				hi = hTable.hash(heap.get(i).index);					//use the hash of the index to break ties
	    		hpi = hTable.hash(heap.get(pi).index);
	    		
	    		//TODO: there is a small probability that the hash will not be unique ~ 4/3000 handle this case here
	    		
	    		if(minNode.test(heap.get(pi), heap.get(i))) {
					swap(pi, i);
					minSiftUp(pi);
				} else {
					if(equNode.test(heap.get(pi), heap.get(i))) {
						if(hpi > hi) {
							swap(pi, i);
							minSiftUp(pi);
						}
					}
				}	
			}
			
			
		}

	}

	public Node peek() {
		return heap.get(0);
	}
	
	public void printHeap() {
		System.out.println();
		
		
		System.out.println(name + ": ");
		for(Node n : heap) {
			System.out.print(heap.indexOf(n) + "|" + n.index + ":" + getKey.apply(n) + ", ");
		}
		System.out.println();
	}
	
	
	private void swap(int i, int n) {
    	
    	Node tn = heap.get(i);
    	heap.set(i, heap.get(n));
    	heap.set(n, tn);
    	
    }
	
	public void resetWorkingNodes() {
		
		if(!working.isEmpty()) {
			for(Node n: working) {
				resetNode.accept(n);;
			}
			working.clear();
		}
		
	}
	
	public boolean isEmpty() {
		return heap.isEmpty();
	}
	
	
	
	
}



class DistPreprocessSmall {
    
    public static void main(String args[]) {
    	
		TableHash tHash = new TableHash(18,3);
    	
    	if(args.length != 0) {

    		
    			//TODO: rewrite this section to process the test graphs without the euclidian coordinates
    			//to run a test file invoke the class with a path parameter ie: java DistPreprocessSmall distTests/02P
    			
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
    				
    				
    				
    				BiGraph g = new BiGraph(n + 1, tHash);
    				
    				System.out.println("Running file: " + fp.toString());
    				System.out.println("Nodes: " + n + " Edges: " + m);
    				
    				int points = 0;
    				
    				
    				for(int i = 1; i <= n; ++i){	//Nodes 1..n
//    					s = br.readLine();
//    					params = s.split(" ");
	
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

    				HashSet<Integer> test = new HashSet<Integer>();
    				
    				System.out.println("Running Hash Table collision test...");
    				for(int i = 1; i < 3000; ++i) {
    					int x = tHash.hash(i);
    					if(test.contains(x))
    						System.out.println("Collision: " + i + " Hash: " + x);
    					test.add(x);
    				}
    				
    				System.out.println("Building Priority Queue: ");
    				g.buildPriorityQueue();
    				System.out.println("Contracting: ");
    				g.contractGraph();
    				System.out.println("Ready");
    				
    				
    				
    				
    				/*
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
    				
    				//TODO: rewrite this section to use the biDijkstraCh() method instead of the biStar
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
    				*/
    				
    				
    			} catch (IOException e) {
    				e.printStackTrace();
    			}
    			

    		
   		
    	} else {		//no initial parameter passed
    	
    	
    	
	        Scanner in = new Scanner(System.in);
	        int n = in.nextInt();
	        int m = in.nextInt();
	        
	        BiGraph g = new BiGraph(n + 1, tHash);
	
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
	
	        
//	        g.preprocess();
	        System.out.println("Ready");
	        
	        int t = in.nextInt();
	
	        for (int i = 0; i < t; i++) {
	            int u, v;
	            u = in.nextInt();
	            v = in.nextInt();
//	            System.out.println(g.biDijkstraCh(u, v));
	        }
	        in.close();
	    }
    }
}
