import java.util.Scanner;
import java.util.Set;
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
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Queue;

class DpsNode {
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
	public DpsNode parent;
	public DpsNode parentR;
	public int level;				//contraction hueristic
	public int priority;			//priority for contraction
	public int neighbors;			//number of contracted neighbors
	public int rank;				//order of node in contracted graph
	public int edgeDiff;			//shortcuts - inDegree - outDegree
	public int hops;				//number of edges from source used to terminate contraction
	public long shortcutDist;		//the minimum distance u->v->w for a shortcut to be created
	public int shortcuts;
	public int cover;				//number of nodes targeted by shortcuts
	
	
	
	public DpsNode(int i) {
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
        this.active = true;
        this.level = 0;
        this.priority = 0;
        this.neighbors = 0;
        this.rank = 0;
        this.edgeDiff = 0;
        this.hops = Integer.MAX_VALUE;
        this.shortcutDist = INFINITY;
        this.shortcuts = 0;
        this.cover = 0;
	}

	public void setPriority() {
		this.priority = this.edgeDiff + this.neighbors + this.level + this.cover;
	}
	
}



class DpsEdge {

	public long length;
	public DpsNode u;
	public DpsNode v;
	public DpsNode cn;				//shortcut u--cn-->v cn is the contracted node
	public boolean shortcut;
	
	
	//shortcut constructor u-->v-->w => u --> v(was w) & cn <- v
	public DpsEdge(DpsNode u, DpsNode v, DpsNode w, long l) {
		
		this.u = u;
		this.v = w;
		this.cn = v;
		this.length = l;
		this.shortcut = true;
	}
	
	//Edge constructor u --> v && l(u,v) = l
	public DpsEdge(DpsNode u, DpsNode v, long l) {
		this.u = u;
		this.v = v;
		this.cn = null;
		this.length = l;
		this.shortcut = false;
	}


	
	
	
}

class BiGraph {
	
	public ArrayList<ArrayList<DpsEdge>> graph;		//adjacency list of edges
	public ArrayList<ArrayList<DpsEdge>> graphR;		//adjacency list of reversed edges
	public PriorityNodeQ heap;						//priority queue with method to update key values
	public PriorityNodeQ heapR;						//priority queue for reversed graph
	public PriorityNodeQ preProc;					//priority queue for graph preprocessing
	public HashMap<Integer,DpsNode> map;				//returns Nodes by vertex (index)
	
	
	public final long INFINITY = Long.MAX_VALUE / 4;

	public int n;
	public long mu = INFINITY;

	public int scEdges;
	
	public int nRank;
	public int maxHop;
	
	TableHash hTable;

	
	public BiGraph(int n, TableHash ht) {
		this.n = n;
		this.hTable = ht;			//inject a hash table object
		this.nRank = 0;
		this.maxHop = 5;			//adjust this for best performance

        this.graph = new ArrayList<ArrayList<DpsEdge>>();
        this.graphR = new ArrayList<ArrayList<DpsEdge>>();
        
        Consumer<DpsNode> resetHeap = (nd) -> {
        	nd.dist = INFINITY;
            nd.processed = false;
            nd.queued = false;
            nd.parent = null;
            nd.hops = Integer.MAX_VALUE;		//this field is used in the shortcut method to limit the dijkstra search
            nd.shortcutDist = INFINITY;			//used in shortcut to determine witness paths
            nd.active = true;					//used to remove the node during the shortcut search
        };
        
        Function<DpsNode, Long> getHeapKey = (nd) -> { return nd.dist; };
        
        Consumer<DpsNode> resetHeapR = (nd) -> {
        	nd.distR = INFINITY;
            nd.processedR = false;
            nd.queuedR = false;
            nd.parentR = null;
        };
        
        Function<DpsNode, Long> getHeapRKey = (nd) -> { return nd.distR; };
        
        Consumer<DpsNode> resetPreProc = (nd) -> {

            nd.queuedP = false;
            nd.processedP = false;
            nd.priority = 0;
            //these Node properties are persisted after buildPriorityQueue and shortcut(true) contraction
//            this.level = 0;
//            this.neighbors = 0;
//            this.rank = 0;
//            this.edgeDiff = 0;
        };
        
        Function<DpsNode, Integer> getPreProcKey = (nd) -> { return nd.priority; };
        
        this.heap = new PriorityNodeQ(ht, (a , b) ->  {return a.dist < b.dist;}, (a, b) -> {return a.dist == b.dist;}, resetHeap, getHeapKey, "heap", false);
        this.heapR = new PriorityNodeQ(ht, (a , b) -> {return a.distR < b.distR;}, (a, b) -> {return a.distR == b.distR;}, resetHeapR, getHeapRKey, "heapR", false);
        this.preProc = new PriorityNodeQ(ht, (a , b) -> {return a.priority < b.priority;}, (a, b) -> {return a.priority == b.priority;}, resetPreProc, getPreProcKey, "preProc", true);
        
        
        
        this.map = new HashMap<Integer,DpsNode>(n+1);
        
        
        
        for (int i = 0; i <= n; i++) {					//graph indexes and nodes are 1 indexed, index 0 not used
        	
            this.graph.add(new ArrayList<DpsEdge>());
            this.graphR.add(new ArrayList<DpsEdge>());          
  
        }

        this.scEdges = 0;
  
	}

	public DpsNode addNode(int i) {
		DpsNode n = new DpsNode(i);
		map.put(i, n);
		return n;
	}
	
	
	
	public void addEdges(int s, int t, int c) {		//(source, target, length in 1 based indexing)

		
		DpsNode source = map.get(s);
		DpsNode target = map.get(t);
		
		
		DpsEdge e = new DpsEdge(source, target, c);
        DpsEdge er = new DpsEdge(target, source, c);

        graph.get(s).add(e);
		graphR.get(t).add(er);

	}
	
	public void printGraph(boolean normal) {
		
		
		if(normal) {
			for(int i = 1; i < graph.size(); ++i) {
				System.out.print("Node: " + (i) + " ");
				System.out.print("[");
				for(DpsEdge e: graph.get(i)) {
					System.out.print(e.v.index + ",");
				}
				System.out.print("]");
				System.out.println();
			}
		} else {
			for(int i = 1; i < graphR.size(); ++i) {
				System.out.print("Node: " + (i + 1) + " ");
				for(DpsEdge e: graphR.get(i)) {
					System.out.print("[");
					System.out.print(e.v.index + ",");
					System.out.print("]");
					System.out.println();
				}
			}
		}
		
	}
	
	
	public DpsNode minPriority(DpsNode a, DpsNode b) {
		//returns minimum priority node with hash tiebreaker
		
		DpsNode min = null;
		
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

		
		heap.initializeQueue();
		heap.resetWorkingNodes();

		
		DpsNode sn = map.get(s);
		DpsNode tn = map.get(t);
		
		sn.dist = 0;
								
	
		heap.enQueue(sn);
		sn.queued = true;
		
		
		if(s == t) {
			return 0;									//I found myself!!
		}
		
		
	
		while(!heap.isEmpty()) {						//process the next node in the forward graph
			
			DpsNode processing = heap.getMin();
			processing.queued = false;
			
//				System.out.println("processing Node: " + processing.index);
			if(processing != tn) {
				for(DpsEdge e : graph.get(processing.index)) {
					
					DpsNode tt = e.v;	
					
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
	
	public int shortcut(DpsNode v, boolean contract, int hops) {
		
		//if create, create shortcuts else return the edge difference shortcuts - ins - outs
		
		System.out.println("Shortcut v: " + v.index + " contract: " + contract);
		
		int shortcuts = 0;
		
		Set<DpsNode> cover = new HashSet<DpsNode>();
		
		//a necessary priority term is the number of processed nodes during the initial search (cost of contraction)
		
		
		heap.initializeQueue();		
		
		ArrayList<DpsEdge> inEdges = graphR.get(v.index);
		ArrayList<DpsEdge> outEdges = graph.get(v.index);
		
		int ins = inEdges.size();
		int outs = outEdges.size();
		
		
		ArrayList<DpsNode> us = new ArrayList<DpsNode>();		//source nodes
		ArrayList<DpsNode> ws = new ArrayList<DpsNode>();		//target nodes
		
		Queue<DpsEdge> sCuts = new LinkedList<DpsEdge>();		//buffer for forward shortcuts to be added
		Queue<DpsEdge> sCutsR = new LinkedList<DpsEdge>();		//buffer for reverse shortcuts to be added

		
		for(DpsEdge e : inEdges) {			//get source nodes ignore contracted nodes
			if(!contract)
				us.add(e.v);
			else {
				if(!e.v.contracted)
					us.add(e.v);
			}
			
		}
		
		for(DpsEdge e : outEdges) {		//get target nodes
			if(!contract)
				ws.add(e.v);
			else {
				if(!e.v.contracted)
					ws.add(e.v);
			}
			
		}	
		
		System.out.println("us: " + us.size() + " ws: " + ws.size());
		
		v.active = false;				//remove v from the active graph
		
		for(DpsNode u : us) {				//check for a witness path to each target
				
			heap.resetWorkingNodes();	
			heap.initializeQueue();
			
			long mu = 0;
			long uvDist = 0;
			
			
			for(DpsEdge e : inEdges) {		//get d(u,v) for this source node
				if(e.v == u) {				
					uvDist = e.length;
				}
			}
			
			System.out.println("Processing u: " + u.index + " u-v dist: " + uvDist); 
			
			long maxShortcut = 0;
			
			System.out.println("Target nodes: ");
			
			for(DpsEdge e : outEdges) {					//for each target set the maximum shortcut distances d(u,v) + d(v,w) from this source
				e.v.shortcutDist = uvDist + e.length;
				maxShortcut = Math.max(maxShortcut, e.v.shortcutDist);
				e.v.dist = INFINITY;
				System.out.print(" " + e.v.index + "|" + e.v.shortcutDist );
			}

			System.out.println();
			
			u.dist = 0;
			u.hops = 0;

			
			long dijkstraStop = maxShortcut;
			
			System.out.println("Dijkstra Stopping Distance: " + dijkstraStop);

			heap.enQueue(u);
			u.queued = true;
			
			while(!heap.isEmpty()) {						//process the next node in the forward graph
				
				DpsNode x = heap.getMin();
				x.queued = false;
				mu = Math.max(mu, x.dist);
				
				System.out.println("popped node: " + x.index + " dist: " + x.dist + " mu: " + mu);
				
				if(mu >= dijkstraStop || x.hops > hops) {	//stopping conditions are max d(u,w) > max (d(u,v)+d(v,w)){if used: - min(x,w)} || hops > hops parameter
					System.out.println("Stopping Dijkstra mu: " + mu + " hops: " + x.hops);
					break;
				}
					

				for(DpsEdge e : graph.get(x.index)) {
					
					DpsNode tt = e.v;	
					
					if(tt.active) {							//ignore v node in witness search
						long td = x.dist + e.length;
						
						if(tt.dist > td) {
							if(ws.contains(tt)) {
								heap.working.add(tt);		//tt is a target node and is never enqueued
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
						}
					}
				}

				x.processed = true;

			}
			

			//Count and optionally generate shortcuts
			for(DpsNode w : ws) {
				
				System.out.println("Node: " + w.index + " u-w dist: " + w.dist + " shortcut dist: " + w.shortcutDist);
				
				if(w.dist > w.shortcutDist) {

					++shortcuts;			//number of shortcuts created for this node's contraction
					cover.add(w);
					
					
					if(contract) {
						System.out.println("Adding shortcut...");
						DpsEdge sc = new DpsEdge(u , v , w , w.shortcutDist);
						sCuts.add(sc);											//don't include the shortcuts in the witness searches!!!
						DpsEdge scr = new DpsEdge(w, v, u, w.shortcutDist);
						sCutsR.add(scr);
						++scEdges;												//total shortcuts in graph
					}
					
				}
			}
		}

		heap.working.add(v);								//size of working is contractionCost if needed to reduce preprocessing time	
		heap.resetWorkingNodes();
		
		v.shortcuts = shortcuts;							
		v.edgeDiff = shortcuts - ins - outs;				//set the edge difference NOTE: call setPriority on a node to update the priority property
		v.cover = cover.size();
		
		
		System.out.println("v Node: " + v.index + " edge diff: " + v.edgeDiff + " cover: "+ v.cover + " v shortcuts: " + v.shortcuts);
		System.out.println("Shortcut Edges: " + scEdges);
		
		if(contract) {
			v.contracted = true;
			v.rank = ++nRank;
			for(DpsNode w : ws) {									//increase contracted neighbors if shortcuts created
				if(shortcuts != 0)
					++w.neighbors;
				if(cover.contains(w))
					w.level = Math.max(w.level, v.level + 1); 		//update level of shortcut neighbors
			}
			for(DpsNode u : us) {
				if(shortcuts != 0)
					++u.neighbors;
			}
			while(!sCuts.isEmpty()) {							//add the shortcuts to the graph
				DpsEdge e = sCuts.remove();
				graph.get(e.u.index).add(e);
			}
			while(!sCutsR.isEmpty()) {
				DpsEdge e = sCutsR.remove();
				graphR.get(e.u.index).add(e);
			}
		}
		
		cover.clear();
		return v.edgeDiff;									//return the edge difference
	}
	
	public void buildPriorityQueue() {
		
//		int priority = 0;
		
		for(int i = 1; i < n; ++i) {
//			System.out.println();
//			System.out.println("Processing Node: " + i);
			DpsNode n = map.get(i);
			shortcut(n, false, maxHop);
			n.setPriority();
//			System.out.println(" Priority: " + n.priority);
			preProc.enQueue(n);
			n.queuedP = true;
//			preProc.printHeap();
		}
	}
	
	public void contractGraph() {
		int loops = 0;
		scEdges = 0;
		DpsNode min = null;
		System.out.println();
		System.out.println("Contracting Graph....");
		System.out.println();
		
		while(!preProc.isEmpty()) {
			
			DpsNode n = preProc.getMin();
			System.out.println();
			System.out.println("Pop: " + n.index + " priority: " + n.priority);
			shortcut(n, false, maxHop);
			n.setPriority();
			System.out.println("Reprioritized: " + n.index + " priority: " + n.priority);
			if(!preProc.isEmpty()) {
				min = preProc.peek();
				System.out.println("preProc min: " + min.index + " priority: " + min.priority);
			}
			System.out.println();
			if(!preProc.isEmpty() && n != minPriority(n, min)) {
				System.out.println("enQueue: " + n.index + " priority: " + n.priority);
				System.out.println();
				preProc.enQueue(n);
				continue;
			} else {
				shortcut(n, true, maxHop);
				n.setPriority();
				n.processedP = true;
				System.out.println("contracted: " + n.index + " priority: " + n.priority);
				System.out.println("-----------------------------------------------------");
			}
			preProc.printHeap();
			++loops;
			if(loops > 100) {
				System.out.println("Terminating with infinite loop...");
				break;
			}
				
				
		}
		System.out.println("Contraction Terminating with empty queue...");
		System.out.println("Contraction added edges: " + scEdges);
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
	
	
	private ArrayList<DpsNode> heap;
	private TableHash hTable;
	public HashSet<DpsNode> working;					//all nodes processed in this queue that must be reset
	private String name;
	
	//these are passed in as a Lambda functions to allow for different Node properties as keys
	private BiPredicate<DpsNode, DpsNode> minNode;		//test(a,b) -> true if a < b	
	private BiPredicate<DpsNode, DpsNode> equNode;		//test(a,b) -> true if a = b

	
	private Consumer<DpsNode> resetNode;				//resets the properties of nodes used in this queue
	private Function<DpsNode, ? extends Number> getKey;
	
	private boolean tieBreak;						//if tieBreak = true use the hash of the node index as a tiebreaker in sifts
	
	public PriorityNodeQ(
			TableHash hTable,
			BiPredicate<DpsNode,DpsNode> minNode,				//(a,b) -> a < b
			BiPredicate<DpsNode,DpsNode> equNode,				//(a,b) -> a == b
			Consumer<DpsNode> resetNode,					// a -> { a.field = x; a.field1 = y; ... }
			Function<DpsNode, ? extends Number> getKey,	// a -> a.key field
			String name,
			boolean tieBreak ) {
		
		this.heap = new ArrayList<DpsNode>();
		this.working = new HashSet<DpsNode>();
		
		this.hTable = hTable;
		this.minNode = minNode;
		this.equNode = equNode;
		this.resetNode = resetNode;
		this.getKey = getKey;
		this.name = name;
		this.tieBreak = tieBreak;
	}
	
	public void decreaseKey(DpsNode dn) {
	
		//this method assumes the property used as a key of the Node has been altered prior to the method call
		int dni = heap.indexOf(dn);
		minSiftUp(dni);
	}
	
	public void enQueue(DpsNode n) {
		
		working.add(n);				
		int end = heap.size();
		heap.add(end, n);
		minSiftUp(end);

	}
	
	public DpsNode getMin() {
		
		DpsNode nm = heap.get(0);
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
		
		if(heap.size() <= 1)
			return;
		
    	int mini = i;
    	int li = heap.size() - 1;	//0 based last index
    	int left = 2*i + 1;
    	int rght = 2*i + 2;
    	int lhi = 0;
    	int rhi = 0;
    	int minhi;
    	
    	if(!tieBreak) {
    		if(left <= li && minNode.test(heap.get(left), heap.get(mini))) {	//minTest returns (a,b) -> a < b Parent must be min
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
		
//		System.out.print(" minSiftUp i: " + i);
//		printHeap();
//		System.out.println();
		
		if(i > 0) {
			pi = (i - 1)/2;
			
			
			if(!tieBreak) {
				
				if(minNode.test(heap.get(i), heap.get(pi))) {			//minTest returns (a,b) -> a < b Parent must be min
//					System.out.println("swapping: " + pi + "-" + i);
					swap(pi, i);
					minSiftUp(pi);
				}
			} else {
				
				hi = hTable.hash(heap.get(i).index);					//use the hash of the index to break ties
	    		hpi = hTable.hash(heap.get(pi).index);
	    		
	    		//TODO: there is a small probability that the hash will not be unique ~ 4/3000 handle this case here
	    		
	    		if(minNode.test(heap.get(i), heap.get(pi))) {
//	    			System.out.println("swapping: " + pi + "-" + i);
					swap(pi, i);
					minSiftUp(pi);
				} else {
					if(equNode.test(heap.get(pi), heap.get(i))) {
						if(hpi > hi) {
//							System.out.println("tiebreak swapping: " + pi + "-" + i);
							swap(pi, i);
							minSiftUp(pi);
						}
					}
				}	
			}
			
			
		}

	}

	public DpsNode peek() {
		return heap.get(0);
	}
	
	public void printHeap() {
		System.out.println();
		System.out.print(name + ": ");
		for(DpsNode n : heap) {
			System.out.print(heap.indexOf(n) + "|" + n.index + ":" + getKey.apply(n) + ", ");
		}
		System.out.println();
	}
	
	
	private void swap(int i, int n) {
    	
    	DpsNode tn = heap.get(i);
    	heap.set(i, heap.get(n));
    	heap.set(n, tn);
    	
    }
	
	public void resetWorkingNodes() {
		
		if(!working.isEmpty()) {
			for(DpsNode n: working) {
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
    				int n = Integer.valueOf(params[0]);			//n=number of nodes
    				int m = Integer.valueOf(params[1]);			//m=number of edges
    				
    				
    				
    				BiGraph g = new BiGraph(n, tHash);			//was n+1
    				
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
    				
    				System.out.println();
    				g.printGraph(true);
    				System.out.println();
    				
    				System.out.println("Building Priority Queue: ");
    				g.buildPriorityQueue();
    				System.out.println("Contracting: =========================================================");
    				g.contractGraph();
    				System.out.println("Ready");
    				
    				
    				
    				//TODO: rewrite this section to run queries and compare results to test files
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
