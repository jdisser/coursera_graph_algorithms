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
	public boolean active;
	public boolean contracted;
	public int pindex;
	public int pindexR;
	public long level;
	public long importance;
	
	
	
	public Node(int i) {
        this.index = i;
        this.dist = INFINITY;
        this.distR = INFINITY;
        this.processed = false;
        this.processedR = false;
        this.queued = false;
        this.queuedR = false;
        this.contracted = false;
        this.pindex = -1;
        this.pindexR  = -1;
        this.active = true;
        this.level = 0;
        this.importance = 0;
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

	
	public static int nodeNumber(int i, int j, int w) {	//1 index position of node in test graph w x h
		return w*(j - 1) + i;
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
        
        this.preProc = new ArrayList<Node>();
        
        
        
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
	
	
	
	
	
	
	
	public Node addNode(int i) {
		Node n = new Node(i);
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
//		System.out.println("BiAStar double processed nodes: " + dblProcessed);
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

class tableHash{
	
	ArrayList<ArrayList<Integer>> hashTable;
	ArrayList<Integer> masks;
	
	int bits;
	int r;
	int k;
	
	
	public tableHash(int bits, int r) {
		
		this.bits = bits;
		this.r = r;
		
		this.hashTable = new ArrayList<ArrayList<Integer>>();
		for(int i = 0; i < r; ++i) {
			this.hashTable.add(new ArrayList<Integer>());
		}
		
		int m = 2 ^ bits;
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
    private static class Impl {
        // See the descriptions of these fields in the starter for friend_suggestion
        int n;
        ArrayList<Integer>[][] adj;
        ArrayList<Long>[][] cost;
        Long[][] distance;
        ArrayList<PriorityQueue<Entry>> queue;
        boolean[] visited;
        ArrayList<Integer> workset;
        final Long INFINITY = Long.MAX_VALUE / 4;
 
        // Position of the node in the node ordering
        Integer[] rank;
        // Level of the node for level heuristic in the node ordering
        Long[] level;

        Impl(int n) {
            this.n = n;
            visited = new boolean[n];
            Arrays.fill(visited, false);
            workset = new ArrayList<Integer>();
            rank = new Integer[n];
            level = new Long[n];
            distance = new Long[][] {new Long[n], new Long[n]};
            for (int i = 0; i < n; ++i) {
                distance[0][i] = distance[1][i] = INFINITY;
                level[i] = 0L;
                rank[i] = 0;
            }
            queue = new ArrayList<PriorityQueue<Entry>>();
            queue.add(new PriorityQueue<Entry>(n));
            queue.add(new PriorityQueue<Entry>(n));
        }

        // Preprocess the graph
        void preprocess() {
            // This priority queue will contain pairs (importance, node) with the least important node in the head
            PriorityQueue<Entry> q = new PriorityQueue<Entry>(n);
            // Implement this method yourself
        }

        void add_edge(int side, int u, int v, Long c) {
            for (int i = 0; i < adj[side][u].size(); ++i) {
                int w = adj[side][u].get(i);
                if (w == v) {
                    Long cc = min(cost[side][u].get(i), c);
                    cost[side][u].set(i, cc);
                    return;
                }
            }
            adj[side][u].add(v);
            cost[side][u].add(c);
        }

        void apply_shortcut(Shortcut sc) {
            add_edge(0, sc.u, sc.v, sc.cost);
            add_edge(1, sc.v, sc.u, sc.cost);
        }

        void clear() {
            for (int v : workset) {
                distance[0][v] = distance[1][v] = INFINITY;
                visited[v] = false;
            }
            workset.clear();
            queue.get(0).clear();
            queue.get(1).clear();
        }

        void mark_visited(int u) {
            visited[u] = true;
            workset.add(u);
        }

        // See the description of this method in the starter for friend_suggestion
        boolean visit(int side, int v, Long dist) {
            // Implement this method yourself
            return false;
        }                

        // Add the shortcuts corresponding to contracting node v. Return v's importance.
        Long shortcut(int v) {
            // Implement this method yourself

            // Compute the node importance in the end
            Long shortcuts = 0;
            Long vlevel = 0L;
            Long neighbors = 0L;
            Long shortcutCover = 0L;
            // Compute the correct values for the above heuristics before computing the node importance
            Long importance = (shortcuts - adj[0][v].size() - adj[1][v].size()) + neighbors + shortcutCover + vlevel;
            return importance;
        }

        // Returns the distance from s to t in the graph
        Long query(int s, int t) {
            if (s == t) {
                return 0L;
            }
            visit(0, s, 0L);
            visit(1, t, 0L);
            Long estimate = INFINITY;
            // Implement the rest of the algorithm yourself
            return estimate == INFINITY ? -1 : estimate;            
        }

        class Entry implements Comparable<Entry>
        {
            Long cost;
            int node;
          
            public Entry(Long cost, int node)
            {
                this.cost = cost;
                this.node = node;
            }
         
            public int compareTo(Entry other)
            {
                if (cost == other.cost) {
                    return node < other.node ? -1 : node > other.node ? 1: 0;
                }
                return cost < other.cost ? -1 : cost > other.cost ? 1 : 0;
            }
        }

        class Shortcut
        {
            int u;
            int v;
            Long cost;

            public Shortcut(int u, int v, Long c)
            {
                this.u = u;
                this.v = v;
                cost = c;
            }
        }
    }

    
    //original main method in starter
    public static void orgMain(String args[]) {
        Scanner in = new Scanner(System.in);
        int n = in.nextInt();
        int m = in.nextInt();
        Impl ch = new Impl(n);
        @SuppressWarnings("unchecked")
        ArrayList<Integer>[][] tmp1 = (ArrayList<Integer>[][])new ArrayList[2][];
        ch.adj = tmp1;
        @SuppressWarnings("unchecked")
        ArrayList<Long>[][] tmp2 = (ArrayList<Long>[][])new ArrayList[2][];
        ch.cost = tmp2;
        for (int side = 0; side < 2; ++side) {
            @SuppressWarnings("unchecked")
            ArrayList<Integer>[] tmp3 = (ArrayList<Integer>[])new ArrayList[n];
            ch.adj[side] = tmp3;
            @SuppressWarnings("unchecked")
            ArrayList<Long>[] tmp4 = (ArrayList<Long>[])new ArrayList[n];
            ch.cost[side] = tmp4;
            for (int i = 0; i < n; i++) {
                ch.adj[side][i] = new ArrayList<Integer>();
                ch.cost[side][i] = new ArrayList<Long>();
            }
        }

        for (int i = 0; i < m; i++) {
            int x, y;
            Long c;
            x = in.nextInt();
            y = in.nextInt();
            c = in.nextLong();
            ch.adj[0][x - 1].add(y - 1);
            ch.cost[0][x - 1].add(c);
            ch.adj[1][y - 1].add(x - 1);
            ch.cost[1][y - 1].add(c);
        }

        ch.preprocess();
        System.out.println("Ready");

        int t = in.nextInt();

        for (int i = 0; i < t; i++) {
            int u, v;
            u = in.nextInt();
            v = in.nextInt();
            System.out.println(ch.query(u-1, v-1));
        }
    }
    
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
