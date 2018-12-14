import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashSet;


class Node {
	private static final long INFINITY = Long.MAX_VALUE/4;
	public int index;
	public long dist;			//this is the distance from the source in the graph
	public long distR;			//this is the distance from the target in the reverse graph
	public boolean visited;
	public boolean visitedR;
	public int pindex;
	
	public Node(int i) {
        this.index = i;
        this.dist = INFINITY;
        this.distR = INFINITY;
        this.visited = false;
        this.visitedR = false;
        this.pindex = -1;
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
	public int length;
	
	public Edge(int s, int t, int l) {
		this.target = t;
		this.source = s;
		this.length = l;
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
            
            Node node = new Node(i);
            map.add(i, node);
        }
        
        this.root = -1;

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
		}
	}
	
	public void enQueue(int v, ArrayList<Node> h) {
		
		int end = h.size();
		Node n = map.get(v);
		h.add(end, n);
		minSiftUp(end, h);
		working.add(n);
		
	}
	
	
	public void addEdges(int s, int t, int c) {		//(source, target, length in 1 based indexing)

		s -= 1;
		t -= 1;								//raw vertexes are 1 based, map indices are 0 based
		
		Edge e = new Edge(s, t, c);
        Edge er = new Edge(t, s, c);

        graph.get(s).add(e);
		graphR.get(t).add(er);

	}

	
	public void genTestEdges(int c, int k, int e) {
		//assumes that a graph and graphR of size c*k exists
		//generates a grid of nodes
		//c = # columns
		//k = # nodes in column
		//e = target nodes in columns c+1 & c-1 +/-e and at same height (e <= k/2)
		//distance between a node in the first column (i <= k) and the node at the same height in the last column (ii = i + (c-1)*k) is c -1
		//let l(e) = 2e + 1
		//then the number of edges processed is < 2( l(e) + 2l(e)^2 + l(e)^3) for c=5 
		//for c=5, k=8, e=2 nodes processed are 2( 5 + 2*25 + 125) < 360
		//for c=5, k=30 it's 2(30 + 2*900 + 27000) < 58000
		
		System.out.println("generating test edges...");
		
		int edges = 0;
		
		//generate the forward facing edges
		for(int i = 1; i <= (c-1)*k - e; ++i) {					//the parameters are index 1 based like the scanned input
			addEdges(i, i + k, 1);							
			for(int j = 1; j<= e; ++j) {
				addEdges(i, i + (k+j), 1);
				addEdges(i, i + (k-j), 1);
			}
			edges += 2*e + 1;
		}
		
		//generate the backward facing edges
		for(int i = c*k; i > k + e; --i) {
			addEdges(i, i - k, 1);
			for(int j = 1; j <= e; ++j) {
				addEdges(i, i - (k-e), 1);
				addEdges(i, i - (k+e), 1);
			}
			edges += 2*e + 1;
		}
		
		System.out.println("generated test edges: " + edges);
		
		
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
	
	public long biDijkstra(int s, int t) {	//s & t are 0 indexed
		
		//implements meet in the middle bidirectional Dijkstra's algorithm

		long start = System.nanoTime();
		
		int processed = 0;
		
		initializeQueues();
		
		resetWorkingNodes();
		
		enQueue(s, heap);
		enQueue(t, heapR);
		
		decreaseKey(s, 0, heap);			//start the algorithm on this node in the forward direction
		decreaseKey(t, 0, heapR);			//and from this one in the reverse direction (Reverse Graph => heapR)
		
		long result = -1;						//if after processing all nodes in the graph this is unchanged the target is unreachable
		
		if(s == t) {
			return 0;							//I found myself!!
		}
		
		while(!heap.isEmpty() || !heapR.isEmpty()) {
			
			
			if(!heap.isEmpty()) {						//process the next node in the forward graph
				
				Node r = getMin(heap);

//				System.out.println("processing Node: " + r.index);

				for(Edge e : graph.get(r.index)) {
					
					if(r.dist == INFINITY)		//this is an unreachable node if it's dist is infinite
						break;
					
					Node tt = map.get(e.target);
					
					working.add(tt);				
					
					if(tt.dist > r.dist + e.length) {
						enQueue(tt.index, heap);
						decreaseKey(tt.index, r.dist + e.length, heap);
						tt.pindex = r.index;				//min path is my daddy
					}
					++processed;
				}
				
				
				r.visited = true;
				
				if(r.visitedR == true) {				//stop when a node has been processed from both ends
					result = r.dist + r.distR;
					break;					
				}
					
			}
			
			if(!heapR.isEmpty()) {						//process the next node in the reverse graph
				
				Node rr = getMin(heapR);
				
//				System.out.println("processing Node: " + rr.index);
				
				for(Edge er : graphR.get(rr.index)) {
					
					if(rr.distR == INFINITY)
						break;
					
					Node ttr = map.get(er.target);
					
					working.add(ttr);
					
					
					if(ttr.distR > rr.distR + er.length) {
						enQueue(ttr.index, heapR);
						decreaseKey(ttr.index, rr.distR + er.length, heapR);
						ttr.pindex = rr.index;
					}
					++processed;
				}
				rr.visitedR = true;
				if(rr.visited == true) {
					result = rr.dist + rr.distR;
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




public class FriendSuggestion {



	
	
    public static void main(String args[]) {
    	
    	if(args.length != 0) {

    		
    		int c = 5;
    		int k = 100000;
    		int e = 3;
    		
    		BiGraph g = new BiGraph(c*k);
    		
    		long start = System.nanoTime();
    		g.genTestEdges(c,k,e);			//see function for description of parameters
    		long finish = System.nanoTime();
    		long elapsed = (finish - start)/1000000;
    		
    		System.out.println("Edge generation time: " + elapsed + " ms" );
    		
    		start = System.nanoTime();
    		
    		for(int i = 100; i < 111; ++i) {
    			long d = g.biDijkstra(i, i + (c-1)*k);
    			System.out.println("i: " + i + " d: " + d);
    		}
    		
    		finish = System.nanoTime();
    		elapsed = (finish - start) / 1000000;
    		
    		System.out.println("Query time: " + elapsed + " ms" );
    		
    		
    		//TODO: run a series of queries
    		//TODO: Print the time it takes for each query and total (including reseting and edge generation)
    		
    		
    	} else {
    		Scanner in = new Scanner(System.in);
            int n = in.nextInt();
            int m = in.nextInt();
            
            BiGraph g = new BiGraph(n);

            for (int i = 0; i < m; i++) {
                int x, y, c;
                x = in.nextInt();
                y = in.nextInt();
                c = in.nextInt();


                g.addEdges(x, y, c);    //edges are added as 1 indexed     

            }
            

            int t = in.nextInt();

            
            //TODO: determine if the output should print after each pair or all together after the last pair
            for (int i = 0; i < t; i++) {
                int u, v;
                u = in.nextInt();
                v = in.nextInt();
              
                System.out.println(g.biDijkstra(u-1, v-1));
            }
            in.close();
    	}
    	
        
    }
}
