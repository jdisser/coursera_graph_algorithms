import java.util.Scanner;
import java.util.ArrayList;


class Node {
	public int index;
	public long dist;			//this is the distance from the source in the graph
	public long distR;			//this is the distance from the target in the reverse graph
	public boolean visited;
	public boolean visitedR;
	public int pindex;
	
	public Node(int i) {
        this.index = i;
        this.dist = Long.MAX_VALUE;
        this.distR = Long.MAX_VALUE;
        this.visited = false;
        this.visitedR = false;
        this.pindex = -1;
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

	
	public int root;
	public int n;

	
	public BiGraph(int n) {
		this.n = n;

        this.graph = new ArrayList<ArrayList<Edge>>();
        this.graphR = new ArrayList<ArrayList<Edge>>();
        
        this.heap = new ArrayList<Node>();
        this.heapR = new ArrayList<Node>();
        
        this.map = new ArrayList<Node>();
        
        for (int i = 0; i < n; i++) {
        	
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
		
		for(Node n : map) {
			n.dist = Long.MAX_VALUE;
			n.distR = Long.MAX_VALUE;
			n.visited = false;
			n.visitedR = false;
			n.pindex = -1;
			
            heap.add(n);			//was (n.index, n)
            heapR.add(n);

			
		}
	}
	
	
	public void addEdges(int x, int y, int c) {		//(source, target, length)

		x -= 1;
		y -= 1;								//raw vertexes are 1 based, map indices are 0 based
		
		Edge e = new Edge(x, y, c);
        Edge er = new Edge(y, x, c);

        graph.get(x).add(e);
		graphR.get(y).add(er);

	}
	
	
	
	private void minSiftDown(int i, ArrayList<Node> h) {
		
//		System.out.println("minSiftDown 1: " + i);
//		printHeap();
		
		
		
    	int mini = i;
    	int li = h.size() - 1;	//0 based last index
    	int left = 2*i + 1;
    	int rght = 2*i + 2;
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
			
			if(h.get(pi).dist > h.get(i).dist) {
				swap(pi, i, h);
				minSiftUp(pi, h);
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
	
	private void decreaseKey(int i, long d, ArrayList<Node> h, boolean reverse) {
		
		System.out.println(" decreaseKey i: " + i + " d: " + d + " reverse: " + reverse);
		Node dn = map.get(i);
		if(reverse) {
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
	
	public long biDijkstra(int s, int t) {
		
		//implements Dijkstra's algorithm

		initializeQueues();
		
		decreaseKey(s, 0, heap, false);			//start the algorithm on this node in the forward direction
		decreaseKey(t, 0, heapR, true);			//and from this one in the reverse direction (Reverse Graph => true)
		
		long result = -1;						//if after processing all nodes in the graph this is unchanged the target is unreachable
		
		if(s == t) {
			return 0;							//I found myself!!
		}
		
		while(!heap.isEmpty() || !heapR.isEmpty()) {
			
			
			if(!heap.isEmpty()) {						//process the next node in the forward graph
				
				Node r = getMin(heap);

//				System.out.println("processing Node: " + r.index);

				for(Edge e : graph.get(r.index)) {
					
					if(r.dist == Long.MAX_VALUE)		//this is an unreachable node if it's dist is infinite
						break;
					
					if(map.get(e.target).dist > r.dist + e.length) {
						decreaseKey(e.target, r.dist + e.length, heap, false);
						map.get(e.target).pindex = r.index;				//min path is my daddy
					}
				}
				
				//TODO: handle the case where friends are only reachable in one direction (distx = infinity)
				r.visited = true;
				if(r.visitedR == true) {				//stop when a node has been processed from both ends
					result = r.dist + r.distR;
					break;					
				}
					
			}
			
			if(!heapR.isEmpty()) {						//process the next node in the reverse graph
				Node rr = getMin(heapR);
				
				for(Edge er : graphR.get(rr.index)) {
					
					if(rr.distR == Long.MAX_VALUE)
						break;
					
					if(map.get(er.target).distR > rr.distR + er.length) {
						decreaseKey(er.target, rr.distR + er.length, heapR, true);
						map.get(er.target).pindex = rr.index;
					}	
				}
				rr.visitedR = true;
				if(rr.visited == true) {
					result = rr.dist + rr.distR;
					break;
				}
			}
			
		}
		
		return result;
    }

}




public class FriendSuggestion {


    public static void main(String args[]) {
        Scanner in = new Scanner(System.in);
        int n = in.nextInt();
        int m = in.nextInt();
        
        BiGraph g = new BiGraph(n);

        for (int i = 0; i < m; i++) {
            int x, y, c;
            x = in.nextInt();
            y = in.nextInt();
            c = in.nextInt();


            g.addEdges(x, y, c);         

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
