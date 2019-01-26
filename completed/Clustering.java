import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.Scanner;

class Node {
	public int index;
	public double x,y;
	public int pindex;
	public int rank;
	
	public Node (int x, int y) {
		this.x = (double) x;
		this.y = (double) y;
		this.index = -1;
		this.pindex = -1;
		this.rank = 0;
	}
}

class Edge {
	public Node v1;
	public Node v2;
	public double len;
	
	public Edge(Node v1, Node v2) {
		this.v1 = v1;
		this.v2 = v2;
		this.len = setLength();
	}
	
	private double setLength() {
		//lengths are stored as sum of squares and must be converted to length
		return pow((v1.x - v2.x), 2.0) + pow((v1.y - v2.y) , 2.0);
	}
}

class Graph {
	
	public ArrayList<Node> graph;		//the disjoint set of all vertices in the graph
	public ArrayList<Edge> heap;		//a min priority queue that is indexed on edge length
	public ArrayList<Edge> xMst;		//the set of edges that are members of the Minimum Spanning Tree
	
	public int trees;
	
	public Graph () {
		this.graph = new ArrayList<Node>();
		this.heap = new ArrayList<Edge>();
		this.xMst = new ArrayList<Edge>();
		this.trees = 0;
	}
	
	public int addNode(Node v) {
		
		int vi;
		
		graph.add(v);
		vi = graph.indexOf(v);
		v.index = vi;
		v.pindex = vi;			//initially all nodes are a tree
		++trees;
		
		return vi;	
	}
	
	public void printGraph() {
		
		System.out.println("Graph: ");
		for(Node n : graph) {
			System.out.print(" " + n.x + "|" + n.y + " i:" + n.index + " p:" + n.pindex + " ");
		}
        System.out.println("");
        System.out.println("--------------------------");
	}
	
	
	public void printEdge(Edge e) {
		System.out.print(" " + e.v1.index + "-" + e.v2.index + ":" + e.len + " ");
	}
	
	
	
	public void setEdges() {
		//generate the set of fully connected edges for the undirected graph
		
		int vi;
		
		for(Node v : graph) {
			vi = graph.indexOf(v);
			for(int j = vi +1; j < graph.size(); ++j) {
				Edge e = new Edge(v, graph.get(j));
				addEdge(e);
//				printEdge(e);
			}
		}
		
		vi = heap.size();
//		System.out.println("generated edges: " + vi);
		
	}
	
	private void minSiftDown(int i) {
		
//		System.out.println("minSiftDown 1: " + i);
//		printHeap();
		
    	int mini = i;
    	int li = heap.size() - 1;	//0 based last index
    	int left = 2*i + 1;
    	int rght = 2*i + 2;
    	if(left <= li && heap.get(left).len < heap.get(mini).len) {
    		mini = left;
    	}
    	if(rght <= li && heap.get(rght).len < heap.get(mini).len) {
    		mini = rght;
    	}
    	if(mini != i) {
	   		swap(i, mini);
	   		minSiftDown(mini);
    	}	
    }
	
	private void swap(int i, int n) {
    	
    	Edge te = heap.get(i);
    	heap.set(i, heap.get(n));
    	heap.set(n, te);
    	
    }
	
	private void minSiftUp(int i) {
		
		int pi;
//		System.out.println(" minSiftUp i: " + i);
//		printHeap();
		if(i > 0) {
			pi = (i - 1)/2;
			
//			System.out.println(" minSiftUp pi: " + pi);
			
			if(heap.get(pi).len > heap.get(i).len) {
				swap(pi, i);
				minSiftUp(pi);
			}
			
		}

	}
	
	private Edge getMin() {
		
		Edge nm = heap.get(0);
		swap(0, heap.size() - 1);
		heap.remove(heap.size() - 1);
		minSiftDown(0);
		
		return nm;
	}
	
	private void addEdge(Edge e) {
		heap.add(e);
		minSiftUp(heap.indexOf(e));
	}
	
	//Disjoint set methods for nodes using path compression and rank heuristics
	
	private int find(int i) {
		
		Node n = graph.get(i);
		
		if(n.pindex != n.index) {
			i = n.pindex;
			n.pindex = find(i);
		}
		
		return n.pindex;		
	}
	
	private void union(Node n1, Node n2) {
		
		
		int p1, p2;
		
		p1 = n1.index;
		p2 = n2.index;
		
		Node t1 = graph.get(find(p1));
		Node t2 = graph.get(find(p2));
		
		
		
		if(t1.rank > t2.rank) {
			t2.pindex = t1.index;
		} else {
			t1.pindex = t2.index;
		}
		
		if(t1.rank == t2.rank) {
			++t2.rank;
		}
		
//		printGraph();

	}
	
	private Edge generateMST(int k) {				//k is number of clusters
		
		setEdges();
		
		Edge e = null;
		
//		System.out.println("Generating MST:");
		
		while(!heap.isEmpty()) {
			e = getMin();
			if(find(e.v1.index) != find(e.v2.index)) {

				if(trees == k)		//if this edge will reduce trees < K it has the maximum cluster distance length
					break;
				
				xMst.add(e);
				union(e.v1,e.v2);
				--trees;			//each union reduces the number of clusters or trees
//				printEdge(e);
			}		
		}
		
		return e;
	}
	
	public double clusterDistance(int k) {
		Edge e = generateMST(k);
      
//        int en = xMst.size();
//        System.out.println("Number of edges: " + en);
        
        return sqrt(e.len);
    }
	
}








public class Clustering {
    

    public static void main(String[] args) {
    	Scanner scanner = new Scanner(System.in);
        int n = scanner.nextInt();

        Graph V = new Graph();
        
        int x,y;
        
        
        for (int i = 0; i < n; i++) {
            x = scanner.nextInt();
            y = scanner.nextInt();
            
            Node v = new Node(x, y);
            V.addNode(v);
        }
        
        
//        V.printGraph();
        int k = scanner.nextInt();
        System.out.println(V.clusterDistance(k));
    }
}

