import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;
import static java.lang.Math.pow;

class Node {
	public int index;
	public double x,y;
	public int pindex;
	
	public Node (int x, int y) {
		this.x = (double) x;
		this.y = (double) y;
		this.index = -1;
		this.pindex = -1;
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
	
	public ArrayList<Node> graph;
	public ArrayList<Edge> heap;
	public ArrayList<Edge> xMst;
	
	public Graph () {
		this.graph = new ArrayList<Node>();
		this.heap = new ArrayList<Edge>();
		this.xMst = new ArrayList<Edge>();
	}
	
	public int addNode(Node v) {
		
		int vi;
		
		graph.add(v);
		vi = graph.indexOf(v);
		graph.get(vi).index = vi;
		
		return vi;	
	}
	
	public void setEdges() {
		
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
	
	
}





public class ConnectingPoints {
    private static double minimumDistance(int[] x, int[] y) {
        double result = 0.;
        //write your code here
        return result;
    }

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
        System.out.println(minimumDistance(x, y));
    }
}

