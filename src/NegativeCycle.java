import java.awt.List;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;

class Node {
	public int index;
	public long dist;
	public boolean visited;
	public int pindex;
}

class Edge {
	public int target;
	public int length;
	
	public Edge(int t, int l) {
		this.target = t;
		this.length = l;
	}
}

class Graph {
	
	public ArrayList<ArrayList<Edge>> graph;		//adjacency list of edges
	public ArrayList<Node> map;						//returns Nodes by vertex (index)
	public Queue<Node> queue;
	
	
	
	public int root;
	public int n;									//|V| number of nodes
	public int m;									//|E| number of edges

	
	public Graph(int n, int m) {
		this.n = n;
		this.m = m;

        this.graph = new ArrayList<ArrayList<Edge>>();
        
        this.queue = new LinkedList<Node>();
        
        this.map = new ArrayList<Node>();
        
        for (int i = 0; i < n; i++) {
        	
            this.graph.add(new ArrayList<Edge>());
            
            Node node = new Node();
            node.index = i;
            node.dist = Long.MAX_VALUE;
            node.visited = false;
            node.pindex = -1;
            
            
            map.add(i, node);
        }
        
        this.root = 0;		//arbitrary selection of the first node as the root

	}
	
	private boolean fordBellman(Node node) {
		
		queue.clear();
		queue.add(node);
		
		boolean noChange = true;
		
		//TODO: add procedure to handle disjoint graphs
		int nodes = 0;
		
		while(!queue.isEmpty()) {
			Node u = queue.remove();
			++nodes;
			int ui = u.index;
			for(Edge e : graph.get(ui)) {
				if(map.get(e.target).dist > u.dist + e.length) {
					map.get(e.target).dist = u.dist + e.length;
					noChange = false;
				}
				queue.add(map.get(e.target));
			}
		}
		
		return noChange;
		
	}
	
	

	public int negativeCycle() {
       
		Node node = map.get(root);
		node.dist = 0;
		boolean noChange = false;
		
		if(n == 1)
			return 0;
		
		int cycles = n - 1;
		
		do {
			noChange = fordBellman(node);
			--cycles;
			
		} while (cycles > 0 && !noChange);
		
		return fordBellman(node)? 1 : 0;

    }

}


public class NegativeCycle {
    

    public static void main(String[] args) {
    	Scanner scanner = new Scanner(System.in);
        
        int n = scanner.nextInt();
        int m = scanner.nextInt();
        
        Graph g = new Graph(n , m);
        
        for (int i = 0; i < m; i++) {
            int x, y, w;
            x = scanner.nextInt() - 1;		//shift indexes to 0 based
            y = scanner.nextInt() - 1;
            w = scanner.nextInt();

            Edge e = new Edge(y, w);
            g.graph.get(x).add(e);
            
            
        }
        System.out.println(g.negativeCycle());
        scanner.close();
    }
}

