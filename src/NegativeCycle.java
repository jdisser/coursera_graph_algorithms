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
	public int source;
	public int length;
	
	public Edge(int t, int s, int l) {
		this.target = t;
		this.source = s;
		this.length = l;
	}
}

class Graph {
	
	public ArrayList<Edge> graph;					//list of edges
	public ArrayList<Node> map;						//returns Nodes by vertex (index)
	public Queue<Node> queue;
	
	
	

	public int n;									//|V| number of nodes
	public int m;									//|E| number of edges

	
	public Graph(int n, int m) {
		this.n = n;
		this.m = m;

        this.graph = new ArrayList<Edge>();
        
        this.queue = new LinkedList<Node>();
        
        this.map = new ArrayList<Node>();
        
        for (int i = 0; i < n; i++) {
 
            Node node = new Node();
            node.index = i;
            node.dist = Long.MAX_VALUE;
            node.visited = false;
            node.pindex = -1;
           
            map.add(i, node);
        }

	}
		

	
	private boolean bellmanFord() {
		
		//this is an adaptation of the bellman-ford algorithm to specifically detect a negative cycle
		//in a disjoint graph, it does not generate a shortest path nor does it have a true source
		//from which to measure distance
		
		boolean noChange;
		
		if(queue.size() == 0)
			return true;		//if there are no negative edges there can not be negative cycles
		
		int cycles = n;
		
		do {
			
			noChange = true;
			
			for(Edge e : graph) {
				Node u = map.get(e.source);		//source node
				Node v = map.get(e.target);		//target node
				
				if(v.dist > u.dist + e.length) {
					v.dist = u.dist + e.length;
					noChange = false;
				}	
				--cycles;
			}
			
		} while (cycles > 0 && !noChange);
		
		return noChange;
	}

	public int negativeCycle() {

		return bellmanFord()? 0 : 1; 	//if there is no change after n-1 + 1 cycles (or less) then there are 0 neg cycles
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

            Edge e = new Edge(y, x, w);
            g.graph.add(e);
            
            if(e.length < 0) {
            	g.map.get(x).dist = 0;		//initialize the source of negative edges to 0
            	g.queue.add(g.map.get(x));	//save to count negative edges
            }
            
        }
        System.out.println(g.negativeCycle());
        scanner.close();
    }
}

