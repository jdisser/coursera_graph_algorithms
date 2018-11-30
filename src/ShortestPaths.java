import java.util.*;

class Node {
	public int index;
	public long dist;
	public boolean shortest;
	public boolean reachable;
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
            node.reachable = false;			//node is reachable
            node.shortest = true;			//shortest path exists
            node.pindex = -1;
           
            map.add(i, node);
        }

	}
		

	
	private boolean bellmanFord(int s) {
		
		Node ns = map.get(s);
		
		ns.dist = 0;	//set source distance to 0
		ns.shortest = true;
		ns.reachable = true;
		
		
		boolean noChange = true;

		int cycles = n;		//relax the matrix
							//on the nth cycle if shortest is false the shortest path does not exist
							//if reachable is false (and dist = MAX_LONG) the node isn't reachable
							//if it is reachable dist contains the shortest path
		
		do {
			
			noChange = true;
			
//			System.out.println("cycle: " + cycles);
			
			for(Edge e : graph) {
				Node u = map.get(e.source);		//source node
				Node v = map.get(e.target);		//target node
				
				v.shortest = true;
				
//				System.out.print("e: [" + e.source + "," + e.target + "," + e.length + "] ");
//				System.out.print("u: [" + u.index + "," + u.dist + "] ");
//				System.out.print("v: [" + v.index + "," + v.dist + "] ");
				
				if(u.dist < Long.MAX_VALUE) {	//if the source is not infinite
					if(v.dist > u.dist + e.length) {
						v.dist = u.dist + e.length;
						v.reachable = true;
						v.shortest = false;
						
						noChange = false;
//						System.out.print("noChange: " + noChange + " --> ");
//						System.out.print("v: [" + v.index + "," + v.dist + "] ");
					}
				}
				
//				System.out.println("");	
			}
			
			--cycles;
			
		} while (cycles > 0 && !noChange);
		
		return noChange;
	}

	
	
	public void shortestPaths(int s) {
	      bellmanFord(s);
	    }
}



public class ShortestPaths {

    

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

        }
        int s = scanner.nextInt() - 1;
        
        scanner.close();
        
        g.shortestPaths(s);
        
        for (int i = 0; i < n; i++) {
            if (g.map.get(i).reachable == false) {
                System.out.println('*');
            } else if (g.map.get(i).shortest == false) {
                System.out.println('-');
            } else {
                System.out.println(g.map.get(i).dist);
            }
        }
    }

}

