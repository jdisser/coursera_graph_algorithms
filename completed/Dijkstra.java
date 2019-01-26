import java.util.*;



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
	public ArrayList<Node> heap;					//priority queue with method to update key values
	public ArrayList<Node> map;						//returns Nodes by vertex (index)
	
	
	public int root;
	public int n;

	
	public Graph(int n) {
		this.n = n;

        this.graph = new ArrayList<ArrayList<Edge>>();
        
        this.heap = new ArrayList<Node>();
        
        this.map = new ArrayList<Node>();
        
        for (int i = 0; i < n; i++) {
        	
            this.graph.add(new ArrayList<Edge>());
            
            Node node = new Node();
            node.index = i;
            node.dist = Long.MAX_VALUE;
            node.visited = false;
            node.pindex = -1;
            
            heap.add(i, node);
            map.add(i, node);
        }
        
        this.root = -1;

	}
	
	private void minSiftDown(int i) {
		
//		System.out.println("minSiftDown 1: " + i);
//		printHeap();
		
    	int mini = i;
    	int li = heap.size() - 1;	//0 based last index
    	int left = 2*i + 1;
    	int rght = 2*i + 2;
    	if(left <= li && heap.get(left).dist < heap.get(mini).dist) {
    		mini = left;
    	}
    	if(rght <= li && heap.get(rght).dist < heap.get(mini).dist) {
    		mini = rght;
    	}
    	if(mini != i) {
	   		swap(i, mini);
	   		minSiftDown(mini);
    	}
   		
    }
	
	private void printHeap() {
		System.out.print("heap: ");
		for(Node n : heap) {
			System.out.print(n.index + ":" + n.dist + " ");
		}
		System.out.println(" ");
	}
	
	
	
	private void minSiftUp(int i) {
		
		int pi;
//		System.out.println(" minSiftUp i: " + i);
//		printHeap();
		if(i > 0) {
			pi = (i - 1)/2;
			
//			System.out.println(" minSiftUp pi: " + pi);
			
			if(heap.get(pi).dist > heap.get(i).dist) {
				swap(pi, i);
				minSiftUp(pi);
			}
			
		}

	}
	
	private Node getMin() {
		
		Node nm = heap.get(0);
		swap(0, heap.size() - 1);
		heap.remove(heap.size() - 1);
		minSiftDown(0);
		
		return nm;
	}
	
	private void decreaseKey(int i, long d) {
		
//		System.out.println(" decreaseKey i: " + i + " d: " + d);
		Node dn = map.get(i);
		dn.dist = d;
		int dni = heap.indexOf(dn);
		minSiftUp(dni);
	}

    
    private void swap(int i, int n) {
    	
    	Node tn = heap.get(i);
    	heap.set(i, heap.get(n));
    	heap.set(n, tn);
    	
    }
	
	public long distance(int s, int t) {
		
		//implements Dijkstra's algorithm
		
		root = s;
		
		decreaseKey(s, 0);							//start the algorithm on this node
		
		while(!heap.isEmpty()) {
			Node r = getMin();
			r.visited = true;
			
//			System.out.println("processing Node: " + r.index);
			
			
			for(Edge e : graph.get(r.index)) {
				
				if(r.dist == Long.MAX_VALUE)		//this is an unreachable node if it's dist is infinite
					break;
				
				if(map.get(e.target).dist > r.dist + e.length) {
					decreaseKey(e.target, r.dist + e.length);
					map.get(e.target).pindex = r.index;				//min path is my daddy
				}
			}
		}
		
		
        return map.get(t).dist == Long.MAX_VALUE? -1 : map.get(t).dist;
    }

}




public class Dijkstra {
    

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        
        int n = scanner.nextInt();
        int m = scanner.nextInt();
        
        Graph g = new Graph(n);
        
        for (int i = 0; i < m; i++) {
            int x, y, w;
            x = scanner.nextInt() - 1;		//shift indexes to 0 based
            y = scanner.nextInt() - 1;
            w = scanner.nextInt();

            Edge e = new Edge(y, w);
            g.graph.get(x).add(e);
            
            
        }
        int x = scanner.nextInt() - 1;
        int y = scanner.nextInt() - 1;
        System.out.println(g.distance(x, y));
        scanner.close();
    }
}

