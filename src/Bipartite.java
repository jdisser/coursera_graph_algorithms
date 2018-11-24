import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Scanner;

class Graph {
	
	public ArrayList<ArrayList<Integer>> graph;
	public Queue<Integer> queue;
	public int[] dist;
	public int[] parent;
	public int[] color;
	public int[] orderIndex;
	public int[] connected;
	public boolean[] linked;
	public int clock;
	public int n;
	public int cc;

	
	public Graph(int n) {
		this.n = n;
		this.dist = new int[n];
		this.parent = new int[n];
		this.orderIndex = new int[n];
		this.linked = new boolean[n];
		this.connected = new int[n];
		
        Arrays.fill(this.linked, false);
        this.graph = new ArrayList<ArrayList<Integer>>();
        
        this.queue = new LinkedList<Integer>();
        
        for (int i = 0; i < n; i++) {
//        	System.out.println("Added v: " + i);
            this.graph.add(new ArrayList<Integer>());
            this.dist[i] = n+1;								//n+1 -> infinity
            this.parent[i] = 0;
        }
        
        this.clock = 1;
        this.cc = 0;
	}
	
	public boolean bfs(int s) {
		dist[s] = 0;
		parent[s] = -1;			//-1 -> no parent
		queue.add(s);
		boolean result = true;
		
		while(!queue.isEmpty()) {
			int u = queue.poll();
			for(int e : this.graph.get(u)) {
				if(dist[e] > n) {
					dist[e] = dist[u] + 1;
					parent[e] = u;
					queue.add(e);
				} else {
					if(dist[e] % 2 == dist[u] % 2)		//test for bipartite graph
						result = false;
				}
			}
		}
		return result;
	}
	
	public int distance(int s, int t) {
        //write your code here
		bfs(s);
        return dist[t] > n ? -1 : dist[t];
    }
	
	public int bipartite() {
        //write your code here
        return bfs(1)? 1 : 0;
    }
}




public class Bipartite {
    

    public static void main(String[] args) {


    	Scanner scanner = new Scanner(System.in);
    	
    	int n;
    	int m;
    	
    	
        n = scanner.nextInt();
        m = scanner.nextInt();

        Graph g = new Graph(n);

        
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            g.graph.get(x - 1).add(y - 1);	//shift all indexes to 0 based
            g.graph.get(y - 1).add(x -1);	//undirected graph

        }

        System.out.println(g.bipartite());
        
        scanner.close();
    }
}

