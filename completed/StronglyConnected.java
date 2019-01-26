import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;


class Graph {
	
	public ArrayList<ArrayList<Integer>> graph;
	public int[] pre;
	public int[] post;
	public int[] order;
	public int[] orderIndex;
	public int[] connected;
	public boolean[] linked;
	public int clock;
	public int n;
	public int cc;
	private Random random = new Random();
	
	public Graph(int n) {
		this.n = n;
		this.pre = new int[n];
		this.post = new int[n];
		this.orderIndex = new int[n];
		this.linked = new boolean[n];
		this.connected = new int[n];
        Arrays.fill(this.linked, false);
        this.graph = new ArrayList<ArrayList<Integer>>();
        
        for (int i = 0; i < n; i++) {
//        	System.out.println("Added v: " + i);
            this.graph.add(new ArrayList<Integer>());
            this.orderIndex[i] = i;
            this.connected[i] = 0;
        }
        
        this.clock = 1;
        this.cc = 0;
	}
	
	public void explore(int x) {
    	
    	if (pre[x] > 0)					//allow 2nd visit for singleton vertexes with negative pre values
    		return;
    	
    	pre[x] = clock;
		++clock;
		
		ArrayList<Integer> al = graph.get(x);
		
		if(!al.isEmpty())
			linked[x] = true;
		
		if(al.isEmpty() && linked[x] == false)	//if this vertex is linked don't allow more visits
			this.pre[x] *= -1;					//This vertex might be connected and included in a 2nd visit
		
//		System.out.println("explore x: " + x + " v->al: " + Arrays.toString(al.toArray()));
		
		for(int vx : al) {
			linked[vx] = true;
//			System.out.println("vx: " + vx + " pre[vx]: " + pre[vx]);
			if(pre[vx] <= 0)
				explore(vx);
		}
		post[x] = clock;				//the pre and post number are overwritten on the 2nd visit
		++clock;
//		System.out.println("x: " + x + " pre: " + pre[x] + " post: " + post[x]);
	}
	
	public void exploreConnected(int x, int ci) {
    	
//		System.out.println("exploreConnected x: " + x + " ci: " + ci + " connected[x]: " + connected[x]);
		
    	if (connected[x] != 0)					
    		return;
    	
    	connected[x] = ci;
		
		ArrayList<Integer> al = graph.get(x);
		
		if(al.isEmpty())
			return;
		
		for(int vx : al) {

			if(connected[vx] == 0)
				exploreConnected(vx, ci);
		}

	}
	
	public void dfs() {
    	for(int i = 0; i < n; ++i) {
    		pre[i] = 0;
    		post[i] = 0;
    	}
    	
    	for(int j = 0; j < n; ++j) {
    		if(pre[j] == 0)
    			explore(j);
    	}
    		
    }
	
	public void ccs() {
		
//		System.out.println("Initial in ccs cc: " + cc + " connected[]: " + Arrays.toString(connected));
		
		
		for (int i = 0; i < n; ++i) {
			if(connected[orderIndex[i]] == 0) {
				++cc;
				exploreConnected(orderIndex[i], cc);
//				System.out.println("cc: " + cc + " connected[]: " + Arrays.toString(connected));
				
			}
		}
	}
	
	public void toposort() {
		order = Arrays.copyOf(post, n);
		randomizedQuickTopoSort(0, n - 1);
	}
	
	private void randomizedQuickTopoSort(int l, int r) {
        if (l >= r) {
            return;
        }
        int k = random.nextInt(r - l + 1) + l;
        int t = order[l];
        order[l] = order[k];
        order[k] = t;
        
        t = orderIndex[l];			//mirror the operations on the index array
        orderIndex[l] = orderIndex[k];
        orderIndex[k] = t;
        
        //use partition2
        
        int m = partition2(l, r);
        randomizedQuickTopoSort(l, m - 1);
        randomizedQuickTopoSort(m + 1, r);
        

    }
	
	private int partition2(int l, int r) {
        int x = order[l];
        int j = l;
        for (int i = l + 1; i <= r; i++) {
            if (order[i] <= x) {
                j++;
                int t = order[i];
                order[i] = order[j];
                order[j] = t;
                
                t = orderIndex[i];				//mirror the operations on the index array
                orderIndex[i] = orderIndex[j];
                orderIndex[j] = t;
            }
        }
        int t = order[l];
        order[l] = order[j];
        order[j] = t;
        
        t = orderIndex[l];				//mirror the operations on the index array
        orderIndex[l] = orderIndex[j];
        orderIndex[j] = t;
        return j;
    }
	
	
}


public class StronglyConnected {
	

	
    private static int numberOfStronglyConnectedComponents(Graph g, Graph gr) {
    	
    	gr.dfs();
        gr.toposort();
        
//        for(int i = 0; i < g.n; ++i) {
//        	System.out.println("v: " + i + " " + Arrays.toString(g.graph.get(i).toArray()));
//        }

        
        for(int e = 0; e <  gr.n; ++e) {
        	g.orderIndex[gr.n - 1 - e] = gr.orderIndex[e];
        }
       
//        System.out.println("Reversed gr orderIndex: " + Arrays.toString(g.orderIndex));
        
        g.ccs();
      
        return g.cc;
    }
    
    
    
    
    
    
    
    

    public static void main(String[] args) {
    	Scanner scanner = new Scanner(System.in);
    	
    	int n;
    	int m;
    	
    	
        n = scanner.nextInt();
        m = scanner.nextInt();

        Graph g = new Graph(n);
        Graph gr = new Graph(n);
        
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            g.graph.get(x - 1).add(y - 1);	//shift all indexes to 0 based
            gr.graph.get(y - 1).add(x -1);	//reversed graph to find sinks

        }
        System.out.println(numberOfStronglyConnectedComponents(g, gr));
        scanner.close();
    }
}

