import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;

public class StronglyConnected {
	
	static int n = 0;
	static int m = 0;
	static int[] pre;
	static int[] post;
	static int[] order;
	static int[] postIndex;
	static boolean[] linked;
	static int clock = 1;
	static ArrayList<ArrayList<Integer>> graph = null;
	static ArrayList<ArrayList<Integer>> graphR = null;
	private static Random random = new Random();
	
	
    private static int numberOfStronglyConnectedComponents(ArrayList<ArrayList<Integer>> adj) {
    	
    	dfs();

        order = Arrays.copyOf(post, n);
        
        randomizedQuickTopoSort(order, postIndex, 0, n - 1);
        
        //TODO: implement the CC count algorithm here
        
        
        ArrayList<Integer> result = new ArrayList<Integer>();
        
        for(int e = n - 1; e >= 0; --e) {
        	result.add(postIndex[e]);
        }
        return 0;
    }
    
    private static void explore(int x) {
    	
    	if (pre[x] > 0)					//allow 2nd visit for singleton vertexes with negative pre values
    		return;
    	
		pre[x] = clock;
		++clock;
		
		ArrayList<Integer> al = graph.get(x);
		
		if(!al.isEmpty())
			linked[x] = true;
		
		if(al.isEmpty() && linked[x] == false)	//if this vertex is linked don't allow more visits
				pre[x] *= -1;					//This vertex might be connected and included in a 2nd visit
		
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
    
    private static void dfs() {
    	for(int i = 0; i < n; ++i) {
    		pre[i] = 0;
    		post[i] = 0;
    	}
    	
    	for(int j = 0; j < n; ++j) {
    		if(pre[j] == 0)
    			explore(j);
    	}
    		
    }
    
    private static void randomizedQuickTopoSort(int[] a, int[] ai, int l, int r) {
        if (l >= r) {
            return;
        }
        int k = random.nextInt(r - l + 1) + l;
        int t = a[l];
        a[l] = a[k];
        a[k] = t;
        
        t = ai[l];			//mirror the operations on the index array
        ai[l] = ai[k];
        ai[k] = t;
        
        //use partition2
        
        int m = partition2(a, ai, l, r);
        randomizedQuickTopoSort(a, ai, l, m - 1);
        randomizedQuickTopoSort(a, ai, m + 1, r);
        

    }
    
    private static int partition2(int[] a, int[] ai, int l, int r) {
        int x = a[l];
        int j = l;
        for (int i = l + 1; i <= r; i++) {
            if (a[i] <= x) {
                j++;
                int t = a[i];
                a[i] = a[j];
                a[j] = t;
                
                t = ai[i];				//mirror the operations on the index array
                ai[i] = ai[j];
                ai[j] = t;
            }
        }
        int t = a[l];
        a[l] = a[j];
        a[j] = t;
        
        t = ai[l];				//mirror the operations on the index array
        ai[l] = ai[j];
        ai[j] = t;
        return j;
    }

    public static void main(String[] args) {
    	Scanner scanner = new Scanner(System.in);
    	
        n = scanner.nextInt();
        m = scanner.nextInt();
        pre = new int[n];
        post = new int[n];
        postIndex = new int[n];
        linked = new boolean[n];
        Arrays.fill(linked, false);
        graph = new ArrayList<ArrayList<Integer>>();
        graphR = new ArrayList<ArrayList<Integer>>();
        
        for (int i = 0; i < n; i++) {
//        	System.out.println("Added v: " + i);
            graph.add(new ArrayList<Integer>());
            graphR.add(new ArrayList<Integer>());
            postIndex[i] = i;
        }
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            graph.get(x - 1).add(y - 1);	//shift all indexes to 0 based
            graphR.get(y - 1).add(x -1);	//reversed graph to find sinks

        }
        System.out.println(numberOfStronglyConnectedComponents(graph));
    }
}

