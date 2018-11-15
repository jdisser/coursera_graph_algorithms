import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

public class Toposort {
	
	static int n = 0;
	static int m = 0;
	static int[] pre;
	static int[] post;
	static int[] order;
	static int[] postIndex;
	static int clock = 1;
	static ArrayList<ArrayList<Integer>> adj = null;
	private static Random random = new Random();
	
	
    private static ArrayList<Integer> toposort(ArrayList<ArrayList<Integer>> adj) {

    	dfs();

        order = Arrays.copyOf(post, n);
        
        randomizedQuickTopoSort(order, postIndex, 0, n - 1);
        
        ArrayList<Integer> result = new ArrayList<Integer>();
        
        for(int e = n - 1; e >= 0; --e) {
        	result.add(postIndex[e]);
        }

        return result;
    }

	 private static void explore(int x) {
	    	
	    	if (pre[x] != 0)
	    		return;
	    	
			pre[x] = clock;
			++clock;
			
			ArrayList<Integer> al = adj.get(x);
			
	//		System.out.println("explore v->al: " + Arrays.toString(al.toArray()));
			
			for(int vx : al) {
	//			System.out.println("vx: " + vx + " pre[vx]: " + pre[vx]);
				if(pre[vx] == 0)
					explore(vx);
			}
			post[x] = clock;
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
        adj = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < n; i++) {
//        	System.out.println("Added v: " + i);
            adj.add(new ArrayList<Integer>());
            postIndex[i] = i;
        }
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            adj.get(x - 1).add(y - 1);	//shift all indexes to 0 based

        }
        ArrayList<Integer> order = toposort(adj);
        for (int x : order) {
            System.out.print((x + 1) + " ");
        }
        scanner.close();
    }
}

