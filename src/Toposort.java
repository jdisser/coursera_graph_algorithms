import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class Toposort {
	
	static int n = 0;
	static int m = 0;
	static int[] pre;
	static int[] post;
	static int clock = 1;
	static ArrayList<ArrayList<Integer>> adj = null;
	
	
    private static ArrayList<Integer> toposort(ArrayList<ArrayList<Integer>> adj) {
        
    	dfs();
    	
        ArrayList<Integer> order = new ArrayList<Integer>();
        
        for(int i = n - 1; n >= 0; --i) {
        	order.add(post[i]);
        }

        return order;
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

    public static void main(String[] args) {
    	Scanner scanner = new Scanner(System.in);
        n = scanner.nextInt();
        m = scanner.nextInt();
        pre = new int[n];
        post = new int[n];
        adj = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < n; i++) {
//        	System.out.println("Added v: " + i);
            adj.add(new ArrayList<Integer>());
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

