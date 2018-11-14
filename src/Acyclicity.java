import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class Acyclicity {
	
	static int n = 0;
	static int m = 0;
	static int[] pre;
	static int[] post;
	static int clock = 1;
	static ArrayList<ArrayList<Integer>> adj = null;
	
    private static int acyclic(ArrayList<ArrayList<Integer>> adj) {
        //write your code here
    	dfs();
    	int el;
    	int elPost = 0;
    	int result = 0;
    	
    	for(ArrayList<Integer> al : adj) {
    		el = adj.indexOf(al);
//    		System.out.println("el: " + el);
    		elPost = post[el];
    		for(int v : al) {
    			if(post[v] > elPost)
    				result = 1;
    		}
    	}

        return result;
    }
    
    private static void explore(int x) {
    	
    	if (pre[x] != 0)
    		return;
    	
		pre[x] = clock;
		++clock;
		
		ArrayList<Integer> al = adj.get(x);
		
		System.out.println("explore v->al: " + Arrays.toString(al.toArray()));
		
		for(int vx : al) {
			System.out.println("vx: " + vx + " pre[vx]: " + pre[vx]);
			if(pre[vx] == 0)
				explore(vx);
		}
		post[x] = clock;
		++clock;
		System.out.println("x: " + x + " pre: " + pre[x] + " post: " + post[x]);
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
        	System.out.println("Added v: " + i);
            adj.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            adj.get(x - 1).add(y - 1);	//shift all indexes to 0 based

        }
        
        for(ArrayList<Integer> al : adj) {
        	System.out.println("v: " + adj.indexOf(al) + " " + Arrays.toString(al.toArray()));
        }
        
        System.out.println(acyclic(adj));
        scanner.close();
    }
}

