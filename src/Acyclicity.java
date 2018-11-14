import java.util.ArrayList;
import java.util.Scanner;

public class Acyclicity {
	
	static int n = 0;
	static int m = 0;
	static int[] pre = new int[n];
	static int[] post = new int[n];
	static int clock = 0;
	static ArrayList<ArrayList<Integer>> adj = null;
	
    private static int acyclic(ArrayList<ArrayList<Integer>> adj) {
        //write your code here
    	dfs();
    	int el;
    	int elPost = 0;
    	int result = 0;
    	
    	for(ArrayList<Integer> al : adj) {
    		el = adj.indexOf(al);
    		elPost = post[el];
    		for(int v : al) {
    			if(post[v] > elPost)
    				result = 1;
    		}
    	}

        return result;
    }
    
    private static void explore(int x) {
		pre[x] = clock;
		++clock;
		for(int vx : adj.get(x)) {
			if(pre[vx] != 0)
				explore(vx);
		}
		post[x] = clock;
		++clock;	
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
        int n = scanner.nextInt();
        int m = scanner.nextInt();
        adj = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < n; i++) {
            adj.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            adj.get(x - 1).add(y - 1);

        }
        
        System.out.println(acyclic(adj));
        scanner.close();
    }
}

