import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class ConnectedComponents {
	
	static int[] visited;
	static ArrayList<ArrayList<Integer>> adj;
	
	private static void explore(int x, int c) {
		visited[x] = c;
		for(int vx : adj.get(x)) {
			if(visited[vx] != c)
				explore(vx, c);
		}
		
	}
	
    private static int numberOfComponents() {
//        int result = 0;
        
    	visited = new int[adj.size()];
    	Arrays.fill(visited, 0);
    	
    	int cc = 0;
    	int i = 0;
        while( i < adj.size()) {
        	if(visited[i] == 0) {
           		++cc;
        		explore(i, cc);
        	}
        	++i;
        }
        return cc;
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
            adj.get(y - 1).add(x - 1);
        }
        System.out.println(numberOfComponents());
        scanner.close();
    }
}

