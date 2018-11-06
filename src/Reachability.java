import java.util.ArrayList;
import java.util.Scanner;

public class Reachability {
	
	static boolean[] visited;
	static ArrayList<ArrayList<Integer>> adj;
	
	private static void explore(int x) {
		visited[x] = true;
		for(int vx : adj.get(x)) {
			if(!visited[vx])
				explore(vx);
		}
		
	}
	
    private static int reach(int x, int y) {
    	
    	visited = new boolean[adj.size()];
    	explore(x);
    	if(visited[y])
    		return 1;
    	else
    		return 0;
    }


    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int n = scanner.nextInt();
        int m = scanner.nextInt();
        adj = new ArrayList<ArrayList<Integer>>();
//        for (int i = 0; i < n; i++) {
//            adj.get(i) = new ArrayList<Integer>();
//        }
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            adj.get(x -1).add(y - 1);
            adj.get(y - 1).add(x - 1);
        }
        
        int x = scanner.nextInt() - 1;
        int y = scanner.nextInt() - 1;
        System.out.println(reach(x, y));
        scanner.close();
    }
}

