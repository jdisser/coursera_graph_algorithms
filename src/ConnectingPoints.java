import java.util.ArrayList;
import java.util.Scanner;
import static java.lang.Math.pow;

class Node {
	public int index;
	public double x,y;
	public int pindex;
	
	public Node (int x, int y) {
		this.x = (double) x;
		this.y = (double) y;
		this.index = -1;
		this.pindex = -1;
	}
}

class Edge {
	public Node v1;
	public Node v2;
	public double length;
	
	public Edge(Node v1, Node v2) {
		this.v1 = v1;
		this.v2 = v2;
		this.length = setLength();
	}
	
	private double setLength() {
		//lengths are stored as sum of squares and must be converted to length
		return pow((v1.x - v2.x), 2.0) + pow((v1.y - v2.y) , 2.0);
	}
}

class Graph {
	
	public ArrayList<Node> graph;
	public ArrayList<Edge> eQueue;
	public ArrayList<Edge> xMst;
	
	public Graph () {
		this.graph = new ArrayList<Node>();
		this.eQueue = new ArrayList<Edge>();
		this.xMst = new ArrayList<Edge>();
	}
	
	public int addNode(Node v) {
		
		int vi;
		
		graph.add(v);
		vi = graph.indexOf(v);
		graph.get(vi).index = vi;
		
		return vi;	
	}
	
	
}





public class ConnectingPoints {
    private static double minimumDistance(int[] x, int[] y) {
        double result = 0.;
        //write your code here
        return result;
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int n = scanner.nextInt();

        Graph V = new Graph();
        
        int x,y;
        
        
        for (int i = 0; i < n; i++) {
            x = scanner.nextInt();
            y = scanner.nextInt();
            
            Node v = new Node(x, y);
            V.addNode(v);
        }
        System.out.println(minimumDistance(x, y));
    }
}

