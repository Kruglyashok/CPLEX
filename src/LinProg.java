import java.util.Vector;
import ilog.concert.*;
import ilog.cplex.*;


public class LinProg {

	public static void main(String[] args) {
		//model1();
		//a test matrix
		int mat[][] = {	{0,40,-1,-1,180},
						{40,0,22,6,45},
						{-1,22,0,14,-1},
						{-1,6,14,0,20},
						{180,45,-1,20,0}};  
		leastPath(mat);
	}
	
	public static Vector<Integer> insideNums(int mat[][], int numRow, Vector<Integer> edges) {
		Vector<Integer> inside = new Vector<>();
		int counter = 0;
		while(mat[numRow][counter] != 0) {
			if (mat[numRow][counter] != -1) 
				inside.add(edges.indexOf(mat[numRow][counter]));
			counter++;
		}
		
		System.out.println(inside.toString());
		return inside;
	}
	public static Vector<Integer> outsideNums(int mat[][], int numRow, Vector<Integer> edges) {
		Vector<Integer> outside = new Vector<>();
		int counter = mat.length-1;
		while(mat[numRow][counter] != 0) {
			if (mat[numRow][counter] != -1) 
				outside.add(edges.indexOf(mat[numRow][counter]));
			counter--;
		}
		
		System.out.println(outside.toString());
		return outside;
	}
	
	public static Vector<IloIntVar> leastPath(int mat[][]) {
		Vector<IloIntVar> y = new Vector<IloIntVar>();
		try {
			//create model
			IloCplex cplex = new IloCplex();
			//count edges
			Vector<Integer> edges = new Vector<>();
			for (int i=1; i < mat.length; ++i) {
				for (int j = 0; j < i; ++j) {
					if (mat[i][j] > 0) {
					edges.add(mat[i][j]);
					}
				}
			}
			System.out.println(edges.toString());
			System.out.println("edges = " + edges.size());
			
			//define vars
			for (int i =0; i < edges.size(); ++i) {
				y.add(cplex.intVar(0,Integer.MAX_VALUE));
			}
			//create objective function
			IloLinearNumExpr objective = cplex.linearNumExpr();
			
			//expressions
			for (int i = 0; i < edges.size(); ++i) {
				objective.addTerm(edges.elementAt(i),y.elementAt(i));
			}
			System.out.println(objective.toString());
			//define objective
			cplex.addMinimize(objective);
			
			//define constraints
			
			//count how many y's go out of 0 vertex
			//first constraint (for start)
			IloLinearNumExpr summ = cplex.linearNumExpr();
			for (int i =0; i < mat.length; ++i) {
				if(mat[0][i] > 0) {
					summ.addTerm(1.0, y.elementAt(edges.indexOf(mat[0][i])));
				}
			}
			cplex.addEq(summ, 1.0);
			System.out.println(summ.toString());
			summ.clear();
			//count how many y's go in finish vertex
			//second constraint (for finish)
			for (int i = 0; i < mat.length; ++i) {
				if(mat[mat.length-1][i] > 0) {
					summ.addTerm(1.0, y.elementAt(edges.indexOf(mat[mat.length-1][i])));
				}
			}
			cplex.addEq(summ, 1.0);
			System.out.println(summ.toString());
			
			//third constraint
			Vector<Integer> inside = new  Vector<>();
			Vector<Integer> outside = new  Vector<>();
			
			for (int i = 1; i < mat.length - 1; ++i) {
				IloLinearNumExpr summ2 = cplex.linearNumExpr();
				inside = insideNums(mat,i,edges);
				outside = outsideNums(mat,i,edges);
				//here we need for rows in mat 1..mat.length-2 to create terms
				for (int j =0; j < inside.size(); ++j) {
					summ2.addTerm(1.0, y.elementAt(inside.elementAt(j)));
				}
				for (int j =0; j < outside.size(); ++j) {
					summ2.addTerm(-1.0, y.elementAt(outside.elementAt(j)));
				}
				cplex.addEq(summ2,0.0);
				System.out.println("summ2 : " + summ2.toString());
				summ2.clear();
			}			
			//solve
			if(cplex.solve()) {
				System.out.println("obj value = " + cplex.getObjValue());
				System.out.println("Solution:");
				for (int i =0; i < y.size(); ++i) {
					System.out.print(cplex.getValue(y.elementAt(i)) + "\t");
				}
			}
			// write model to file
			 cplex.exportModel("lpex1.lp");	
		}
		catch(IloException e) {
			e.printStackTrace();
		}
		
		return y;
	}
	
	/*
	public static void model1() {
		try {
			//create model
			IloCplex cplex = new IloCplex();
			
			//define vars
			IloIntVar x = cplex.intVar(0, Integer.MAX_VALUE, "x");
			IloIntVar y = cplex.intVar(0, Integer.MAX_VALUE, "y");
			
			//create objective func
			IloLinearNumExpr objective = cplex.linearNumExpr();
			//IloLinearNumExpr objective = cplex.linearNumExpr();
			
			//expressions	
			objective.addTerm(0.16, x);
			objective.addTerm(0.15, y);
			
			//define objective
			cplex.addMinimize(objective);
			
			//define constraints 
			cplex.addGe(cplex.sum(cplex.prod(61, x), cplex.prod(60, y)),300);
			cplex.addGe(cplex.sum(cplex.prod(12, x), cplex.prod(6, y)), 43);
			cplex.addGe(cplex.sum(cplex.prod(13, x), cplex.prod(30, y)), 111);
			
			//solve
			if(cplex.solve()) {
				System.out.println("obj value = " + cplex.getObjValue());
				System.out.println("x = " + cplex.getValue(x));
				System.out.println("y = " + cplex.getValue(y));
				
				
			}
			
			
		}
		catch(IloException e) {
			e.printStackTrace();
		}
		
	}
	*/
}