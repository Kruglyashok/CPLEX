import java.util.Vector;
import ilog.concert.*;
import ilog.cplex.*;


public class LinProg {

	public static void main(String[] args) {
		//a test matrix
		int mat[][] = {	{0,40,-1,-1,180},
						{40,0,22,6,45},
						{-1,22,0,14,-1},
						{-1,6,14,0,20},
						{180,45,-1,20,0}};  
		leastPath(mat);
		doubleLeastPath(mat);
		int d[] = {1,1,1,1,1,1,1};
		System.out.println("Robust results:");
		for (int i =0; i < d.length; ++i) {
			System.out.println("Gamma = " + i + "  robust obj func = " + robust(mat,i,d));
			}
	}
	
	public static Vector<Integer> insideNums(int mat[][], int numRow, Vector<Integer> edges) {
		Vector<Integer> inside = new Vector<>();
		int counter = 0;
		while(mat[numRow][counter] != 0) {
			if (mat[numRow][counter] != -1) 
				inside.add(edges.indexOf(mat[numRow][counter]));
			counter++;
		}
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
			summ.clear();
			//count how many y's go in finish vertex
			//second constraint (for finish)
			for (int i = 0; i < mat.length; ++i) {
				if(mat[mat.length-1][i] > 0) {
					summ.addTerm(1.0, y.elementAt(edges.indexOf(mat[mat.length-1][i])));
				}
			}
			cplex.addEq(summ, 1.0);			
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
	
	public static Vector<IloIntVar> doubleLeastPath(int mat[][]) {
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
			
			//define vars
			for (int i =0; i < mat.length; ++i) {
				y.add(cplex.intVar(0,Integer.MAX_VALUE));
			}
			//create objective function
			IloLinearNumExpr objective = cplex.linearNumExpr();
			
			//expressions
			objective.addTerm(-1.0, y.elementAt(0));
			objective.addTerm(1.0, y.elementAt(y.size()-1));
			
			//define objective
			cplex.addMaximize(objective);
			
			//define constraints
			for (int i=0; i < mat.length; ++i) {
				for (int j =0; j < i; ++j) {
					if (mat[i][j] == -1 || mat[i][j] == 0) continue;
					cplex.addGe(cplex.sum(mat[i][j], 
							cplex.sum(cplex.prod(y.elementAt(i), -1.0),y.elementAt(j))), 0);
				}
			}
			//solve
			if(cplex.solve()) {
				System.out.println("obj value = " + cplex.getObjValue());
			}
			// write model to file
			 cplex.exportModel("lpex2.lp");	
		}
		catch(IloException e) {
			e.printStackTrace();
		}
		
		
		return y;		
	}	

	public static double robust(int mat[][], int gamma, int d[]) {
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
			//defining y
			Vector<IloIntVar> y = new Vector<>();
			for(int i=0; i < edges.size(); ++i) {
				y.add(cplex.intVar(0,1));
			}
			//definig mu
			IloIntVar mu = cplex.intVar(0, Integer.MAX_VALUE);
			//defining q
			Vector<IloIntVar> q = new Vector<>();
			for(int i=0; i < edges.size(); ++i) {
				q.add(cplex.intVar(0, Integer.MAX_VALUE));
			}
			//create objective function
			IloLinearNumExpr objective = cplex.linearNumExpr();
			//expressions
			for (int i = 0; i < edges.size(); ++i) {
				objective.addTerm(edges.elementAt(i),y.elementAt(i));
			}
			objective.addTerm(gamma, mu);
			for (int i=0; i < edges.size(); ++i) {
				objective.addTerm(1, q.elementAt(i));
			}
			System.out.println("Objective in 3d task is : " + objective.toString());
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
			summ.clear();
			//count how many y's go in finish vertex
			//second constraint (for finish)
			for (int i = 0; i < mat.length; ++i) {
				if(mat[mat.length-1][i] > 0) {
					summ.addTerm(1.0, y.elementAt(edges.indexOf(mat[mat.length-1][i])));
				}
			}
			cplex.addEq(summ, 1.0);			
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
				summ2.clear();
			}
			//constraint for undefined variables
			for (int i=0; i < edges.size(); ++i) {
				IloLinearNumExpr summ2 = cplex.linearNumExpr();
				summ2.addTerm(mu,1);
				summ2.addTerm(q.elementAt(i), 1);
				summ2.addTerm(y.elementAt(i), -1*d[i]);
				cplex.addGe(summ2, 0);		
				summ2.clear();
			}		
			//solve
			if(cplex.solve()) {
				System.out.println("Robust obj value = " + cplex.getObjValue());
				System.out.println("Solution:");
				for (int i =0; i < y.size(); ++i) {
					System.out.print(cplex.getValue(y.elementAt(i)) + "\t");
				}
				// write model to file
				 cplex.exportModel("lpex3.lp");
				 return cplex.getObjValue();
			}
 		}
		catch(IloException e) {
			e.printStackTrace();
		}
		
		return -1;
 	}
}
