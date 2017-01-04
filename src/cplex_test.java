/**
 * Created by loujian on 11/28/16.
 */
/*
Maximize x1+2x2+3x3

subject to
-x1+x2+x3 <= 20
x1 - 3x2 + x3 <= 30

with these bounds
0 <= x1 <= 40
0 <= x2 <= +INF
0 <= x3 <= +INF
*/



import ilog.concert.*;
import ilog.cplex.*;

class cplex_test {

    static public void main(String[] args) {
        try {
            IloCplex cplex = new IloCplex();
            // create model and solve it
            double[] lb = {0.0, 0.0, 0.0};
            double[] ub = {40.0, Double.MAX_VALUE, Double.MAX_VALUE};
            IloNumVar[] x = cplex.numVarArray(3, lb, ub);

            double[] objvals = {1.0, 2.0, 3.0};

            cplex.addMaximize(cplex.scalProd(x, objvals));  //scalProd():Creates and returns a linear expression representing the scalar product of the given variables.

            cplex.addLe(cplex.sum(cplex.prod(-1.0, x[0]),
                    cplex.prod( 1.0, x[1]),
                    cplex.prod( 1.0, x[2])), 20.0);
            cplex.addLe(cplex.sum(cplex.prod( 1.0, x[0]),
                    cplex.prod(-3.0, x[1]),
                    cplex.prod( 1.0, x[2])), 30.0);

            if ( cplex.solve() ) {
                System.out.println("Solution status = " + cplex.getStatus());
                System.out.println("Solution value = " + cplex.getObjValue());
                double[] val = cplex.getValues(x);
                int ncols = cplex.getNcols();
                for (int j = 0; j < ncols; ++j)
                    System.out.println("Column: " + j + " Value = " + val[j]);
            }




            cplex.end();







        } catch (IloException e) {
            System.err.println("Concert exception caught: " + e);
        }





    }


}
