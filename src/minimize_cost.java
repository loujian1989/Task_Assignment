/**
 * Created by loujian on 11/28/16.
 */

import ilog.concert.*;
import ilog.cplex.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;



public class minimize_cost {

    int T, R, N, K;
    double[][] resource_need;  //R_{t,r}
    double[][] lost_for_resource;  //l_{i,r}
    double[][] upper_bound_player;

    minimize_cost(int T, int R, int N, int K, double[][] resource, double[][] lost, double[][] upper_bound)
    {
        this.T= T;
        this.R= R;
        this.N= N;
        this.K= K;

        resource_need= new double[T][R];
        lost_for_resource= new double[N][R];
        upper_bound_player= new double[N][R];


        for(int t=0; t < T; t++)
            for(int r=0; r < R; r++)
                resource_need[t][r]= resource[t][r];

        for(int i=0; i < N; i++)
            for(int r=0; r<R; r++)
                lost_for_resource[i][r]= lost[i][r];

        for(int i=0; i < N; i++)
            for(int r=0; r < R; r++)
                upper_bound_player[i][r] = upper_bound[i][r];
    }


    public double[] solve_problem() {

        try {
            IloCplex cplex = new IloCplex();
            // create model and solve it

            //lower bound and upper bound of x
            double[] lb_x = new double[(N)*(R)];
            double[] ub_x = new double[(N)*(R)];
            for(int i=0; i<N; i++)
                for(int r=0; r<R; r++)
                {
                    lb_x[i*R + r]= 0;
                    ub_x[i*R + r]= upper_bound_player[i][r];
                }
            IloNumVar[] x= cplex.numVarArray(N*R, lb_x, ub_x);

            //lower bound and upper bound of y
            double[] lb_y= new double[(T+1)*(N+1)*(R+1)];
            double[] ub_y= new double[(T+1)*(N+1)*(R+1)];

            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                    {
                        lb_y[t*(N*R)+ i*R + r] = 0;
                        ub_y[t*(N*R)+ i*R + r] = Double.MAX_VALUE;
                    }
            IloNumVar[] y= cplex.numVarArray(T*N*R, lb_y, lb_y);

            //p is binary variables
            IloIntVar[]p= cplex.boolVarArray(T*N);

            double[] objvals= new double[T*N*R];
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                    {
                        objvals[ t*(N*R) + i*R + r ] = lost_for_resource[i][r];
                    }

            //We would like to minimize the overall cost
            cplex.addMinimize(cplex.scalProd(y, objvals));

            //add constraint: \sum_{t\in T} p_{t,i} \leq 1, \forall i\in I
            for(int i=0; i<N; i++)
            {
                IloIntVar[] p_i= cplex.boolVarArray(T);
                double[] constant_1= new double[T];
                for(int t=0; t<T; t++)
                {
                    constant_1[t]=1;
                    p_i[t] = p[t * N + i];
                }
                cplex.addLe( cplex.scalProd(p_i, constant_1), 1.0);
            }

            //add constraint: \sum_{i\in I} p_{t, i}\leq K, \forall t\in T
            for(int t=0; t<T; t++)
            {
                IloIntVar[] p_t= cplex.boolVarArray(N);
                double[] constant_1= new double[N];
                for(int i=0; i<N; i++)
                {
                    constant_1[i]=1;
                    p_t[i] = p[t* N + i];
                }
                cplex.addLe( cplex.scalProd(p_t, constant_1), K );
            }

            //add constraint: \sum_{i\in I} y_{t, i, r} \geq R_{t, r}, \forall, t\in T, r\in R
            for(int t=0; t<T; t++)
                for(int r=0; r<R; r++)
                {
                    IloNumVar[] y_tr= cplex.numVarArray(N, 0, Double.MAX_VALUE);
                    double[] constant_2= new double[N];
                    for(int i=0; i<N; i++)
                    {
                        constant_2[i]=1;
                        y_tr[i]= y[t*N*R+ i*R + r];
                    }
                    cplex.addGe(cplex.scalProd(y_tr, constant_2), resource_need[t][r]);
                }
            //add constraint: y_{t, i, r} - U_{i, r} p_{t, i} \leq 0, \forall t\in T, i\in I, r\in R
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                        cplex.addLe(cplex.sum( y[t*(N*R)+ i*R+ r], cplex.prod( -upper_bound_player[i][r], p[t*N + i] )), 0);

            //add constraint: y_{t, i, r}- x_{i, r} \leq 0
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                        cplex.addLe( cplex.sum( y[t*(N*R)+ i*R+ r], cplex.negative(x[i*R + r])), 0 );

            //add constraint: x_{i, r} - y_{t, i, r} + U_{i, r} p_{t, i} \leq U_{i, r}
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                        cplex.addLe( cplex.sum( x[i*R+ r], cplex.negative(y[t*(N*R) + i*R + r]), cplex.prod( upper_bound_player[i][r], p[t*N+ i] ) ) , upper_bound_player[i][r]);


            if ( cplex.solve() ) {
                if ( cplex.solve() ) {
                    System.out.println("Solution status = " + cplex.getStatus());
                    System.out.println("Solution value = " + cplex.getObjValue());
                    double[] x_val = cplex.getValues(x);
                    double[] p_val= cplex.getValues(p);
                    for (int i = 0; i < N ; i++)
                        for(int r=0; r<R ; r++)
                            System.out.println("Player "+ i + " should produce" + r + "Resource: " + x_val[i*R+r] + "\r\n");
                    for(int t=0; t<T; t++)
                    {
                        System.out.println("The task " + t + "consists the following players: \r\n");
                        for(int i=0; i<N; i++)
                        {
                            if(p_val[t*N+ i]==1.0)
                                System.out.println(i+ "\r\n");
                        }
                    }
                }


            }


            cplex.end();

        } catch (IloException e) {
            System.err.println("Concert exception caught: " + e);
        }

        return null;
    }




}
