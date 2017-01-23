/**
 * Created by loujian on 1/20/17.
 For here, we assume the x, y, p,
 */



import ilog.concert.*;
import ilog.cplex.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Scanner;

public class Task {

    public static void main(String[] args) throws Exception
    {
        //Here we read the data from a file, we can revise it later
        Scanner cin=new Scanner(new File("tmp_test"));
        File writename = new File("test_out");
        writename.createNewFile();
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));


        int T= cin.nextInt(); //number of task
        int R= cin.nextInt(); //number of resource types
        int N= cin.nextInt(); //number of players
        int K= cin.nextInt(); //the maximum size of a team

        double[][] resource_need= new double[T][R];  //R_{t,r}
        double[][] lost_for_resource= new double[N][R];  //l_{i,r}
        double[][] utility = new double[N][N]; //u_i(j);
        double[] Omega = new double[N]; //Omega_i, denote the upper bound of epsilon
        double[] alpha = new double[N]; //alpha_i, denote the weight of epsilon for each player

        for(int i=0; i<N; i++)
            Omega[i]=10000; //here we just set Omega to some value

        //input the data from the file
        for(int t=0; t<T; t++)
            for(int r=0; r<R; r++)
                resource_need[t][r]= cin.nextDouble();

        for(int i=0; i<N; i++)
            for(int r=0; r<R; r++)
                lost_for_resource[i][r]= cin.nextDouble();

        for(int i=0; i<N; i++)
            for(int j=0; j<N; j++)
                utility[i][j]= cin.nextDouble();

        for(int i=0; i<N; i++)
            alpha[i]= cin.nextDouble();

        //For now we can start the optimization

        try {
            IloCplex task_solver = new IloCplex();

            IloRange[][]  rng = new IloRange[15][]; // Here we denote 30 kinds of constraints, we can define all these kinds of constraints later


            IloIntVar[] x= task_solver.boolVarArray(N*R);
            IloNumVar[] y= task_solver.boolVarArray(T*N*R);

            //lower bound and upper bound of epsilon
            double[] lb_epsilon= new double[N];
            double[] ub_epsilon= new double[N];

            for(int i=0; i<N; i++)
            {
                lb_epsilon[i]=0;
                ub_epsilon[i]=Omega[i];
            }
            //define the variable epsilon
            IloNumVar[] epsilon = task_solver.numVarArray(N, lb_epsilon, ub_epsilon);


            //lower bound and upper bound of delta
            double[] lb_delta= new double[T*N];
            double[] ub_delta= new double[T*N];

            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                {
                    lb_delta[t*N + i]= 0;
                    ub_delta[t*N + i]=Double.MAX_VALUE;
                }
            //define the variable delta
            IloNumVar[] delta = task_solver.numVarArray(T*N, lb_delta, ub_delta);

            //p is binary variables
            IloIntVar[] p= task_solver.boolVarArray(T*N);

            //z is binary variables
            IloIntVar[] z = task_solver.boolVarArray(T*T*N*N);


            double[] objvals= new double[T*N*R]; //the coefficient of y
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                    {
                        objvals[ t*(N*R) + i*R + r ] = lost_for_resource[i][r];
                    }

            //We would like to minimize the overall cost and the
            task_solver.addMinimize(task_solver.sum(task_solver.scalProd(y, objvals), task_solver.scalProd(alpha, epsilon)));


            //We would like to add constraints

            //add constraint 0: \sum_{t\in T} p_{t,i} \leq 1, \forall i\in I
            rng[0] = new IloRange[N];
            for(int i=0; i<N; i++)
            {
                int[] binary_vector = new int[T*N];
                for(int j=0; j<T*N; j++)
                    binary_vector[j]=0;

                for(int t=0; t<T; t++)
                    binary_vector[t*N+i]=1; //It is a way to denote the left equation of the constraint

                rng[0][i]= task_solver.addLe( task_solver.scalProd(p, binary_vector), 1);
            }

            //add constraint 1: \sum_{i\in I} p_{t, i}\leq K, \forall t\in T
            rng[1]= new IloRange[T];
            for(int t=0; t<T; t++)
            {
                int[] binary_vector = new int[T*N];
                for(int j=0; j<T*N; j++)
                    binary_vector[j]=0;

                for(int i=0; i<N; i++)
                    binary_vector[t*N+i]=1; //It is a way to denote the left equation of the constraint

                rng[1][t] = task_solver.addLe( task_solver.scalProd(p, binary_vector), K );
            }

            //add constraint 2: \sum_{i\in I} y_{t, i, r} \geq R_{t, r}, \forall, t\in T, r\in R
            rng[2]= new IloRange[T*R];
            for(int t=0; t<T; t++)
                for(int r=0; r<R; r++)
                {
                    int[] binary_vector = new int[T*N*R];
                    for(int j=0; j<T*N*R; j++)
                        binary_vector[j]=0;

                    for(int i=0; i<N; i++)
                        binary_vector[t*(N*R) + i*R + r]= 1;

                    rng[2][t*R+r] = task_solver.addGe(task_solver.scalProd(y, binary_vector), resource_need[t][r]);
                }

            //Here we miss constraint 3, which will be added in the constraint generation
            rng[3]= new IloRange[T*T*N];

            //add constraint 4 : y_{t, i, r} - p_{t, i} \leq 0, \forall t\in T, i\in I, r\in R
            rng[4] = new IloRange[T*N*R];
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                        rng[4][t*(N*R)+ i*R+ r] = task_solver.addLe(task_solver.sum( y[t*(N*R)+ i*R+ r], task_solver.negative(p[t*N + i])), 0);

            //add constraint 5: y_{t, i, r}- x_{i, r} \leq 0
            rng[5] = new IloRange[T*N*R];
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                        rng[5][t*(N*R)+ i*R+ r] = task_solver.addLe( task_solver.sum( y[t*(N*R)+ i*R+ r], task_solver.negative(x[i*R + r])), 0 );

            //add constraint 6: x_{i, r} - y_{t, i, r} + p_{t, i} \leq U_{i, r}
            rng[6]= new IloRange[T*N*R];
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    for(int r=0; r<R; r++)
                        rng[6][t*(N*R)+ i*R+ r] = task_solver.addLe( task_solver.sum( x[i*R+ r], task_solver.negative(y[t*(N*R) + i*R + r]), p[t*N+ i]  ) , 1);
/*
            //add constraint 7: z_{t, t', i, i'}- p_{t, i} \leq 0
            rng[7]= new IloRange[T*T*N*N];
            for(int t1=0; t1<T; t1++)
                for(int t2=0; t2<T; t2++)
                    for(int i1=0; i1< N; i1++)
                        for(int i2=0; i2< N; i2++)
                            rng[7][t1*(T*N*N) + t2*(N*N) + i1* N + i2] = task_solver.addLe( task_solver.sum( z[t1*(T*N*N) + t2*(N*N) + i1* N + i2], task_solver.negative(p[t1*N+ i1])), 0);
*/
            //add constraint 8: z_{t, t', i, i'}- p_{t', i'} \leq 0
            rng[8]= new IloRange[T*T*N*N];
            for(int t1=0; t1<T; t1++)
                for(int t2=0; t2<T; t2++)
                    for(int i1=0; i1< N; i1++)
                        for(int i2=0; i2< N; i2++)
                            rng[8][t1*(T*N*N) + t2*(N*N) + i1* N + i2] = task_solver.addLe( task_solver.sum( z[t1*(T*N*N) + t2*(N*N) + i1* N + i2], task_solver.negative(p[t2*N+ i2])), 0);

            //add constraint 9: p_{t, i}+ p_{t', i'} - z_{t, t', i, i'} \leq 1
            rng[9]= new IloRange[T*T*N*N];
            for(int t1=0; t1<T; t1++)
                for(int t2=0; t2<T; t2++)
                    for(int i1=0; i1< N; i1++)
                        for(int i2=0; i2< N; i2++)
                            rng[9][t1*(T*N*N) + t2*(N*N) + i1* N + i2] = task_solver.addLe( task_solver.sum(p[t1*N + i1], p[t2*N + i2], task_solver.negative(z[t1*(T*N*N) + t2*(N*N) + i1* N + i2])), 1 );

            //add constraint 10: \delta_{t, i} - \Omega_i p_{t, i} \leq 0
            rng[10]= new IloRange[T*N];
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    rng[10][t*N+ i] = task_solver.addLe( task_solver.sum( delta[t*N+i], task_solver.negative(task_solver.prod(Omega[i], p[t*N+i])) ), 0);

            //add constraint 11: \delta_{t, i} - \epsilon_i \leq 0
            rng[11]= new IloRange[T*N];
            for(int t=0; t< T; t++)
                for(int i=0; i<N; i++)
                    rng[11][t*N+i] = task_solver.addLe( task_solver.sum( delta[t*N+i], task_solver.negative(epsilon[i]) ), 0 );

            //add constraint 12: \epsilon_i - \delta_{t, i} + \Omega_i p_{t, i} \leq \Omega
            rng[12]= new IloRange[T*N];
            for(int t=0; t<T; t++)
                for(int i=0; i<N; i++)
                    rng[12][t*N+i] = task_solver.addLe( task_solver.sum( epsilon[i], task_solver.negative(delta[t*N+i]), task_solver.prod(Omega[i], p[t*N+i]) ), Omega[i] );

            boolean[][][] flag= new boolean[T][T][N];
            for(int t1=0; t1<T; t1++)
                for(int t2=0; t2<T; t2++)
                    for(int i=0; i<N; i++)
                        flag[t1][t2][i]=true;

            for (;;)
            {
                if(task_solver.solve())
                {
                    double[] z_current = task_solver.getValues(z); //Note here we define z_current as double
                    double[] delta_current= task_solver.getValues(delta);

                    int[] violated_index= violate(T, R, N, utility, z_current , delta_current , flag);
                    if(violated_index[0]==-1) //it means all no constraint violated in the oracle
                    {
                        out.write("Here is the x matrix: \r\n");
                        double[]x_final= task_solver.getValues(x);
                        double[]epsilon_final= task_solver.getValues(epsilon);
                        for(int i=0; i<N; i++)
                        {
                            for(int r=0; r<R; r++)
                            {
                                out.write(x_final[i*R + r]+ " ");
                            }
                            out.write("\r\n");
                        }
                        out.write("The epsilon is: \r\n");
                        for(int i=0; i<N; i++)
                            out.write(epsilon_final[i] + " ");
                        out.write("\r\n");



                        break;
                    }
                    else
                    {
                        int t1= violated_index[0];
                        int t2= violated_index[1];
                        int i= violated_index[2];

                        double[] temp_vector1= new double[T*T*N*N];
                        double[] temp_vector2= new double[T*T*N*N];
                        for(int j=0; j<T*T*N*N; j++)
                        {
                            temp_vector1[j]=0;
                            temp_vector2[j]=0;
                        }
                        for(int j=0; j<N; j++)
                        {
                            temp_vector1[t1 *(T*N*N) + t1*(N*N) + i*N + j]= utility[i][j];
                            temp_vector2[t1 *(T*N*N) + t2*(N*N) + i*N + j]= utility[i][j];
                        }

                        rng[3][t1*(T*N)+ t2*N+ i]= task_solver.addGe(task_solver.sum(task_solver.scalProd(temp_vector1, z), task_solver.negative(task_solver.scalProd(temp_vector2, z)), delta[t1*N+ i]), 0);
                        flag[t1][t2][i]=false;

                    }
                }

                else
                {
                    out.write("The problem is infeasible.\r\n");
                    break;
                }
            }

        }catch (IloException e) {
            System.err.println("Concert exception '" + e + "' caught");
        }

        out.flush();
        out.close();


    }

    static int[] violate (int T, int R, int N, double[][]utility, double[]z, double[]delta, boolean[][][]flag)
    {   //flag means whether a constraint is violated, if it is not violated, then we can set it false
        //here we check the constraint 3
        int[] min_index= new int[3]; //here we can denote the three index of t, t', and i;
        double min_violation=9999999999.0;
        int false_sum=0;

        for(int t1=0; t1<T; t1++)
            for(int t2=0; t2<T; t2++)
                for(int i=0; i<N; i++)
                {
                    if(flag[t1][t2][i] == false) {
                        false_sum++;
                        continue;
                    }
                    double current_sum=0;
                    for(int j=0; j<N; j++)
                        current_sum+= utility[i][j] * (z[t1*(T*N*N)+ t1*(N*N)+i*N+j] - z[t1*(T*N*N)+ t2*(N*N)+i*N+j]) + delta[t1*N + i];

                    if(current_sum<min_violation)
                    {
                        min_index[0]= t1;
                        min_index[1]= t2;
                        min_index[2]= i;

                        min_violation=current_sum;
                    }

                }

        if(min_violation>=0 || false_sum== T*T*N)
            min_index[0]=-1; //which denote there is no violation
        return min_index;
    }



}
