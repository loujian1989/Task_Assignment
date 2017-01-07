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


public class main_minimize_cost {

    public static void main(String[] args)throws Exception
    {
        Scanner cin=new Scanner(new File("test_in"));
        File writename = new File("test_out");
        writename.createNewFile();
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));


        int T= cin.nextInt(); //number of task
        int R= cin.nextInt(); //number of resource types
        int N= cin.nextInt(); //number of players
        int K= cin.nextInt(); //the maximum size of a team

        double[][] resource_need= new double[T][R];  //R_{t,r}
        double[][] lost_for_resource= new double[N][R];  //l_{i,r}
        double[][] upper_bound_player= new double[N][R];  //U_{i,r}

        //input the data from the file
        for(int t=0; t<T; t++)
            for(int r=0; r<R; r++)
                resource_need[t][r]= cin.nextDouble();

        for(int i=0; i<N; i++)
            for(int r=0; r<R; r++)
                lost_for_resource[i][r]= cin.nextDouble();

        for(int i=0; i<N; i++)
            for(int r=0; r<R; r++)
                upper_bound_player[i][r] = cin.nextDouble();

        minimize_cost MC= new minimize_cost(T, R, N, K, resource_need, lost_for_resource, upper_bound_player);

        MC.solve_problem();

        out.flush();
        out.close();

    }

}
