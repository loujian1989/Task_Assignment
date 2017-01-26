/**
 * Created by loujian on 1/24/17.
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Random;


public class Generate_Cases {

    public static void main(String[] args)throws Exception
    {
        File writename = new File("N20T4R3.txt");
        writename.createNewFile();
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));

        int max_resource_needed=4; //the resource_needed range in the experiment
        int max_lost= 10; //the lost range assumed in the experiment
        int max_utility=10; //the utility range assumed in the experiment
        int weight_player=1; //the weight for each player

        int N=20;
        int T=4;
        int R=3;

        out.write(T+ "\r\n");
        out.write(R+ "\r\n");
        out.write(N+ "\r\n");
        out.write("\r\n");

        int[][] resource_need= new int[T][R];  //R_{t,r}
        int[][] lost_for_resource= new int[N][R];  //l_{i,r}
        int[][] utility = new int[N][N]; //u_i(j);
        int[] Omega = new int[N]; //Omega_i, denote the upper bound of epsilon
        int[] alpha = new int[N]; //alpha_i, denote the weight of epsilon for each player
        int[][]upper_bound_player = new int[N][R];


        Random rd = new Random();

        for(int t=0; t<T; t++)
        {
            for(int r=0; r<R; r++) //Here we generate the resource needed for some task
            {
                resource_need[t][r]= rd.nextInt(max_resource_needed);
                out.write(resource_need[t][r] + " ");
            }
            out.write("\r\n");
        }

        out.write("\r\n");

        for(int i=0; i<N; i++)
        {
            for(int r=0; r<R; r++)
            {
                lost_for_resource[i][r]= rd.nextInt(max_lost);
                out.write(lost_for_resource[i][r]+ " ");
            }
            out.write("\r\n");
        }

        out.write("\r\n");

        for(int i=0; i<N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if(i==j)
                    utility[i][j]=0;
                else
                    utility[i][j]=rd.nextInt(max_utility);
                out.write(utility[i][j]+ " ");
            }
            out.write("\r\n");
        }

        out.write("\r\n");

        for(int i=0; i<N; i++)
        {
            alpha[i]=weight_player;
            out.write(alpha[i]+ " ");
        }
        out.write("\r\n");

        for(int i=0; i<N; i++)
        {
            for(int r=0; r<R; r++)
            {
                upper_bound_player[i][r] = max_resource_needed;
                out.write(upper_bound_player[i][r]+ " ");
            }
            out.write("\r\n");
        }

        out.flush();
        out.close();
    }



}
