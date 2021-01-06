#include <fstream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{

    MPI_Init(NULL, NULL);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ifstream stream("../test/mpi_project_dev0.tsv");

    int P;
    stream >> P;

    int N;
    int A;
    int M;
    int T;
    stream >> N >> A >> M >> T;
    int num_of_instance = N / (P - 1);

    int attribute_per_processor = (A + 1) * (N / (P - 1));
    int instance_per_processor = N / (P - 1);

    double allAttributes[N * (A + 1)];
    double attributes[instance_per_processor][A+1];

    if (rank == 0)
    {
        for (int i = 0; i < N * (A + 1); i++)
        {
            stream >> allAttributes[i];
            // for (int k = 0; k < A + 1; k++)
            // {
            //     cout << attributes[k] << " ";
            // }
            // cout << endl;
        }
        // for(int i=0; i<N*A+N; i++) {
        //     cout << allAttributes[i] << " ";
        //     if((i+1)%(A+1) == 0) {
        //         cout << endl;
        //     }
        // }
    }

    cout << "---------------------------------------------" << endl;
    MPI_Scatter(allAttributes, attribute_per_processor, MPI_DOUBLE, attributes, attribute_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank != 0)
    {
        // for(int i=0; i<instance_per_processor; i++) {
        //     for(int j=0; j<A+1; j++) {
        //         cout << attributes[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // for(int i=0; i<attribute_per_processor; i++) {
        //     cout << attributes[i] << " ";
        //     if((i+1)%(A+1) == 0) {
        //         cout << endl;
        //     }
        // }
        double W[A];
        int nearest_hit = INT32_MAX;
        int nearest_miss = INT32_MAX;
        int nearest_hit_instance;
        int nearest_miss_instance;
        int nearest_value;
        for(int i=0; i<M; i++) {
            for(int j=0; j<instance_per_processor; j++) {
                if(i==j) 
                    continue;
                for(int k=0; k<A+1; k++) {
                    nearest_value += abs(attributes[i][k] - attributes[j][k]);
                }
                if(attributes[i][10] == attributes[j][10] ) 
                {
                    if(nearest_value < nearest_hit) 
                        nearest_hit = nearest_value;
                        nearest_hit_instance = j;
                } 
                else
                {
                    if(nearest_value < nearest_miss)
                        nearest_miss = nearest_value;
                        nearest_miss_instance = j;
                }
            }
             

        }
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
