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

    ifstream stream("../test/test.txt");

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
    double attributes[instance_per_processor][A + 1];

    if (rank == 0)
    {
        for (int i = 0; i < N * (A + 1); i++)
        {
            stream >> allAttributes[i];
        }
    }

    int send[P] = {attribute_per_processor};
    send[0] = 0;
    int displs[P];
    for(int i=1; i<P; i++) {
        displs[i] = (i-1)*attribute_per_processor;
    }
    displs[0] = 0;

    MPI_Scatter(allAttributes, attribute_per_processor, MPI_DOUBLE, attributes, attribute_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        cout << "hayda" << endl;
        // for(int i=0; i<instance_per_processor; i++) {
        //     for(int j=0; j<A+1; j++) {
        //         cout << attributes[i][j] << " ";
        //     }
        //     cout << endl;
        // }

        double max[A];
        double min[A];
        double W[A];
        for (int u = 0; u < A; u++)
        {
            W[u] = 0;
            max[u] = INTMAX_MIN;
            min[u] = INTMAX_MAX;
        }

        for (int i = 0; i < M; i++)
        {
            double nearest_hit = INTMAX_MAX;
            double nearest_miss = INTMAX_MAX;
            int nearest_hit_instance;
            int nearest_miss_instance;

            for (int j = 0; j < instance_per_processor; j++)
            {
                if (i == j)
                    continue;
                int nearest_value = 0;
                for (int k = 0; k < A; k++)
                {
                    if (max[k] < attributes[j][k])
                    {
                        max[k] = attributes[j][k];
                    }

                    if (min[k] > attributes[j][k])
                    {
                        min[k] = attributes[j][k];
                    }
                    nearest_value += abs(attributes[i][k] - attributes[j][k]);
                    // cout << nearest_value << " ";
                }
                // cout << endl;
                // cout << attributes[i][A] << "==" << attributes[j][A] << endl;

                if (attributes[i][A] == attributes[j][A])
                {
                    if (nearest_value < nearest_hit)
                    {
                        nearest_hit = nearest_value;
                        nearest_hit_instance = j;
                    }
                }
                else
                {
                    if (nearest_value < nearest_miss)
                    {
                        nearest_miss = nearest_value;
                        nearest_miss_instance = j;
                    }
                }
            }
            cout << "-----------------------------------------------" << endl;
            cout << "Nearest Miss Instance: " << nearest_miss_instance << endl;
            cout << "Nearest Hit Instance:  " << nearest_hit_instance << endl;
            cout << "-----------------------------------------------" << endl;

            for (int a = 0; a < A; a++)
            {
                W[a] += abs(attributes[i][a] - attributes[nearest_miss_instance][a]);
                W[a] -= abs(attributes[i][a] - attributes[nearest_hit_instance][a]);
                W[a] /= max[a]-min[a];
                W[a] /= M;
            }
        }



        for (int i = 0; i < A; i++)
        {
            cout << W[i] << " ";
        }
        cout << endl;
        
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
