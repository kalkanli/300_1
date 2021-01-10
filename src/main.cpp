#include <fstream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <queue>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
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
    stream >> N >> A;
    int num_of_instance = N / (P - 1);

    int attribute_per_processor = (A + 1) * (N / (P - 1));
    int instance_per_processor = N / (P - 1);
    int receiveCount = attribute_per_processor;

    double allAttributes[N * (A + 1)];
    double attributes[instance_per_processor][A + 1];

    if (rank == 0)
    {
        stream >> M >> T;
        for (int i = 0; i < N * (A + 1); i++)
        {
            stream >> allAttributes[i];
        }
    }

    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int sendcounts[P];
    int displs[P];
    for (int x = 0; x < P; x++)
    {
        sendcounts[x] = attribute_per_processor;
        displs[x] = (x - 1) * attribute_per_processor;
    }
    sendcounts[0] = 0;
    displs[0] = 0;
    MPI_Scatterv(allAttributes, sendcounts, displs, MPI_DOUBLE, attributes, attribute_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 1)
    {

        double max[A];
        double min[A];
        double W[A];
        for (int u = 0; u < A; u++)
        {
            W[u] = 0.0;
            max[u] = -5000;
            min[u] = 5000;
        }
        for (int j = 0; j < instance_per_processor; j++)
            {
                
                for (int k = 0; k < A; k++) {
                    if(attributes[j][k] > max[k])
                        max[k] = attributes[j][k];
                    if(attributes[j][k] < min[k])
                        min[k] = attributes[j][k];
                       
                }
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
                    // if (max[k] < attributes[j][k])
                    // {
                    //     max[k] = attributes[j][k];
                    // }

                    // if (min[k] > attributes[j][k])
                    // {
                    //     min[k] = attributes[j][k];
                    // }
                    nearest_value += abs(attributes[i][k] - attributes[j][k]);
                }

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
    
            // cout << "-----------------------------------------------" << endl;
            // printf("%dP nearest miss for iteration %d : %d\n", rank, i, nearest_miss_instance);
            // printf("%dP nearest hit for iteration %d : %d\n", rank, i, nearest_hit_instance);
            // cout << "-----------------------------------------------" << endl;
            for (int a = 0; a < A; a++)
            {
                // double diff = 0;
                // diff += abs(attributes[i][a] - attributes[nearest_miss_instance][a]);
                // diff -= abs(attributes[i][a] - attributes[nearest_hit_instance][a]);
                // diff /= (max[a] - min[a]);
                // diff /= M;
                // W[a] += diff;
                //printf("((%f - %f )/%f)/%d\n", attributes[i][a],attributes[nearest_miss_instance][a],max[a] - min[a],M);
                W[a] += ((abs(attributes[i][a] - attributes[nearest_miss_instance][a])) / (max[a] - min[a]))/ M;
                W[a] -= ((abs(attributes[i][a] - attributes[nearest_hit_instance][a])) / (max[a] - min[a]))/ M;
                cout << W[a] << " ";
            }
            cout << endl;
        }

        // for (int a = 0; a < A; a++)
        // {
        //     cout << min[a] << " ";
        // }
        // cout << " P" << rank << endl;

        // int result[T];

        // priority_queue<pair<double, int>> queue;
        // for (int i = 0; i < A; i++)
        // {
        //     queue.push(make_pair(W[i], i));
        // }

        // printf("Slave P%d: ", rank);
        // for (int i = 1; i <= T; i++)
        // {
        //     result[i - 1] = queue.top().second;
        //     queue.pop();
        // }

        // bool swapped;
        // for (int i = 0; i < T - 1; i++)
        // {
        //     swapped = false;
        //     for (int j = 0; j < T - i - 1; j++)
        //     {
        //         if (result[j] > result[j + 1])
        //         {
        //             int temp = result[j];
        //             result[j] = result[j + 1];
        //             result[j + 1] = temp;
        //             swapped = true;
        //         }
        //     }
        //     if (swapped == false)
        //             break;
        // }

        // for (int i = 1; i <= T; i++)
        // {
        //     cout << result[i-1] << " ";
        // }
        // cout << endl;
    }
   // MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
