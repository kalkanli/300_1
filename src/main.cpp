#include <fstream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <queue>
#include <mpi.h>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ifstream stream(argv[1]);

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

    /*
    *   Scatters the array into slave processors. 54-63
    */
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

    int maximumIndices[T];
    int maxAttributes[T*P];
    if (rank != 0)
    {

        double max[A];
        double min[A];
        double W[A];
        for (int u = 0; u < A; u++)
        {
            W[u] = 0.0;
            max[u] = INTMAX_MIN;
            min[u] = INTMAX_MAX;
        }
        for (int j = 0; j < instance_per_processor; j++)
        {
            for (int k = 0; k < A; k++)
            {
                if (attributes[j][k] > max[k])
                    max[k] = attributes[j][k];
                if (attributes[j][k] < min[k])
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

            for (int a = 0; a < A; a++)
            {
                W[a] += ((abs(attributes[i][a] - attributes[nearest_miss_instance][a])) / (max[a] - min[a])) / M;
                W[a] -= ((abs(attributes[i][a] - attributes[nearest_hit_instance][a])) / (max[a] - min[a])) / M;
            }
        }

        //------------------------------------------------------------------------------------------------------------
        double minimumWeight = INTMAX_MAX;
        int minimumWeightIndex = 0;

        for (int t = 0; t < T; t++)
        {
            maximumIndices[t] = t;
            if (W[t] < minimumWeight)
            {
                minimumWeightIndex = t;
                minimumWeight = W[t];
            }
        }
        for (int a = T; a < A; a++)
        {
            if (W[a] > minimumWeight)
            {
                maximumIndices[minimumWeightIndex] = a;
                minimumWeight = W[a];

                for (int t = 0; t < T; t++)
                {
                    if (W[maximumIndices[t]] < minimumWeight)
                    {
                        minimumWeight = W[maximumIndices[t]];
                        minimumWeightIndex = t;
                    }
                }
            }
        }
        sort(maximumIndices, maximumIndices + T);
        printf("Slave P%d: ", rank);
        for (int b = 0; b < T; b++)
        {
            printf("%d ", maximumIndices[b]);
        }
        cout << endl;
    }
    MPI_Gather(maximumIndices, T, MPI_INT, maxAttributes, T, MPI_INT, 0, MPI_COMM_WORLD);



    if (rank == 0)
    {
        printf("Master P0: ");
        sort(maxAttributes+T, maxAttributes+P*T);
        for(int i = T; i < (P)*T-1; i++) {
            if(maxAttributes[i] != maxAttributes[i+1])
                printf("%d ", maxAttributes[i]);
        }
        printf("%d\n", maxAttributes[P*T-1]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
