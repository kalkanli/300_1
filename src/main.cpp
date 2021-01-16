/*
 *
 * Student Name: A. Taha Kalkanlı
 * Student Number: 2017400189
 * Compile Status: Compiling
 * Program Status: Working
 * Notes: in order to run this program after compiling you should use the terminal code:
 *        "mpirun -np <processor_number> --oversubscribe main <input_file_path>"
 * 
 */
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

    int N, A, M, T;
    stream >> N >> A;

    int attributePerProcessor = (A + 1) * (N / (P - 1));
    int instancePerProcessor = N / (P - 1);

    double allAttributes[N * (A + 1)];
    double attributes[instancePerProcessor][A + 1];

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
        sendcounts[x] = attributePerProcessor;
        displs[x] = (x - 1) * attributePerProcessor;
    }
    sendcounts[0] = 0;
    displs[0] = 0;
    MPI_Scatterv(allAttributes, sendcounts, displs, MPI_DOUBLE, attributes, attributePerProcessor, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int highestWeightedAttributes[T];
    int gatheredHighestWeightedAttributes[T*P];
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
        for (int j = 0; j < instancePerProcessor; j++)
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
            double nearestHit = INTMAX_MAX;
            double nearestMiss = INTMAX_MAX;
            int nearestHitInstance;
            int nearestMissInstance;

            for (int j = 0; j < instancePerProcessor; j++)
            {
                if (i == j) continue;
                int nearestValue = 0;
                for (int k = 0; k < A; k++)
                {
                    nearestValue += abs(attributes[i][k] - attributes[j][k]);
                }

                if (attributes[i][A] == attributes[j][A])
                {
                    if (nearestValue < nearestHit)
                    {
                        nearestHit = nearestValue;
                        nearestHitInstance = j;
                    }
                }
                else
                {
                    if (nearestValue < nearestMiss)
                    {
                        nearestMiss = nearestValue;
                        nearestMissInstance = j;
                    }
                }
            }

            for (int a = 0; a < A; a++)
            {
                W[a] += ((abs(attributes[i][a] - attributes[nearestMissInstance][a])) / (max[a] - min[a])) / M;
                W[a] -= ((abs(attributes[i][a] - attributes[nearestHitInstance][a])) / (max[a] - min[a])) / M;
            }
        }

        double minimumWeight = INTMAX_MAX;
        int minimumWeightIndex = 0;

        for (int t = 0; t < T; t++)
        {
            highestWeightedAttributes[t] = t;
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
                highestWeightedAttributes[minimumWeightIndex] = a;
                minimumWeight = W[a];

                for (int t = 0; t < T; t++)
                {
                    if (W[highestWeightedAttributes[t]] < minimumWeight)
                    {
                        minimumWeight = W[highestWeightedAttributes[t]];
                        minimumWeightIndex = t;
                    }
                }
            }
        }
        sort(highestWeightedAttributes, highestWeightedAttributes + T);
        printf("Slave P%d: ", rank);
        for (int b = 0; b < T; b++)
        {
            printf("%d ", highestWeightedAttributes[b]);
        }
        cout << endl;
    }
    MPI_Gather(highestWeightedAttributes, T, MPI_INT, gatheredHighestWeightedAttributes, T, MPI_INT, 0, MPI_COMM_WORLD);



    if (rank == 0)
    {
        printf("Master P0: ");
        sort(gatheredHighestWeightedAttributes+T, gatheredHighestWeightedAttributes+P*T);
        for(int i = T; i < (P)*T-1; i++) {
            if(gatheredHighestWeightedAttributes[i] != gatheredHighestWeightedAttributes[i+1])
                printf("%d ", gatheredHighestWeightedAttributes[i]);
        }
        printf("%d\n", gatheredHighestWeightedAttributes[P*T-1]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
