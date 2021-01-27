/*
 *
 * Student Name:    A. Taha Kalkanlı
 * Student Number:  2017400189
 * Compile Status:  Compiling
 * Program Status:  Working
 * Notes:           in order to run this program after compiling you should use the terminal code:
 *                      "mpirun -np <processor_number> --oversubscribe main <input_file_path>"
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

    int attributePerProcessor = (A + 1) * (N / (P - 1)); // number of attributes + class value that will be processed by each processor.
    int instancePerProcessor = N / (P - 1); // number of instances value that will be processed by each processor.

    double allAttributes[N * (A + 1)]; // all atributes + class values. this is used as a send buffer.
    double attributes[instancePerProcessor][A + 1]; // used as a recieve buffer.

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
    sendcounts[0] = 0; // number of elements that will be sent to each processor from allAttributes[]
    displs[0] = 0; // starting index of elements that will be sent to each processors
    MPI_Scatterv(allAttributes, sendcounts, displs, MPI_DOUBLE, attributes, attributePerProcessor, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

    int highestWeightedAttributes[T]; // send buffer for gatterv method. data that will be sent to master processor.
    int gatheredHighestWeightedAttributes[T*P]; // recieve buffer for master processor.
    
    if (rank != 0) // ıf rank != 0 means not the master processor.
    {
        double max[A]; // holds the maximum features along scattered instances in every process.
        double min[A]; // holds the minimum features along scattered instances in every process.
        double W[A]; // weight array.
        
        for (int u = 0; u < A; u++) // initializes the arrays above before using them.
        {
            W[u] = 0.0;
            max[u] = INTMAX_MIN;
            min[u] = INTMAX_MAX;
        }

        for (int j = 0; j < instancePerProcessor; j++) // finds the minimum and maximum of each feature.
        {
            for (int k = 0; k < A; k++)
            {
                if (attributes[j][k] > max[k])
                    max[k] = attributes[j][k];
                if (attributes[j][k] < min[k])
                    min[k] = attributes[j][k];
            }
        }

        for (int i = 0; i < M; i++) // target instance is M.
        {
            double nearestHit = INTMAX_MAX; // nearest manhattan distance to hit instance.
            double nearestMiss = INTMAX_MAX; // nearest manhattan distance to miss instance.
            int nearestHitInstance; // instance # of nearest hit.
            int nearestMissInstance; // instance # of nearest miss.

            for (int j = 0; j < instancePerProcessor; j++) // finds the nearest hit and nearest miss instances. then updates the W[] 
            {
                if (i == j) continue; // if the same instance just skip it.
                int manhattanDistance = 0; // Manhattan Distance
                for (int k = 0; k < A; k++)
                {
                    manhattanDistance += abs(attributes[i][k] - attributes[j][k]);
                }

                if (attributes[i][A] == attributes[j][A]) // if it is a hit.
                {
                    if (manhattanDistance < nearestHit)
                    {
                        nearestHit = manhattanDistance;
                        nearestHitInstance = j;
                    }
                }
                else // if it is a miss.
                {
                    if (manhattanDistance < nearestMiss)
                    {
                        nearestMiss = manhattanDistance;
                        nearestMissInstance = j;
                    }
                }
            }

            for (int a = 0; a < A; a++) // updates the weight array with the found nearest miss and nearest hit instances.
            {
                W[a] += ((abs(attributes[i][a] - attributes[nearestMissInstance][a])) / (max[a] - min[a])) / M;
                W[a] -= ((abs(attributes[i][a] - attributes[nearestHitInstance][a])) / (max[a] - min[a])) / M;
            }
        }

        double minimumWeight = INTMAX_MAX; // used as a lower bound for minimum weights.
        int minimumWeightIndex = 0; // index of the minimum weighted instance in highestWeightedAttributes[]

        for (int t = 0; t < T; t++) // puts the first T elements of the W[] into highestWeightedAttributes[]
        {
            highestWeightedAttributes[t] = t;
            if (W[t] < minimumWeight) // holds the minimum weighted instance and the weight of it in the case of replacement.
            {
                minimumWeightIndex = t;
                minimumWeight = W[t];
            }
        }

        // after the Tth element it checks the other elements of the W array and 
        // replaces the minimum element in highestWeightedAttributes with the  
        // element of W[] larger than the minimumWeight. 
        // After that it finds the new minimumWeight and updates the minimumWeightIndex.
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

        sort(highestWeightedAttributes, highestWeightedAttributes + T); // sorts the highestWeightedAttributes.
        printf("Slave P%d: ", rank); // prints the slave part of the output.
        for (int b = 0; b < T; b++)
        {
            printf("%d ", highestWeightedAttributes[b]);
        }
        cout << endl;
    }
    MPI_Gather(highestWeightedAttributes, T, MPI_INT, gatheredHighestWeightedAttributes, T, MPI_INT, 0, MPI_COMM_WORLD);



    if (rank == 0)
    {
        printf("Master P0: "); // prints the master part of the output.
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
