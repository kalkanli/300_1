#include <fstream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{

    ifstream stream("../test/mpi_project_dev0.tsv");

    MPI_Init(NULL, NULL);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int P;
    stream >> P;

    int N;
    int A;
    int M;
    int T;
    int num_of_instance = N / (P - 1);
    stream >> N >> A >> M >> T;

    if (rank == 0)
    {
        double *instances[N];

        for (int i = 1; i <= N; i++)
        {
            double *attributes = new double[A + 1];
            for (int j = 1; j <= A; j++)
            {
                stream >> attributes[j];
            }
            stream >> attributes[0];
            instances[i - 1] = attributes;
            // for (int k = 0; k < A + 1; k++)
            // {
            //     cout << attributes[k] << " ";
            // }
            // cout << endl;
        }
        for (int a = 0; a < 50; a++)
        {
            double *ptr = instances[a];
            for (int b = 0; b < 11; b++)
            {
                cout << ptr[b] << " ";
            }
            cout << endl;
        }
        cout << sizeof(instances) / sizeof(instances[0]) << endl;
        for (int i = 0; i < N; i++)
        {
             
            // for(int r=0; r<A+1; r++) {
            //     cout << send[r] << " ";
            // }
            // cout << "################ " << i/num_of_instance+1 << endl;
            MPI_Send(
                instances[i],
                A + 1,
                MPI_DOUBLE,
                /*i/num_of_instance+*/1,
                0,
                MPI_COMM_WORLD);
        }
    }
    else
    {
        double *instances[N];
        for (int i = 0; i < N; i++)
        {
            double *attributes = new double[A + 1];
            MPI_Recv(
                attributes,
                A + 1,
                MPI_DOUBLE,
                0,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
            instances[0] = attributes;

            for (int j = 0; j < A + 1; j++)
            {
                cout << instances[i][j] << " ";
            }
            cout << "---------------" << rank << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
