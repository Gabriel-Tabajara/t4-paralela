#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include "mpi.h"
#include <omp.h>

int RESULT_TAG = 0;
int REQUEST_TAG = 1;
int KILL_TAG = 2;
int WORK_TAG = 3;

int __coordinates[][2] = {
    {10, 5},
    {29, 24},
    {21, 33},
    {33, 95},
    {49, 18},
    {50, 2},
    {62, 22},
    {73, 86},
    {29, 68},
    {99, 33},
    {9,20},
    {12,98},
    {15, 68},
    {18, 45},
    {20, 10},
    {25, 30},
    {30, 50},
    {35, 70},
    {40, 90},
    {45, 15},
    {50, 25},
    {55, 35},
    {60, 45},
    {65, 55},
    {70, 65},
    {75, 75},
    {80, 85},
    {85, 95},
    {90, 5},
    {95, 15}
};

int *minPath;
int minCost = INT_MAX;
int pathCount = 0;

int **allocateDistanceMatrix(int n)
{
    int **matrix = (int **)malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++)
    {
        matrix[i] = (int *)malloc(n * sizeof(int));
    }
    return matrix;
}

int *mergeCompleteRoute(int *route, int n, int *baseRoute)
{
    int *completeRoute = (int *)malloc((n) * sizeof(int)); // n já inclui as 3 fixas
    completeRoute[0] = baseRoute[0];
    completeRoute[1] = baseRoute[1];
    completeRoute[2] = baseRoute[2];
    for (int i = 3; i < n; i++)
    {
        completeRoute[i] = route[i - 3];
    }
    return completeRoute;
}


int calculateCost(int *route, int **distanceMatrix, int n, int *baseRoute)
{
    int *completeRoute = mergeCompleteRoute(route, n, baseRoute);
    int cost = 0;

    for (int i = 0; i < n - 1; i++)
    {
        assert(completeRoute[i] >= 0 && completeRoute[i] < n);
        assert(completeRoute[i + 1] >= 0 && completeRoute[i + 1] < n);
        cost += distanceMatrix[completeRoute[i]][completeRoute[i + 1]];
    }

    assert(completeRoute[n - 1] >= 0 && completeRoute[n - 1] < n);
    cost += distanceMatrix[completeRoute[n - 1]][completeRoute[0]];

    free(completeRoute);
    return cost;
}


void printPath(int *route, int n, int cost, int *baseRoute)
{
    int *completeRoute = mergeCompleteRoute(route, n, baseRoute);

    for (int i = 0; i < n; i++)
    {
        printf("%d -> ", completeRoute[i]);
    }
    printf("%d | Cost: %d\n", completeRoute[0], cost);

    free(completeRoute);
}

void tryAllRoutes(int *route, int **distanceMatrix, int start, int n, int *baseRoute, int *iterationsLeft, int threadId, int numThreads)
{
    if (*iterationsLeft <= 0) {
        return;
    }

    if (start == n-1)
    {
        (*iterationsLeft)--;
        if (*iterationsLeft < 0) return;

        pathCount++;
        int cost = calculateCost(route, distanceMatrix, n + 3, baseRoute);
        // printPath(route, n + 3, cost, baseRoute);
        if (cost < minCost)
        {
            int *completeRoute = mergeCompleteRoute(route, n+3, baseRoute);
            #pragma omp critical
            {
                memcpy(minPath, completeRoute, sizeof(int) * (n+3));
                minCost = cost;
            }
            free(completeRoute);
        }
        return;
    }

    for (int i = start + threadId; i < n; i += numThreads)
    {
        if (*iterationsLeft <= 0) {
            return;
        }

        int temp = route[start];
        route[start] = route[i];
        route[i] = temp;

        tryAllRoutes(route, distanceMatrix, start + 1, n, baseRoute, iterationsLeft, threadId, numThreads);

        temp = route[start];
        route[start] = route[i];
        route[i] = temp;
    }
}

double calculateDistance(int x1, int y1, int x2, int y2)
{
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

void generateDistanceMatrix(int coordinates[][2], int n, int **distanceMatrix)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                distanceMatrix[i][j] = 0;
            }
            else
            {
                distanceMatrix[i][j] = (int)calculateDistance(coordinates[i][0], coordinates[i][1], coordinates[j][0], coordinates[j][1]);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int CITIES = argc > 1 ? atoi(argv[1]) : 10;

    int ITERATIONS = argc > 2 ? atoi(argv[2]) : 1000;

    int NUM_THREADS = argc > 3 ? atoi(argv[3]) : 8;

    int coordinates[CITIES][2];

    for (int i = 0; i < CITIES; i++)
    {
        coordinates[i][0] = __coordinates[i][0];
        coordinates[i][1] = __coordinates[i][1];
    }


    int n = sizeof(coordinates) / sizeof(coordinates[0]);

    int my_rank;
    int proc_n;
    int message;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);

    double start_time = MPI_Wtime();

    int **distanceMatrix = allocateDistanceMatrix(n);

    minPath = (int *)malloc((n + 3) * sizeof(int));

    generateDistanceMatrix(coordinates, n, distanceMatrix);

    //printf("Number of processes: %d\n", proc_n);
    if (my_rank == 0)
    {

        printf("Master has started\n");

        int qnt_restante = n - 1;

        int totalCombinations = (n - 1) * (n - 2); // combinações de i != j com i,j ≠ 0
        int currentPos = 0;

        int interationsLeft = ITERATIONS;
        int iterationsPerCombination = interationsLeft / totalCombinations;
        if (iterationsPerCombination == 0)
        {
            iterationsPerCombination = 1;
        } else if (interationsLeft % totalCombinations != 0)
        {
            iterationsPerCombination += 1;
        }

        // printf("Total combinations: %d, iterations left %d, iterations per combination: %d\n", totalCombinations, interationsLeft, iterationsPerCombination);
        while (interationsLeft > 0)
        {
            int *message = (int *)malloc((n + 1) * sizeof(int));
            message[0] = iterationsPerCombination < interationsLeft ? iterationsPerCombination : interationsLeft; 
            message[1] = 0;

            int i = (currentPos / (n - 2)) + 1;
            int j = (currentPos % (n - 2)) + 1;
            if (j >= i) j++;  // pular quando j == i

            if (j >= n) {
                currentPos++;
                free(message);
                continue;
            }

            message[2] = i;
            message[3] = j;

            int idx = 4;
            for (int k = 1; k < n; k++) {
                if (k != i && k != j)
                    message[idx++] = k;
            }

            // Espera worker pedir trabalho
            MPI_Status status;
            MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
            // printf("Master received request from %d \n", status.MPI_SOURCE);
            int sender_rank = status.MPI_SOURCE;
            MPI_Send(message, n+1, MPI_INT, sender_rank, WORK_TAG, MPI_COMM_WORLD);
            // printf("Master sent message to %d: ", sender_rank);
            // for (int k = 0; k < n + 1; k++)
            // {
            //     printf("%d ", message[k]);
            // }
            // printf("\n");
            free(message);
            currentPos++;
            interationsLeft -= iterationsPerCombination;
        }


        for (int i = 1; i < proc_n; i++)
        {
            // printf("Master sending kill signal to process %d\n", i);
            MPI_Send(0, 0, MPI_INT, i, KILL_TAG, MPI_COMM_WORLD);
        }

        for (int i = 1; i < proc_n; i++)
        {
            int *resultMessage = (int *)malloc((n + 1) * sizeof(int));

            MPI_Status status;
            MPI_Recv(resultMessage, n+1, MPI_INT, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
            // printf("Master received result from %d \n", status.MPI_SOURCE);

            int cost = resultMessage[n];
            if (cost < minCost)
            {
                minCost = cost;
                memcpy(minPath, resultMessage, sizeof(int) * n);
            }
            free(resultMessage);
        }

        printf("Master Best Route with min cost %d: ", minCost);
        for (int i = 0; i < n; i++)
        {
            printf("%d -> ", minPath[i]);
        }
        printf("%d ", minPath[0]);
        printf("\n");

        printf("Master finished\n");
    }
    else
    {
        printf("Process %d started\n", my_rank);
        while (true)
        {
            // printf("Process %d asking for work\n", my_rank);
            MPI_Send(0, 0, MPI_INT, 0, REQUEST_TAG, MPI_COMM_WORLD);   
    
            int *message = (int *)malloc((n+1) * sizeof(int));
            
            MPI_Status status;
            MPI_Recv(message, n+1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // printf("Process %d received message from master with tag %d\n", my_rank, status.MPI_TAG);

            if (status.MPI_TAG == KILL_TAG)
            {
                int *finalMessage = (int *)malloc((n + 1)* sizeof(int));

                memcpy(finalMessage, minPath, sizeof(int) * (n + 1));
                finalMessage[n] = minCost;

                MPI_Send(finalMessage, n+1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);

                free(message);
                free(finalMessage);
                printf("Process %d finished\n", my_rank);
                break;
            }
    
            // printf("Process %d received message from %d message: ", my_rank, status.MPI_SOURCE);
    
            // for (int i = 0; i < n + 1; i++)
            // {
            //     printf("%d ", message[i]);
            // }
    
            // printf("\n");
    
            int cities[n];
            for (int i = 3; i < n; i++)
                cities[i - 3] = message[i+1];

            int baseRoute[3] = {message[1], message[2], message[3]};
            int iterationsLeft = message[0];
            // printf("Process %d will send cities: ", my_rank);
            // for (int i = 0; i < n - 3; i++)
            // {
            //     printf("%d ", cities[i]);
            // }
            // printf("\n");
            // printf("Process %d will calculate routes with base route: %d -> %d -> %d\n", my_rank, baseRoute[0], baseRoute[1], baseRoute[2]);
            // printf("Process %d has %d iterations left\n", my_rank, iterationsLeft);

            omp_set_num_threads(NUM_THREADS);

            int iterationsPerThread = iterationsLeft / NUM_THREADS;
            int rest = iterationsLeft % NUM_THREADS;
            int iterationsForThread[NUM_THREADS];
            for (int i = 0; i < NUM_THREADS; i++)
            {
                iterationsForThread[i] = iterationsPerThread + (i < rest ? 1 : 0);
            }

            #pragma omp parallel 
            {
                int threadId = omp_get_thread_num();
                int iterations = iterationsForThread[threadId];
                if (iterations > 0)
                {
                    int localCities[n];
                    for (int k = 0; k < n; k++) localCities[k] = cities[k];
                    tryAllRoutes(localCities, distanceMatrix, 0, n - 3, baseRoute, &iterations, threadId, NUM_THREADS);
                }
            }

            free(message);
        }
        
        //printf("Process %d Best Route with min cost %d: ", my_rank, minCost);
        // for (int i = 0; i <= n; i++)
        // {
            //printf("%d ", minPath[i]);
            // if (i < n)
                //printf("-> ");
        // }
        //printf("\n");
    }

    MPI_Finalize();

    double end_time = MPI_Wtime();

    if (my_rank == 0) {
        printf("Tempo de execução: %f s\n", end_time - start_time);
    }

    return 0;
}
