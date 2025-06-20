#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdbool.h>
#include <assert.h> 
#include <time.h> 

int* minPath;
int minCost = INT_MAX;
int pathCount = 0; 

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

int** allocateDistanceMatrix(int n) {
    int** matrix = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        matrix[i] = (int*)malloc(n * sizeof(int));
    }
    return matrix;
}

void freeDistanceMatrix(int** matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

int calculateCost(int* route, int** distanceMatrix, int n) {
    int cost = 0;
    for (int i = 0; i < n - 1; i++) {
        assert(route[i] >= 0 && route[i] < n);
        assert(route[i + 1] >= 0 && route[i + 1] < n);
        cost += distanceMatrix[route[i]][route[i + 1]];
    }
    assert(route[n - 1] >= 0 && route[n - 1] < n);
    assert(route[0] >= 0 && route[0] < n);
    cost += distanceMatrix[route[n - 1]][route[0]]; 
    return cost;
}

void printPath(int* route, int n, int cost) {
    for (int i = 0; i < n; i++) {
        printf("%d -> ", route[i]);
    }
    printf("%d | Cost: %d\n", route[0], cost);
}


void tryAllRoutes(int* route, int** distanceMatrix, int start, int n) {
    if (start == n - 1) {
        pathCount++; 
        int cost = calculateCost(route, distanceMatrix, n);
        printPath(route, n, cost); 
        if (cost < minCost) {
            minCost = cost;
            memcpy(minPath, route, sizeof(int) * n);
            minPath[n] = route[0];
        }
        return;
    }

    for (int i = start; i < n; i++) {
        int temp = route[start];
        route[start] = route[i];
        route[i] = temp;

        tryAllRoutes(route, distanceMatrix, start + 1, n);

        temp = route[start];
        route[start] = route[i];
        route[i] = temp;
    }
}

double calculateDistance(int x1, int y1, int x2, int y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

void generateDistanceMatrix(int coordinates[][2], int n, int** distanceMatrix) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                distanceMatrix[i][j] = 0;
            } else {
                distanceMatrix[i][j] = (int)calculateDistance(coordinates[i][0], coordinates[i][1], coordinates[j][0], coordinates[j][1]);
            }
        }
    }
}

void saveResultToFile(const char* fileName, int coordinates[][2], int* bestPath, int n) {
    FILE* file = fopen(fileName, "w");
    if (file == NULL) {
        printf("Error: Could not create result file\n");
        return;
    }

    for (int i = 0; i <= n; i++) {
        int id = bestPath[i];
        fprintf(file, "%d %d\n", coordinates[id][0], coordinates[id][1]);
    }

    fclose(file);
}

int main(int argc, char *argv[]) {
    int CITIES = argc > 1 ? atoi(argv[1]) : 10;

    int coordinates[CITIES][2];

    for (int i = 0; i < CITIES; i++)
    {
        coordinates[i][0] = __coordinates[i][0];
        coordinates[i][1] = __coordinates[i][1];
    }


    int n = sizeof(coordinates) / sizeof(coordinates[0]);

    int** distanceMatrix = allocateDistanceMatrix(n);

    minPath = (int*)malloc((n + 1) * sizeof(int));

    generateDistanceMatrix(coordinates, n, distanceMatrix);

    int cities[n];
    for (int i = 0; i < n; i++) cities[i] = i;

    tryAllRoutes(cities, distanceMatrix, 2, n);
    
    printf("Minimum Cost: %d\n", minCost);
    printf("Best Route: ");
    for (int i = 0; i <= n; i++) {
        printf("%d ", minPath[i]);
        if (i < n) printf("-> ");
    }
    printf("\n");

    printf("Number of paths tried: %d\n", pathCount);

    saveResultToFile("result.txt", coordinates, minPath, n);

    freeDistanceMatrix(distanceMatrix, n);
    free(minPath);

    return 0;
}
