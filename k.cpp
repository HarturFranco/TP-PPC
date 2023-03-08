#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <cstring>

#include "mpi.h"

using namespace std;

vector<vector<float>> readCsvFile(const string& filename) {
    vector<vector<float>> data;

    ifstream infile(filename);
    string line;

    // skip cabeçalho
    getline(infile, line);

    // loop pelas linhas do arquivo
    while (getline(infile, line)) {
        istringstream iss(line);
        string token;

        

        // vetor de float com valores da linha
        vector<float> row;
        while (getline(iss, token, ',')) {
            row.push_back(stof(token));
        }

        // add linha em data
        data.push_back(row);
    }

    return data;
}


struct dataPoint
{
    int cluster;
    vector<float> features;

};

// distancia euclidiana
float distance(float* p1, float* p2, int n_features) {
    float sum = 0;
    for (int i = 0; i < n_features; i++) {
        sum += pow(p1[i] - p2[i], 2);
    }
    return sqrt(sum);
}

void assignClusters(vector<dataPoint>& dataPoints, vector<dataPoint>& centroids) {
    
    for (dataPoint& point : dataPoints) {
        float minDist = numeric_limits<float>::max();
        int cluster = 0;
        for (int i = 0; i < centroids.size(); i++) {
            float dist = distance(point.features.data(), centroids[i].features.data(), centroids[i].features.size());
            if (dist < minDist) {
                minDist = dist;
                cluster = i;
            }
        }
        point.cluster = cluster;
    }
}

// calcula centroid dos clusters
vector<dataPoint> calculateCentroids(vector<dataPoint>& dataPoints, int k) {
    vector<dataPoint> centroids(k);
    vector<int> clusterSizes(k);
    
    for(int i = 0; i < centroids.size(); i++){
        centroids[i].features = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    }
    
    for (dataPoint& point : dataPoints) {
        int cluster = point.cluster;
        clusterSizes[cluster]++;
        for (int i = 0; i < point.features.size(); i++) {
            centroids[cluster].features[i] += point.features[i];
        }
    }
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < centroids[i].features.size(); j++) {
            centroids[i].features[j] /= clusterSizes[i];
        }
    }
    return centroids;
}

vector<dataPoint> kMeans(vector<dataPoint>& dataPoints, int k, long maxIterations){

    vector<dataPoint> centroids(k);

    srand(time(0));
    for (int i = 0; i < k; i++) {
        centroids[i] = dataPoints[rand() % dataPoints.size()];
    }
//     centroids[0] = dataPoints[0];
//     centroids[1] = dataPoints[8];
//     centroids[2] = dataPoints[64];

    int iterations = 0;
    // while (iterations < maxIterations) {
    while (true) {
        assignClusters(dataPoints, centroids);
        vector<dataPoint> newCentroids = calculateCentroids(dataPoints, k);
        // se centroid não muda - termina
        bool convergence = true;
        for (int i = 0; i < k; i++) {
            if (distance(newCentroids[i].features.data(), centroids[i].features.data(), centroids[i].features.size()) > 0.0001) {
                convergence = false;
                break;
            }
        }
        if (convergence) {
            cout << "convergiu: " << iterations;
            break;
            
        }
        centroids = newCentroids;
        // cout << "IT: " << iterations << endl;
        iterations++;
    }
    return centroids;

}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    double start, end;
    start = MPI_Wtime();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // vector<vector<float>> data = readCsvFile("./archive/diabetes_binary_5050split_health_indicators_BRFSS2015.csv");
    vector<vector<float>> data = readCsvFile("./archive/diabetes_012_health_indicators_BRFSS2015.csv");
    vector<dataPoint> points;

    int n_rows = data.size();
    int n_cols = data[0].size();
    int block_size = n_rows * n_cols / size;


    vector<float> flattened_data(n_rows * n_cols);
    int index = 0;
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            flattened_data[index++] = data[i][j];
        }
    }

    vector<float> block(block_size);
    MPI_Scatter(flattened_data.data(), block_size, MPI_FLOAT, block.data(), block_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Receive the block of data
    vector<vector<float>> received_data(n_rows / size, vector<float>(n_cols));
    index = 0;
    for (int i = 0; i < n_rows / size; i++) {
        for (int j = 0; j < n_cols; j++) {
            received_data[i][j] = block[index++];
        }
    }
    cout << endl << "received_data Size: "<< received_data[0].size() << endl;

    for (int i = 0; i < received_data.size(); i++)
    {
        dataPoint point = {-1, {}};
        for (int j = 1; j < received_data[i].size(); j++)
        {
            
            point.features.push_back(received_data[i][j]);
        }
        points.push_back(point);
    }
    cout << "p_id: " << rank << " | " << points.size() << " | " << points[0].features.size() << endl;
    // cout << endl;

    int k = 3; // Number of clusters
    long maxIterations = 100000; // Maximum number of iterations for k-means
    vector<dataPoint> centroids = kMeans(points, k, maxIterations);
    
//     cout << endl << "p_id: " << rank << ": " << endl;
//     for (int i = 0; i < centroids.size(); i++)
//     {
//         cout << "Centroid " << i << ": " << endl;
//         for (int j = 0; j < centroids[i].features.size(); j++)
//         {
//              cout << centroids[i].features[j] << " ";
//         }
//         cout << endl;
//     }

   
    int n_features = centroids[0].features.size();
    cout << endl << "Num Features: " << n_features;
    vector<float> fltt_data(k*n_features);
    index = 0;
    for (int i = 0; i < centroids.size(); i++) {
        for (int j = 0; j < centroids[i].features.size(); j++){
            fltt_data[index++] = centroids[i].features[j];
        }
        
    }
    // MPI_Barrier( MPI_COMM_WORLD );
    block_size = k*n_features*size;
    vector<float> blockRC(block_size);

    MPI_Gather(fltt_data.data(), k*n_features, MPI_FLOAT, blockRC.data(), k*n_features, MPI_FLOAT, 0, MPI_COMM_WORLD); 

    if (rank == 0){

        cout << endl << endl;
        vector<vector<float>> allPCentroids(size*k, vector<float>(n_features));
        index = 0;
        for (int i = 0; i < size * k; i++) {
            for (int j = 0; j < n_features; j++) {
                allPCentroids[i][j] = blockRC[index++];
                // cout << allPCentroids[i][j] << " ";
            }
            cout << endl;
        }

        double start_merge, end_merge;
        start_merge = MPI_Wtime();
        // Loop até set ter tamanho K
        while (allPCentroids.size() > k) {
            
            float min_dist = numeric_limits<float>::max();
            
            int closest1 = 0, closest2 = 0;

            // Loop pelos pares de pontos
            for (int i = 0; i < allPCentroids.size(); i++) {
                for (int j = i + 1; j < allPCentroids.size(); j++) {
                    // Calculo da distancia entre os pontos
                    float dist = distance(allPCentroids[i].data(), allPCentroids[j].data(), n_features);
                    // Update distancia minima e pontos mais proximo
                    if (dist < min_dist) {
                        min_dist = dist;
                        closest1 = i;
                        closest2 = j;
                    }
                }
            }

            // clacula ponto medio entre pontos mais proximos
            vector<float> midpoint(n_features);
            for (int i = 0; i < n_features; i++) {
                midpoint[i] = (allPCentroids[closest1][i] + allPCentroids[closest2][i]) / 2;
            }

            // Remove os pontos mais proximos
            allPCentroids.erase(allPCentroids.begin() + max(closest1, closest2));
            allPCentroids.erase(allPCentroids.begin() + min(closest1, closest2));

            // Adiciona o ponto médio entre os pontos
            allPCentroids.push_back(midpoint);
        }
        
        end_merge = MPI_Wtime();

        for (int i = 0; i < allPCentroids.size(); i++)
        {
            for (int j = 0; j < allPCentroids[i].size(); j++)
            {
                cout << allPCentroids[i][j] << " ";
            }
            cout << endl;
            
        }
        end = MPI_Wtime();
        cout << endl << "Exec. Time: " << end - start << endl;
        cout << endl << "Merge Exec. Time: " << end_merge - start_merge << endl;
    }
    

    
    MPI_Finalize();
    

    return 0;
} 
