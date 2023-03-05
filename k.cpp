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

    // Skip the first line (headers)
    getline(infile, line);

    // Loop through each line of the file
    while (getline(infile, line)) {
        istringstream iss(line);
        string token;

        

        // Create a vector of floats and fill it with the values from the line
        vector<float> row;
        while (getline(iss, token, ',')) {
            row.push_back(stof(token));
        }

        // Add the vector to the data vector
        data.push_back(row);
    }

    return data;
}


struct dataPoint
{
    int cluster;
    vector<float> features;

};

// Euclidean distance function
float distance(const dataPoint& p1, const dataPoint& p2) {
    float sum = 0;
    for (int i = 0; i < p1.features.size(); i++) {
        sum += pow(p1.features[i] - p2.features[i], 2);
    }
    return sqrt(sum);
}

void assignClusters(vector<dataPoint>& dataPoints, vector<dataPoint>& centroids) {
    
    for (dataPoint& point : dataPoints) {
        float minDist = numeric_limits<float>::max();
        int cluster = 0;
        for (int i = 0; i < centroids.size(); i++) {
            float dist = distance(point, centroids[i]);
            if (dist < minDist) {
                minDist = dist;
                cluster = i;
            }
        }
        point.cluster = cluster;
    }
}

// Function to calculate centroids of clusters
vector<dataPoint> calculateCentroids(vector<dataPoint>& dataPoints, int k) {
    vector<dataPoint> centroids(k);
    vector<int> clusterSizes(k);
    
    for(int i = 0; i < centroids.size(); i++){
        centroids[i].features = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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

    int iterations = 0;
    while (iterations < maxIterations) {
    // while (true) {
        assignClusters(dataPoints, centroids);
        vector<dataPoint> newCentroids = calculateCentroids(dataPoints, k);
        // If centroids have not moved, terminate early
        bool convergence = true;
        for (int i = 0; i < k; i++) {
            if (distance(newCentroids[i], centroids[i]) > 0.0001) {
                convergence = false;
                break;
            }
        }
        if (convergence) {
            cout << "convergiu: " << iterations;
            break;
            
        }
        centroids = newCentroids;
        cout << "IT: " << iterations << endl;
        iterations++;
    }
    return centroids;

}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<vector<float>> data = readCsvFile("./archive/diabetes_012_health_indicators_BRFSS2015.csv");
    vector<dataPoint> points;

    int n_rows = data.size();
    int n_cols = data[0].size();
    int block_size = n_rows * n_cols / size;


    float array[20000][22];
    vector<float> flattened_data(n_rows * n_cols);
    int index = 0;
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            flattened_data[index++] = data[i][j];
            // cout << flattened_data[index -1] << " ";
        }
    }
    // cout << "size: " <<flattened_data.size() << endl;

    std::vector<float> block(block_size);
    MPI_Scatter(flattened_data.data(), block_size, MPI_FLOAT, block.data(), block_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Receive the block of data
    std::vector<std::vector<float>> received_data(n_rows / size, std::vector<float>(n_cols));
    index = 0;
    for (int i = 0; i < n_rows / size; i++) {
        for (int j = 0; j < n_cols; j++) {
            received_data[i][j] = block[index++];
        }
    }

    // cout << "p_id: " << rank << " | " << received_data.size() << endl;
    

    for (int i = 0; i < received_data.size(); i++)
    {
        dataPoint point = {-1, {}};
        for (int j = 1; j < received_data[i].size(); j++)
        {
            
            point.features.push_back(received_data[i][j]);
        }
        points.push_back(point);
    }
    cout << "p_id: " << rank << " | " << points.size() << endl;
    // cout << endl;

    int k = 3; // Number of clusters
    long maxIterations = 100000; // Maximum number of iterations for k-means
    vector<dataPoint> centroids = kMeans(points, k, maxIterations);

    cout << "p_id: " << rank << ": " << endl;
    for (int i = 0; i < centroids.size(); i++)
    {
        cout << "Centroid " << i << ": " << endl;
        for (int j = 0; j < centroids[i].features.size(); j++)
        {
             cout << centroids[i].features[j] << " ";
        }
        cout << endl;
    }
    
    MPI_Finalize();
    // memset(array, 0, sizeof(array));

    // for (int i = 0; i < data.size(); i++)
    // {
    //     dataPoint point = {-1, {}};
    //     for (int j = 1; j < data[i].size(); j++)
    //     {
            
    //         point.features.push_back(data[i][j]);
    //     }
    //     points.push_back(point);
    // }
    
    
    // int k = 3; // Number of clusters
    // long maxIterations = 100000; // Maximum number of iterations for k-means
    // vector<dataPoint> centroids = kMeans(points, k, maxIterations);
    
    
    // for (int i = 0; i < centroids.size(); i++)
    // {
    //     cout << "Centroid " << i << ": " << endl;
    //     for (int j = 0; j < centroids[i].features.size(); j++)
    //     {
    //          cout << centroids[i].features[j] << " ";
    //     }
    //     cout << endl;
    // }
    

    //  for (dataPoint& point : points) {
    //     cout << "Point Cluster " << point.cluster << endl << "Feat: ";
    //     for (int i = 0; i < point.features.size(); i++)
    //     {
    //         cout << point.features[i] << " ";
    //     }
    //     cout << endl;
        
    //  }
    

    // cout << "N Features in point: " << points.size() << endl;


    // cout << "data 0: " << endl;
    // for (int j = 1; j < data[0].size(); j++)
    // {
    //     cout << data[0][j] << " ";
    // }
    // cout << endl << points[0].cluster << " point features: " << endl;
    // for (int j = 0; j < points[0].features.size(); j++)
    // {
    //     cout << points[0].features[j] << " ";
    // }
    // // Print the data
    // for (const auto& row : data) {
    //     for (const auto& value : row) {
    //         cout << value << " ";
    //     }
    //     cout << endl;
    // }

    return 0;
} 