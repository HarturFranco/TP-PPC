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
    // for (int i = 0; i < dataPoints.size(); i++){
    //     float minDist = numeric_limits<float>::max();
    //     int closestCentroid= -1;
    //     for (int j = 0; j < centroids.size(); j++){
    //         float dist = distance(dataPoints[i], centroids[j]);
    //         if (dist < minDist){
    //             minDist = dist;
    //             closestCentroid = i;
    //         }
    //     }
    //     dataPoints[i].cluster = closestCentroid;

    // }
    
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
            break;
            cout << "convergiu";
        }
        centroids = newCentroids;
        iterations++;
    }
    return centroids;

}


int main() {
    vector<vector<float>> data = readCsvFile("./archive/diabetes_012_health_indicators_BRFSS2015.csv");
    vector<dataPoint> points;

    
    for (int i = 0; i < data.size(); i++)
    {
        dataPoint point = {-1, {}};
        for (int j = 1; j < data[i].size(); j++)
        {
            
            point.features.push_back(data[i][j]);
        }
        points.push_back(point);
    }
    
    
    int k = 3; // Number of clusters
    long maxIterations = 99999999999999999; // Maximum number of iterations for k-means
    vector<dataPoint> centroids = kMeans(points, k, maxIterations);
    

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