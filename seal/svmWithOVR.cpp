//
// Created by ferhatyaman on 13.10.2020.
//
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

void ReadMatrixFile(std::vector<std::vector<double>> & dataColumns,
                    string dataFileName, size_t N, size_t M)
{

    string fileName = dataFileName + ".csv";

    std::cerr << "file name = " << fileName << std::endl;

    ifstream file(fileName);
    string line, value;

    size_t counter = 0;
    while((file.good()) && (counter < N)) {
        getline(file, line);
//        uint32_t curCols = std::count(line.begin(), line.end(), ',');
        stringstream ss(line);
        std::vector<double> row(M);
        for(uint32_t i = 0; i < M; i++) {
            string substr;
            getline(ss, substr, ',');
            double val;
            val = std::stod(substr);
            row[i] = val;
        }
        dataColumns.push_back(row);
        counter++;
    }

    file.close();

    std::cout << "Read in data: ";
    std::cout << dataFileName << std::endl;
}

/**
 * Read model parameters file values into different matrix or .
 *
 * @param dataColumns Column vector where the data placed as row of double values.
 * @param dataFileName Path of the file
 * @param M Number of features
 * @param num_class Number of classes
 * @return nothing.
 * @note This functions assumes file does not have header line.
 */
template <class T>
void ReadVectorFile(std::vector<T> & dataColumn,
                    string dataFileName, size_t M)
{

    string fileName = dataFileName + ".csv";

    std::cerr << "file name = " << fileName << std::endl;

    ifstream file(fileName);
    string line;

    for(uint32_t i = 0; i < M && file.good(); i++) {
        getline(file, line);
        stringstream ss(line);
        T val;
        ss >> val;
        dataColumn.push_back(val);
    }


    file.close();

    std::cout << "Read in data: ";
    std::cout << dataFileName << std::endl;
}

double dot(vector<double> px, vector<double> py) {
    double sum = 0;
    int i = 0;
    int j = 0;
    while(i != px.size() && j != py.size()) {
        if(i == j) {
            sum += px[i] * py[j];
            ++i;
            ++j;
        }
        else {
            if(i > j)
                ++j;
            else
                ++i;
        }
    }
    return sum;
}
int main(){
    //// DATA READ

    std::vector<double> yData;
    std::vector< std::vector<double>> xData;

    std::vector< std::vector<double>> coefData;
    std::vector<double> Rho;

    size_t testN = 272;  // how many sample for predicting
    size_t M = 1000;
    size_t nr_class = 11;
    size_t dual_class = (nr_class * (nr_class - 1))/2;

    //Test Data
    ReadMatrixFile(xData, "../data/X_test",testN,M);
    ReadVectorFile(yData,"../data/Y_pred", testN);

    // Precomputed Trained Data
    ReadMatrixFile(coefData, "../data/coef", nr_class, M);
    ReadVectorFile(Rho,"../data/rho",nr_class);


    for (int j = 0; j < testN; ++j) {
        vector<double> dec_values(nr_class);
        for(int i = 0; i<nr_class; i++){
            dec_values[i] = dot(xData[j],coefData[i]);
            dec_values[i] += Rho[i];
        }
        auto it = max_element(dec_values.begin(), dec_values.end());
//        cout << "MaxElement= " << *it <<" and its Label=" << << endl;
        cout << "i= "<< j <<", Predicted in C++: " << it - dec_values.begin() << ", Predicted in Python: "<< yData[j] << endl;
        dec_values.clear();
    }


//    for(int i = 0; i<nr_class; i++)
//        cout << "Decision Value["<< i <<"]= " << dec_values[i] << endl;


    return 0;
}