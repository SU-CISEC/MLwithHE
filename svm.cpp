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

void ReadMatrixFile(vector<vector<double>> & dataColumns,
                    string dataFileName, size_t N, size_t M)
{

    string fileName = dataFileName + ".csv";

    cerr << "file name = " << fileName << endl;

    ifstream file(fileName);
    string line, value;

    size_t counter = 0;
    while((file.good()) && (counter < N)) {
        getline(file, line);
        stringstream ss(line);
        vector<double> row(M);
        for(uint32_t i = 0; i < M; i++) {
            string substr;
            getline(ss, substr, ',');
            double val;
            val = stod(substr);
            row[i] = val;
        }
        dataColumns.push_back(row);
        counter++;
    }

    file.close();

    cout << "Read in data: ";
    cout << dataFileName << endl;
}

template <class T>
void ReadVectorFile(vector<T> & dataColumn,
                    string dataFileName, size_t M)
{

    string fileName = dataFileName + ".csv";

    cerr << "file name = " << fileName << endl;

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

    cout << "Read in data: ";
    cout << dataFileName << endl;
}

double dot(vector<double> px, vector<double> py) {
    double sum = 0;
    size_t i = 0;
    size_t j = 0;
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

    vector<double> yData;
    vector< vector<double>> xData;

    vector< vector<double>> coefData;
    vector<double> Rho;

    size_t testN = 2713;  // how many sample for predicting
    size_t M = 141;
    size_t nr_class = 11;



    //Test Data
    ReadMatrixFile(xData, "../data/SampleData/X_test", testN, M);
    ReadVectorFile(yData,"../data/SampleData/Y_test", testN);

    // Precomputed Trained Data
    ReadMatrixFile(coefData, "../model/coef",nr_class,M);
    ReadVectorFile(Rho,"../model/rho",nr_class);

    int count_comp = 0;
    for (size_t j = 0; j < testN; ++j) {
        vector<double> dec_values(nr_class);
        for(size_t i = 0; i<nr_class; i++){
            dec_values[i] = dot(xData[j],coefData[i]);
            dec_values[i] += Rho[i];
//            cout << "Dec_value["<< i << "]= " << dec_values[i] << endl;
        }
        auto it = max_element(dec_values.begin(), dec_values.end());
//        cout << "MaxElement= " << *it <<" and its Label=" << << endl;
//        cout << "i= "<< j <<", Predicted in C++: " << it - dec_values.begin() << ", Predicted in Python: "<< yData[j] << endl;
        if((it - dec_values.begin()) == yData[j])
            count_comp +=1;
        dec_values.clear();
    }
    cout << "Count in Palisade vs Python: "<< count_comp << endl;

//    for(int i = 0; i<nr_class; i++)
//        cout << "Decision Value["<< i <<"]= " << dec_values[i] << endl;


    return 0;
}