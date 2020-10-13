//
// Created by phicoder on 28.09.2020.
//
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


/**
 * Read X test file values into a matrix.
 *
 * @param dataColumns Column vector where the data placed as row of double values.
 * @param dataFileName Path of the file
 * @param M Number of features
 * @param N Number of sample size
 * @return nothing.
 * @note This functions assumes file does not have header line.
 */
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

double predict_values(int l, vector<vector<double>> SV,vector<vector<double>> sv_coef, vector<int>nSV, vector<double> rho, vector<double> x, vector<double> &dec_values)
{
    int i;
    int nr_class = 11;
    vector<double> kvalue(l);
    for(i=0;i<l;i++)
        kvalue[i] = dot(x,SV[i]);

    vector<int> start(nr_class);
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+nSV[i-1];

    vector<int> vote(nr_class);
    for(i=0;i<nr_class;i++)
        vote[i] = 0;

    int p=0;
    for(i=0;i<nr_class;i++)
        for(int j=i+1;j<nr_class;j++)
        {
            double sum = 0;
            int si = start[i];
            int sj = start[j];
            int ci = nSV[i];
            int cj = nSV[j];

            int k;
            vector<double> coef1 = sv_coef[j-1];
            vector<double> coef2 = sv_coef[i];
            for(k=0;k<ci;k++)
                sum += coef1[si+k] * kvalue[si+k];
            for(k=0;k<cj;k++)
                sum += coef2[sj+k] * kvalue[sj+k];
            sum -= rho[p];
            dec_values[p] = sum;

            if(dec_values[p] > 0)
                ++vote[i];
            else
                ++vote[j];
            p++;
        }

    int vote_max_idx = 0;
    for(i=1;i<nr_class;i++){
        if(vote[i] > vote[vote_max_idx])
            vote_max_idx = i;
    }
    start.clear();
    kvalue.clear();
    dec_values.clear();
    vote.clear();

    return vote_max_idx;

}

int main(){
    //// DATA READ

    std::vector<double> yData;
    std::vector< std::vector<double>> xData;
    std::vector< std::vector<double>> svData;
    std::vector< std::vector<double>> svCoefData;
    std::vector<double> svRho;

    std::vector<int> svIndex;
    std::vector<int> svClass;

    size_t trainN = 263; // how many sample used for training
    size_t testN = 272;  // how many sample for predicting
    size_t M = 1000;
    size_t nr_class = 11;
    size_t dual_class = (nr_class * (nr_class - 1))/2;

    //Test Data
    ReadMatrixFile(xData, "../data/XScaled",testN,M);
    ReadVectorFile(yData,"../data/Y_test", testN);

    // Precomputed Trained Data
    ReadMatrixFile(svData, "../data/SV",trainN,M);
    ReadMatrixFile(svCoefData, "../data/sv_coef", nr_class-1, trainN);
    ReadVectorFile(svRho,"../data/rho",dual_class);

    //Support Vector Indices for Trained Model
    ReadVectorFile(svIndex,"../data/SV_index",trainN);
    ReadVectorFile(svClass,"../data/SV_class",nr_class);

    vector<double> dec_values(dual_class);
    for(int i = 0; i<testN; i++){
        int predict_in_c = predict_values(trainN, svData,svCoefData,svClass,svRho,xData[i],dec_values);
        dec_values.clear();
        cout << "i= "<< i <<", Predicted in C++: " << predict_in_c << ", Predicted in Python: "<< yData[i] << endl;
    }
    for(int i = 0; i<dual_class; i++)
        cout << "Decision Value["<< i <<"]= " << dec_values[i] << endl;

    return 0;
}