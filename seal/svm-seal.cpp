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

void ReadVectorFile(std::vector<double> & dataColumn,
                    string dataFileName, size_t M)
{

    string fileName = dataFileName + ".csv";

    std::cerr << "file name = " << fileName << std::endl;

    ifstream file(fileName);
    string line;

    for(uint32_t i = 0; i < M && file.good(); i++) {
        getline(file, line);
        stringstream ss(line);
        double val;
        ss >> val;
        dataColumn.push_back(val);
    }


    file.close();

    std::cout << "Read in data: ";
    std::cout << dataFileName << std::endl;
}



int main(){
    //// DATA READ

    std::vector<double> yData;
    std::vector< std::vector<double>> xData;
    std::vector< std::vector<double>> svData;
    std::vector< std::vector<double>> svCoefData;

    string xFilename = "../data/XScaled";
    string yFilename = "../data/Y";
    string svFilename = "../data/SV";
    string svCoefFilename = "../data/sv_coef";

    // TODO: Define N for test and train separately
    size_t N = 270;
    size_t M = 1000;
    size_t n_classes = 10;
    ReadMatrixFile(xData, xFilename,N,M);
    ReadVectorFile(yData,yFilename, N);
    ReadMatrixFile(svData, svFilename,N,M);
    ReadMatrixFile(svCoefData, svCoefFilename,n_classes,N);
    return 0;
}