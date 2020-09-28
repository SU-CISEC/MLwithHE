//
// Created by ferhatyaman on 25.09.2020.
//

#include <fstream>

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
void ReadXFile(std::vector<std::vector<double>> & dataColumns,
                 string dataFileName, size_t N, size_t M)
{

    uint32_t cols = 0;

    string fileName = dataFileName + ".csv";

    std::cerr << "file name = " << fileName << std::endl;

    ifstream file(fileName);
    string line, value;


    size_t counter = 0;
    while((file.good()) && (counter < N)) {
        getline(file, line);
        uint32_t curCols = std::count(line.begin(), line.end(), ',');
        std::vector<double> row(cols);
        for(uint32_t i = 0; i < cols; i++) {
            string substr;
            getline(ss, substr, ',');
            double val;
            val = std::stod(substr);
            row[i] = val;
        }
        dataColumns.push_back(row);
        }
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

void ReadModelFile(std::vector<std::vector<double>> & dataColumns, std::vector<double> &y,
                  string dataFileName, size_t num_class, size_t M)
{

    uint32_t cols = 0;

    string fileName = dataFileName + ".csv";

    std::cerr << "file name = " << fileName << std::endl;

    ifstream file(fileName);
    string line, value;

    if(file.good()) {

        getline(file, line);
        cols = std::count(line.begin(), line.end(), ',') + 1;
        stringstream ss(line);
        vector<string> result;

        size_t tempcounter = 0;

        for(uint32_t i = 0; i < cols; i++) {
            string substr;
            getline(ss, substr, ',');
            if ((substr != "") && (i>4) && (i<M+5)) {
                headers.push_back(substr);
                tempcounter++;
            }
        }

        cols = tempcounter;

    }

    size_t counter = 0;
    while((file.good()) && (counter < N)) {
        getline(file, line);
        uint32_t curCols = std::count(line.begin(), line.end(), ',') + 1;
        if (curCols > 2) {
            stringstream ss(line);
            for(uint32_t i = 0; i < 5; i++) {
                string substr;
                getline(ss, substr, ',');
                if ((i==1))
                    y.push_back(std::stod(substr));
            }
            std::vector<double> row(cols);
            for(uint32_t i = 5; i < cols + 5; i++) {
                string substr;
                getline(ss, substr, ',');
                if (i < M+5)
                {
                    double val;
                    val = std::stod(substr);
                    row[i-5] = val;
                }
            }
            dataColumns.push_back(row);
        }
        counter++;
    }

    file.close();

    std::cout << "Read in data: ";
    std::cout << dataFileName << std::endl;
}
