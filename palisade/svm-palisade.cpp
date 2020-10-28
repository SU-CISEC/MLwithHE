//
// Created by Ferhat Yaman on 26.09.2020.
//

#include <fstream>
#include "palisade.h"
#include <cmath>
#include <algorithm>

using namespace std;
using namespace lbcrypto;


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

void assignDataToSlots(vector<vector<std::vector<std::complex<double>>>> &arrayData, std::vector<std::vector<double>> data,
                       size_t n_array, usint m);
Ciphertext<DCRTPoly> BinaryTreeAdd(std::vector<Ciphertext<DCRTPoly>> &vector);

void print_results(vector<vector<double>> res);

void softmax(vector<double> & temp){
    int sum = 0;
    for (size_t i = 0; i < temp.size(); ++i)
        sum += exp(temp[i]);
    for (size_t i = 0; i < temp.size(); ++i)
        temp[i] =  exp(temp[i])/ sum;
}

int main(int argc, char **argv) {

    TimeVar t;
    TimeVar tAll;

    TIC(tAll);

    double keyGenTime(0.0);
    double encryptionTime(0.0);
    double computationTime(0.0);
    double decryptionTime(0.0);
    double endToEndTime(0.0);


    cout << "\n======SVM PALISADE Solution========\n" << std::endl;

    //// DATA READ

    std::vector<double> yData,yTrue;
    std::vector< std::vector<double>> xData;
    std::vector< std::vector<double>> coefData;
    std::vector<double> rhoData;


    size_t testN = 2713;  // how many sample for predicting
    size_t M = 141;
    size_t nr_class = 11;
    //Test Data
    ReadMatrixFile(xData, "../data/feature141/X_test", testN, M);
    ReadVectorFile(yData,"../data/finalData/Y_pred", testN);
    ReadVectorFile(yTrue,"../data/feature141/Y_test", testN);

    // Precomputed Trained Data
    ReadMatrixFile(coefData, "../data/feature141/coef",nr_class,M);
    ReadVectorFile(rhoData, "../data/feature141/rho", nr_class);

    //// KEY GENERATION
    TIC(t);
//    double scalingFactor = 5e-1;
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 20;

    SecurityLevel securityLevel = HEStd_128_classic;
    // 0 means the library will choose it based on securityLevel
    uint32_t m = 2048;
    uint32_t batchSize = 256;

//    CryptoContext<DCRTPoly> cc =
//            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
//                    multDepth, scaleFactorBits, batchSize, securityLevel, m);
    CryptoContext<DCRTPoly> cc =
            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(multDepth,
                                                                 scaleFactorBits,
                                                                 batchSize,
                                                                 securityLevel,
                                                                 m,
                                                                 EXACTRESCALE,
                                                                 BV,
                                                                 1,
                                                                 2,
                                                                 20,
                                                                 0,
                                                                 SPARSE);

    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    cc->Enable(LEVELEDSHE);

    std::cout << "\nNumber of Samples = " << xData.size() << std::endl;
    std::cout << "Number of Features = " << xData[0].size() << std::endl;


    auto keyPair = cc->KeyGen();
    cc->EvalMultKeysGen(keyPair.secretKey);
    cc->EvalSumKeyGen(keyPair.secretKey);
    keyGenTime = TOC(t);
    //// ENCODE and ENCRYPTION
    TIC(t);
    size_t sizeF = (size_t)std::ceil((double)xData[0].size() / (m / 4));

    std::vector<std::vector<std::vector<std::complex<double>>>> xDataArray(sizeF);
    std::vector<std::vector<std::vector<std::complex<double>>>> coefDataArray(sizeF);

    // Divide data into matrix to fit inside to slot size = n/2
    assignDataToSlots(xDataArray, xData , sizeF, m);
    assignDataToSlots(coefDataArray, coefData , sizeF, m);


    std::vector<std::vector<Ciphertext<DCRTPoly>>> X(sizeF);
    std::vector<std::vector<Plaintext>> W(sizeF);
    std::vector<Plaintext> B(nr_class);

    for (size_t i = 0; i < sizeF; i++) {
        X[i] = std::vector<Ciphertext<DCRTPoly>>(testN);
        W[i] = std::vector<Plaintext>(nr_class);
    }


#pragma omp parallel for
    for (size_t i=0; i<testN; i++){
        for (size_t x=0; x < sizeF; x++){
            Plaintext xTemp = cc->MakeCKKSPackedPlaintext(xDataArray[x][i]);
            X[x][i] = cc->Encrypt(keyPair.publicKey, xTemp);
        }
    }
    cout << "X is encrypted." << endl;
#pragma omp parallel for
    for (size_t i = 0; i < nr_class; ++i) {
        for (size_t x=0; x < sizeF; x++){
            W[x][i] = cc->MakeCKKSPackedPlaintext(coefDataArray[x][i]);
            W[x][i]->SetFormat(EVALUATION);
        }
        B[i] = cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(nr_class,rhoData[i]));
    }
    cout << "W is packed." << endl;
    encryptionTime = TOC(t);

    //// SVM OPERATIONS
    TIC(t);
    vector<vector<Ciphertext<DCRTPoly>>> dec_values(testN);
    for (size_t s = 0; s < testN; ++s) {
        dec_values[s] = vector<Ciphertext<DCRTPoly>>(nr_class);
    }

#pragma omp parallel for
    for (size_t s = 0; s < testN; ++s) {

        for (size_t i = 0; i < nr_class; ++i) {
            std::vector<Ciphertext<DCRTPoly>> temp(sizeF);
            for (size_t x=0; x < sizeF; x++){
                temp[x] = cc->EvalMult(X[x][s],W[x][i]);
                temp[x] = cc->EvalSum(temp[x],M);
            }
//            dec_values[s][i] = BinaryTreeAdd(temp);
            dec_values[s][i] = cc->Rescale(temp[0]);
            temp.clear();
            dec_values[s][i] = cc->EvalAdd(dec_values[s][i], rhoData[i]);
        }
    }

    computationTime = TOC(t);
    //// DECRYPTION
    TIC(t);
    int count_comp = 0;
    int count_python = 0;
    int count_palisade = 0;

    vector<vector<double>> results(testN);
    for (uint32_t i = 0; i < testN; ++i)
        results[i] = vector<double>(nr_class);

#pragma omp parallel for
    for (size_t s = 0; s < testN; ++s) {
        Plaintext result_dec_values;
        vector<double> temp;
        for (size_t i = 0; i < nr_class; ++i) {
            cc->Decrypt(keyPair.secretKey, dec_values[s][i], &result_dec_values);
            result_dec_values->SetLength(1);
            auto value = result_dec_values->GetCKKSPackedValue();
            temp.push_back(value[0].real());
            results[s][i] = value[0].real();
        }
//        softmax(temp);
        auto it = max_element(temp.begin(), temp.end());
//        cout << "Max Element = " << *it << endl;
//        cout << "s= "<< s <<", Predicted in Palisade: " << it - temp.begin() << ", Predicted in Python: "<< yData[s] << ", True Value of Sample: "<< yTrue[s] << endl;
        if((it - temp.begin()) == yData[s])
            count_comp +=1;
        if((it - temp.begin()) == yTrue[s])
            count_palisade +=1;
        if(yData[s] == yTrue[s])
            count_python +=1;
        temp.clear();
    }
    print_results(results);

    cout << "Count in Palisade vs Python: "<< count_comp << endl;
    cout << "Count in Palisade vs True: "<< count_palisade << endl;
    cout << "Count in Python vs True: "<< count_python << endl;

    decryptionTime = TOC(t);
    cout << "\nKey Generation Time: \t\t" << keyGenTime/1000 << " s" << endl;
    cout << "Encoding and Encryption Time: \t" << encryptionTime/1000 << " s" << endl;
    cout << "Computation Time: \t\t" << computationTime/1000 << " s" << endl;
    cout << "Decryption & Decoding Time: \t" << decryptionTime/1000 << " s" << endl;

    endToEndTime = TOC(tAll);
    cout << "\nEnd-to-end Runtime: \t\t" << endToEndTime/1000 << " s" << endl;


    return 0;
}

void print_results(vector<vector<double>> res) {
    string fileName =  "../data/feature141/results.csv";

    std::cerr << "file name = " << fileName << std::endl;

    ofstream file(fileName);
    file.precision(8);
    uint32_t M = res.size();
    uint32_t nr_class = res[0].size();
    size_t counter = 0;
    while((file.good()) && (counter < M)) {
        for(uint32_t i = 0; i < nr_class-1; i++) {
            file << res[counter][i] << ',';
        }
        file << res[counter][nr_class-1] << "\n";
        counter++;
    }

    file.close();

    std::cout << "Read in data: " << fileName << std::endl;
}

void assignDataToSlots(vector<vector<std::vector<std::complex<double>>>> &arrayData, std::vector<std::vector<double>> data,
                       size_t n_array, usint m) {
    for(size_t x = 0; x < n_array; x++)
        arrayData[x] = std::vector<std::vector<std::complex<double>>>(data.size());

    for (size_t i=0; i < data.size(); i++){

        for(size_t x = 0; x < n_array; x++)
            arrayData[x][i] = std::vector<std::complex<double>>(data[i].size());

        size_t counter = 0;

        for (size_t j=0; j<data[i].size(); j++) {
            if ((j>0) && (j%(m/4)==0))
                counter++;
            arrayData[counter][i][j%(m/4)] = data[i][j];
        }
    }

}
Ciphertext<DCRTPoly> BinaryTreeAdd(std::vector<Ciphertext<DCRTPoly>> &vector) {

    auto cc = vector[0]->GetCryptoContext();

    for(size_t j = 1; j < vector.size(); ++j) {
        vector[0] = cc->EvalAdd(vector[0],vector[j]);
    }

    Ciphertext<DCRTPoly> result(new CiphertextImpl<DCRTPoly>(*(vector[0])));

    return result;

}