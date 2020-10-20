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

    std::vector<double> yData;
    std::vector< std::vector<double>> xData;
    std::vector< std::vector<double>> coefData;
    std::vector<double> rhoData;


    size_t testN = 272;  // how many sample for predicting
    size_t M = 1000;
    size_t nr_class = 11;
    //Test Data
    ReadMatrixFile(xData, "../data/X_test", testN, M);
    ReadVectorFile(yData,"../data/Y_pred", testN);

    // Precomputed Trained Data
    ReadMatrixFile(coefData, "../data/coef",nr_class,M);
    ReadVectorFile(rhoData, "../data/rho", nr_class);

    //// KEY GENERATION
    TIC(t);
//    double scalingFactor = 5e-1;
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 50;

    SecurityLevel securityLevel = HEStd_128_classic;
    // 0 means the library will choose it based on securityLevel
    uint32_t m = 8192;
    uint32_t batchSize = m/4;

//    CryptoContext<DCRTPoly> cc =
//            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
//                    multDepth, scaleFactorBits, batchSize, securityLevel, m);
    CryptoContext<DCRTPoly> cc =
            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(multDepth,
                                                                 scaleFactorBits,
                                                                 batchSize,
                                                                 securityLevel,
                                                                 m,
                                                                 APPROXAUTO,
                                                                 BV,
                                                                 3,
                                                                 2,
                                                                 50,
                                                                 10,
                                                                 OPTIMIZED);

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


//#pragma omp parallel for
    for (size_t i=0; i<testN; i++){
        for (size_t x=0; x < sizeF; x++){
            Plaintext xTemp = cc->MakeCKKSPackedPlaintext(xDataArray[x][i]);
            X[x][i] = cc->Encrypt(keyPair.publicKey, xTemp);
        }
    }
    cout << "X is encrypted." << endl;
//#pragma omp parallel for
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

//#pragma omp parallel for
    for (size_t s = 0; s < 10; ++s) {

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
//#pragma omp parallel for
    for (size_t s = 0; s < 10; ++s) {
        vector<Plaintext> result_dec_values(nr_class);
        vector<double> temp;
        for (size_t i = 0; i < nr_class; ++i) {
            cc->Decrypt(keyPair.secretKey, dec_values[s][i], &result_dec_values[i]);
            result_dec_values[i]->SetLength(1);
            auto value = result_dec_values[i]->GetCKKSPackedValue();
            temp.push_back(value[0].real());
            cout << temp[i] << endl;
        }
        auto it = max_element(temp.begin(), temp.end());
        cout << "Max Element = " << *it << endl;
        cout << "s= "<< s <<", Predicted in Palisade: " << it - temp.begin() << ", Predicted in Python: "<< yData[s] << endl;
        temp.clear();
    }
    decryptionTime = TOC(t);
    cout << "\nKey Generation Time: \t\t" << keyGenTime/1000 << " s" << endl;
    cout << "Encoding and Encryption Time: \t" << encryptionTime/1000 << " s" << endl;
    cout << "Computation Time: \t\t" << computationTime/1000 << " s" << endl;
    cout << "Decryption & Decoding Time: \t" << decryptionTime/1000 << " s" << endl;

    endToEndTime = TOC(tAll);
    cout << "\nEnd-to-end Runtime: \t\t" << endToEndTime/1000 << " s" << endl;


    return 0;
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

    for(size_t j = 1; j < vector.size(); j=j*2) {
        for(size_t i = 0; i<vector.size(); i = i + 2*j) {
            if ((i+j)<vector.size())
                vector[i] = cc->EvalAdd(vector[i],vector[i+j]);
        }
    }

    Ciphertext<DCRTPoly> result(new CiphertextImpl<DCRTPoly>(*(vector[0])));

    return result;

}