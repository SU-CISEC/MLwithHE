//
// Created by Ferhat Yaman on 26.09.2020.
//

#include <fstream>
#include "palisade.h"
#include <cmath>

using namespace std;
using namespace lbcrypto;

/**
 * Read X test file values into a matrix.
 *
 * @param dataColumns Column vector where the data placed as row of double values.
 * @param patientList Column vector where the data placed as row of string value.
 * @param dataFileName Path of the file
 * @param M Number of features
 * @return nothing.
 * @note This functions assumes file does not have header line.
 */
void ReadTestData(vector<vector<complex<double>>> & dataColumns, vector<string> & patientList,string dataFileName, size_t M)
{
    string fileName = dataFileName + ".txt";

    cerr << "file name = " << fileName << endl;

    ifstream file(fileName);
    string line;

    while(file.good() && getline(file, line)) {
        stringstream ss(line);

        string substr;
        getline(ss, substr, ' ');
        patientList.push_back(substr);

        vector<complex<double>> row(M);
        for(uint32_t i = 0; i < M; i++) {
            getline(ss, substr, ' ');
            double val;
            val = stod(substr);
            row[i] = val;
        }
        dataColumns.push_back(row);
    }

    file.close();

    cout << "Read in data: ";
    cout << dataFileName << endl;
}

/**
 * Read Coef model values into a matrix.
 *
 * @param dataColumns Column vector where the data placed as row of double values.
 * @param dataFileName Path of the file
 * @param M Number of features
 * @param N Number of classes
 * @return nothing.
 * @note This functions assumes file does not have header line.
 */
void ReadMatrixFile(vector<vector<complex<double>>> & dataColumns, string dataFileName, size_t N, size_t M)
{
    string fileName = dataFileName + ".csv";

    ifstream file(fileName);
    string line, value;

    size_t counter = 0;
    while((file.good()) && (counter < N)) {
        getline(file, line);
        stringstream ss(line);
        vector<complex<double>> row(M);
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

/**
 * Read model parameters file rho values into column vector.
 *
 * @param dataColumns Column vector where the data placed as row of double value.
 * @param dataFileName Path of the file
 * @param M Number of features
 * @param num_class Number of classes
 * @return nothing.
 * @note This functions assumes file does not have header line.
 */
template <class T>
void ReadVectorFile(vector<T> & dataColumn, string dataFileName, size_t M)
{
    string fileName = dataFileName + ".csv";

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

/**
 * Print out result values into csv file.
 *
 * @param dataColumns Column vector which the data placed as row of double values.
 * @param patientList Column vector where the data places as row of string value.
 * @return nothing.
 * @note This functions assumes file does not have header line.
 */
void print_results(vector<vector<double>> res, vector<string> patientList) {
    string fileName =  "../data/results.csv";

    uint32_t N = res.size();
    uint32_t nr_class = res[0].size();
    size_t counter = 0;

    ofstream file(fileName);
    file.precision(8);

    file << "Patient-Id,Bronchusandlung,Bladder,Colon,Skin,Stomach,Corpusuteri,Liverandintrahepaticbileducts,Ovary,Kidney,Cervixuteri,Breast\n";
    while((file.good()) && (counter < N)) {
        file << patientList[counter] << ',';

        for(uint32_t i = 0; i < nr_class-1; i++)
            file << res[counter][i] << ',';

        file << res[counter][nr_class-1] << "\n";
        counter++;
    }

    file.close();

    cout << "Print out data: " << fileName << endl;
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


    cout << "\n======SVM PALISADE Solution========\n" << endl;

    //// DATA READ
    vector<string> patient_Ids;
    vector<vector<complex<double>>> xData;
    vector<vector<complex<double>>> coefData;
    vector<double> rhoData;

    size_t M = 141;
    size_t nr_class = 11;
    //Test Data
    ReadTestData(xData, patient_Ids, "../data/data.test", M);

    // Precomputed Trained Data
    ReadMatrixFile(coefData, "../model/coef",nr_class,M);
    ReadVectorFile(rhoData, "../model/rho", nr_class);

    size_t testN = xData.size();  // how many sample for predicting

//// KEY GENERATION
    TIC(t);
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 20;

    SecurityLevel securityLevel = HEStd_128_classic;
    // 0 means the library will choose it based on securityLevel
    uint32_t m = 2048;
    uint32_t batchSize = 141;

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

    cout << "\nNumber of Samples = " << xData.size() << endl;
    cout << "Number of Features = " << xData[0].size() << endl;


    auto keyPair = cc->KeyGen();
    cc->EvalMultKeysGen(keyPair.secretKey);
    cc->EvalSumKeyGen(keyPair.secretKey);
    keyGenTime = TOC(t);

    //// ENCODE and ENCRYPTION
    TIC(t);


    vector<Ciphertext<DCRTPoly>> X(testN);
    vector<Plaintext> W(nr_class);
    vector<Plaintext> B(nr_class);


#pragma omp parallel for
    for (size_t i=0; i<testN; i++){
        Plaintext xTemp = cc->MakeCKKSPackedPlaintext(xData[i]);
        X[i] = cc->Encrypt(keyPair.publicKey, xTemp);
    }
    cout << "X is encrypted." << endl;
#pragma omp parallel for
    for (size_t i = 0; i < nr_class; ++i) {
        W[i] = cc->MakeCKKSPackedPlaintext(coefData[i]);
        W[i]->SetFormat(EVALUATION);
        B[i] = cc->MakeCKKSPackedPlaintext(vector<complex<double>>(nr_class,rhoData[i]));
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
            auto temp = cc->EvalMult(X[s],W[i]);
            temp = cc->EvalSum(temp,M);
            dec_values[s][i] = cc->Rescale(temp);
            dec_values[s][i] = cc->EvalAdd(dec_values[s][i], rhoData[i]);
        }
    }

    computationTime = TOC(t);
    //// DECRYPTION
    TIC(t);
    vector<vector<double>> results(testN);
    for (uint32_t i = 0; i < testN; ++i)
        results[i] = vector<double>(nr_class);

#pragma omp parallel for
    for (size_t s = 0; s < testN; ++s) {
        Plaintext result_dec_values;
        for (size_t i = 0; i < nr_class; ++i) {
            cc->Decrypt(keyPair.secretKey, dec_values[s][i], &result_dec_values);
            result_dec_values->SetLength(1);
            auto value = result_dec_values->GetCKKSPackedValue();
            results[s][i] = value[0].real();
        }
    }
    decryptionTime = TOC(t);

    print_results(results, patient_Ids);
    cout << "\nKey Generation Time: \t\t" << keyGenTime/1000 << " s" << endl;
    cout << "Encoding and Encryption Time: \t" << encryptionTime/1000 << " s" << endl;
    cout << "Computation Time: \t\t" << computationTime/1000 << " s" << endl;
    cout << "Decryption & Decoding Time: \t" << decryptionTime/1000 << " s" << endl;

    endToEndTime = TOC(tAll);
    cout << "\nEnd-to-end Runtime: \t\t" << endToEndTime/1000 << " s" << endl;


    return 0;
}
