//
// Created by Ferhat Yaman on 26.09.2020.
//

#include "../utils.h"

void assignDataToSlots(vector<vector<std::vector<std::complex<double>>>> &arrayData, std::vector<std::vector<double>> data,
                       size_t n_array, usint m);

int main(int argc, char **argv) {


    cout << "\n======SVM PALISADE Solution========\n" << std::endl;

    //// DATA READ

    std::vector<double> yData;
    std::vector< std::vector<double>> xData;
    std::vector< std::vector<double>> svData;
    std::vector< std::vector<double>> svCoefData;

    string xFilename = "../data/XScaled";
    string yFilename = "../data/Y";
    string svFilename = "../data/SV";
    string svCoefFilename = "../data/sv_coef";

    size_t trainN = 270; // how many sample used for training
    size_t testN = 270;  // how many sample for predicting
    size_t M = 1000;
    size_t nr_class = 10;

    //Test Data
    ReadMatrixFile(xData, xFilename,testN,M);
    ReadVectorFile(yData,yFilename, testN);

    // Precomputed Trained Data
    ReadMatrixFile(svData, svFilename,trainN,M);
    ReadMatrixFile(svCoefData, svCoefFilename, nr_class, trainN);

    //// KEY GENERATION
//    double scalingFactor = 5e-1;
    usint m = 16384;

    usint init_size = 4;
    usint dcrtBits = 54;
    // TODO Check Parameters again
    CryptoContext<DCRTPoly> cc =
            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
                    init_size-1,
                    dcrtBits,
                    0,
                    HEStd_128_classic,
                    m/2, /*ringDimension*/
                    APPROXRESCALE,
                    BV,
                    3, /*numLargeDigits*/
                    2, /*maxDepth*/
                    dcrtBits, /*firstMod*/
                    0,
                    OPTIMIZED);

    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    cc->Enable(LEVELEDSHE);

    std::cout << "\nNumber of Samples = " << xData.size() << std::endl;
    std::cout << "Number of Features = " << xData[0].size() << std::endl;


    auto keyPair = cc->KeyGen();

    //// ENCODE and ENCRYPTION
    size_t sizeF = (size_t)std::ceil((double)xData[0].size() / (m / 4));
    size_t sizeC = (size_t)std::ceil((double)svCoefData[0].size() / (m / 4));

    std::vector<std::vector<std::vector<std::complex<double>>>> xDataArray(sizeF);
    std::vector<std::vector<std::vector<std::complex<double>>>> svDataArray(sizeF);
    std::vector<std::vector<std::vector<std::complex<double>>>> svCoefDataArray(sizeC);

    // Divide data into matrix to fit inside to slot size = n/2
    assignDataToSlots(xDataArray, xData , sizeF, m);
    assignDataToSlots(svDataArray, svData , sizeF, m);
    assignDataToSlots(svCoefDataArray, svCoefData , sizeF, m);


    std::vector<std::vector<Ciphertext<DCRTPoly>>> X(sizeF);
    std::vector<std::vector<Ciphertext<DCRTPoly>>> SV(sizeF);
    std::vector<std::vector<Ciphertext<DCRTPoly>>> C(sizeC);
    std::vector<Ciphertext<DCRTPoly>> Y(testN);

    for (size_t i = 0; i < sizeF; i++) {
        X[i] = std::vector<Ciphertext<DCRTPoly>>(testN);
        SV[i] = std::vector<Ciphertext<DCRTPoly>>(trainN);
    }
    for (size_t i = 0; i < sizeC; ++i) {
        C[i] = std::vector<Ciphertext<DCRTPoly>>(nr_class);
    }


//#pragma omp parallel for
    for (size_t i=0; i<testN; i++){
        for (size_t x=0; x < sizeF; x++){
            Plaintext xTemp = cc->MakeCKKSPackedPlaintext(xDataArray[x][i]);
            X[x][i] = cc->Encrypt(keyPair.publicKey, xTemp);
        }
        Plaintext sTemp2 = cc->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(m/4,yData[i]));
        Y[i] = cc->Encrypt(keyPair.publicKey, sTemp2);
    }
    cout << "X and Y encrypted." << endl;
//#pragma omp parallel for
    for (size_t i = 0; i < trainN; ++i) {
        for (size_t x=0; x < sizeF; x++){
            Plaintext svTemp = cc->MakeCKKSPackedPlaintext(svDataArray[x][i]);
            SV[x][i] = cc->Encrypt(keyPair.publicKey, svTemp);
        }
    }
    cout << "SV encrypted." << endl;

//#pragma omp parallel for
    for (size_t i = 0; i < nr_class; ++i) {
        for (size_t x=0; x < sizeC; x++){
            Plaintext cTemp = cc->MakeCKKSPackedPlaintext(svCoefDataArray[x][i]);
            C[x][i] = cc->Encrypt(keyPair.publicKey, cTemp);
        }
    }
    cout << "C encrypted." << endl;
    //// SVM OPERATIONS




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
