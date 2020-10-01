//
// Created by Ferhat Yaman on 26.09.2020.
//

#include "../utils.h"

void assignDataToSlots(vector<vector<std::vector<std::complex<double>>>> &arrayData, std::vector<std::vector<double>> data,
                       size_t n_array, usint m);
Ciphertext<DCRTPoly> BinaryTreeAdd(std::vector<Ciphertext<DCRTPoly>> &vector);

int main(int argc, char **argv) {


    cout << "\n======SVM PALISADE Solution========\n" << std::endl;

    //// DATA READ

    std::vector<double> yData;
    std::vector< std::vector<double>> xData;
    std::vector< std::vector<double>> svData;
    std::vector< std::vector<double>> svCoefData;

    std::vector<int> svIndex;
    std::vector<int> svClass;

    size_t trainN = 268; // how many sample used for training
    size_t testN = 268;  // how many sample for predicting
    size_t M = 1000;
    size_t nr_class = 11;

    //Test Data
    ReadMatrixFile(xData, "../data/XScaled",testN,M);
    ReadVectorFile(yData,"../data/Y", testN);

    // Precomputed Trained Data
    ReadMatrixFile(svData, "../data/SV",trainN,M);
    ReadMatrixFile(svCoefData, "../data/sv_coef", nr_class-1, trainN);

    //Support Vector Indices for Trained Model
    ReadVectorFile(svIndex,"../data/SV_index",trainN);
    ReadVectorFile(svClass,"../data/SV_class",nr_class);

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
    cc->EvalMultKeysGen(keyPair.secretKey);
    cc->EvalSumKeyGen(keyPair.secretKey);

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
    for (size_t i = 0; i < nr_class-1; ++i) {
        for (size_t x=0; x < sizeC; x++){
            Plaintext cTemp = cc->MakeCKKSPackedPlaintext(svCoefDataArray[x][i]);
            C[x][i] = cc->Encrypt(keyPair.publicKey, cTemp);
        }
    }
    cout << "C encrypted." << endl;
    //// SVM OPERATIONS

    std::vector<Ciphertext<DCRTPoly>> dec_values(nr_class*(nr_class-1)/2);

    //  kvalue[i] = dot(x,model->SV[i],model->param,blas_functions);
    std::vector<Ciphertext<DCRTPoly>> kvalue(trainN);
//#pragma omp parallel for
    for (size_t i = 0; i < trainN; ++i) {
        std::vector<Ciphertext<DCRTPoly>> temp(sizeF);
        for (size_t x=0; x < sizeF; x++){
            temp[i] = cc->EvalMult(X[x][i],SV[x][i]);
        }
        kvalue[i] = BinaryTreeAdd(temp);
        kvalue[i] = cc->ModReduce(kvalue[i]);
    }



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