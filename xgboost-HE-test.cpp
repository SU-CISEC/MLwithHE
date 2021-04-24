#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <utils/debug.h>
#include "palisade.h"

using namespace std;
using namespace lbcrypto;


void cleanInput(string & input) {
    input.erase(remove(input.begin(), input.end(), '['), input.end());
    input.erase(remove(input.begin(), input.end(), ']'), input.end());
    input.erase(remove(input.begin(), input.end(), '\"'), input.end());
    input.erase(remove(input.begin(), input.end(), ','), input.end());
}

void readTrees(vector<double> & trees0, vector<double> & trees1, vector<double> & trees2, vector<vector<double>> & trees_w ,string nameOfTheFile) {
    ifstream file(nameOfTheFile);
    if (file.is_open())
    {
        string input;
        int totalNodes = 0;

        vector<double> c1;
        vector<double> c2;
        vector<double> c3;
        vector<double> c4;

        while (file >> input)
        {
            cleanInput(input);
            switch(totalNodes % 7) {
                case 0:
                    trees0.push_back(stof(input));
                    break;
                case 1:
                    trees1.push_back(stof(input));
                    break;
                case 2:
                    trees2.push_back(stof(input));
                    break;
                case 3:
                    c1.push_back(stof(input));
                    break;
                case 4:
                    c2.push_back(stof(input));
                    break;
                case 5:
                    c3.push_back(stof(input));
                    break;
                case 6:
                    c4.push_back(stof(input));
                    break;
            }
            totalNodes += 1;
        }
        trees_w.push_back(c1);
        trees_w.push_back(c2);
        trees_w.push_back(c3);
        trees_w.push_back(c4);
    }
}

template <typename T>
void printVectors(vector<vector<T>> & trees, string name) {
    string value;
    for (int i = 0; i < trees.size(); i++) {
        cout << name + ":" << i << endl;
        value = "";
        for (int j = 0; j < trees[i].size(); j++) {
            value += to_string(trees[i][j]) + " ";
        }
        cout << value.substr(0, value.length() - 1) << endl;
    }
}

template <typename T>
void printVector(vector<T> & arr, string name) {
    //cout << name + ":" << i << endl;
    for (int i = 0; i < arr.size(); i++) {
        cout << arr[i] << " ";
    }
    cout << endl;
}

void readRowTestData(vector<vector<complex<double>>> & data0,vector<vector<complex<double>>> & data1, vector<vector<complex<double>>> & data2, string nameOfTheFile) {
    ifstream file(nameOfTheFile);
    if (file.is_open()) {
        string line;
        int readerInt;
        while (getline(file, line)) {
            vector<complex<double>> currentTestData0;
            vector<complex<double>> currentTestData1;
            vector<complex<double>> currentTestData2;

            istringstream iss(line);
            for(std::string::size_type i = 0; i < line.size(); i=i+3) {
                iss >> readerInt;
                currentTestData0.push_back(readerInt);
                iss >> readerInt;
                currentTestData1.push_back(readerInt);
                iss >> readerInt;
                currentTestData2.push_back(readerInt);
            }

            data0.push_back(currentTestData0);
            data1.push_back(currentTestData1);
            data2.push_back(currentTestData2);
        }
    }
}

void readLabels(vector<int> & labes, string nameOfTheFile) {
    ifstream file(nameOfTheFile);
    if (file.is_open()) {

        string line;
        int readerInt;

        getline(file, line);

        istringstream iss(line);

        while (iss >> readerInt) {
            labes.push_back(readerInt);
        }
    }
}

double firstWay(double x0, double x2, double y) {
    return (1 - x2) * (y * (x0 - 1)) + 1;
}

double secondWay(double x0, double x2, double y) {
    return (1 - y - x0) * (x2 * (y - 1) - y) + 1;
}

double thirdWay(double x0, double x2, double y) {
    return (1 - x0) * (x2 * (y - 1) - y) + 1;
}


double calculateOneTree(double z1, double z2, double z3, double c1, double c2, double c3, double c4) {
    double f1 = z1 * z2 * c1;
    double f2 = z1 * (1 - z2) * c2;
    double f3 = (1 - z1) * z3 * c3;
    double f4 = (1 - z1) * (1 - z3) * c4;
    return f1 + f2 + f3 + f4;
}

int findMax(vector<double> & scores, int size) {
    double currentMax = scores[0], currentScore;
    int currentMaxIndex = 0;
    for (int i = 1; i < size; i++) {
        currentScore = scores[i];
        if (currentScore > currentMax) {
            currentMax = currentScore;
            currentMaxIndex = i;
        }
    }
    return currentMaxIndex;
}

int calculateOneSample(vector<vector<double>> & trees, int treeSizePerClass, int numberOfClasses, vector<int> & testDataX0, vector<int> & testDataX2) {
    vector<double> scores;
    double z1, z2, z3, c1, c2, c3, c4, x0, x2, y, classScore, treeScore;
    int currentTreeIndex = 0, currentLabelIndex = 0;
    for (int i = 0; i < numberOfClasses; i++) {
        classScore = 0;
        for (int k = 0; k < treeSizePerClass; k++) {

            x0 = testDataX0[currentLabelIndex], x2 = testDataX2[currentLabelIndex], y = trees[currentTreeIndex][0];
            z1 = thirdWay(x0, x2, y); currentLabelIndex++;

            x0 = testDataX0[currentLabelIndex], x2 = testDataX2[currentLabelIndex], y = trees[currentTreeIndex][1];
            z2 = thirdWay(x0, x2, y); currentLabelIndex++;

            x0 = testDataX0[currentLabelIndex], x2 = testDataX2[currentLabelIndex], y = trees[currentTreeIndex][2];
            z3 = thirdWay(x0, x2, y); currentLabelIndex++;

            c1 = trees[currentTreeIndex][3], c2 = trees[currentTreeIndex][4], c3 = trees[currentTreeIndex][5], c4 = trees[currentTreeIndex][6];

            treeScore = calculateOneTree(z1, z2, z3, c1, c2, c3, c4);
            classScore += treeScore;

            currentTreeIndex++;
        }
        scores.push_back(classScore);
    }

    return findMax(scores, numberOfClasses);
}

double evaluateAllSamples(vector<vector<double>> & trees, int treeSizePerClass, int numberOfClasses, vector<vector<int>> & testDataX0, vector<vector<int>> & testDataX2, vector<int> & labels) {
    int predictedClass, correctGuess = 0;
    int sampleSize = labels.size();

    for (int labelIndex = 0; labelIndex < sampleSize; labelIndex++) {
        predictedClass = calculateOneSample(trees, treeSizePerClass, numberOfClasses, testDataX0[labelIndex], testDataX2[labelIndex]);
        if (predictedClass == (labels[labelIndex])) {
            correctGuess++;
        }
    }
    return double(correctGuess) / double(sampleSize);

}


void compareCipher(CryptoContext<DCRTPoly> &cc,vector<Ciphertext<DCRTPoly>> &z, vector<Ciphertext<DCRTPoly>> &x0, vector<Ciphertext<DCRTPoly>> &x2,
                   vector<double> &y, usint testN) {
    // (1 - x0) * (x2 * (y - 1) - y) + 1
    vector<double> ones(y.size(),1);
    vector<double> y_minus_one(y.size());
    for (usint i = 0; i <y.size(); i++)
        y_minus_one[i] = y[i] - ones[i];

    Plaintext ones_pt = cc->MakeCKKSPackedPlaintext(ones);
    Plaintext y_minus_one_pt = cc->MakeCKKSPackedPlaintext(y_minus_one);
    Plaintext y_pt = cc->MakeCKKSPackedPlaintext(y);
    //    #pragma omp parallel for
    for (usint i = 0; i <testN; i++){
        Ciphertext<DCRTPoly> one_minus_x0 = cc->EvalSub(ones_pt,x0[i]);
        cc->EvalMultMutable(x2[i],y_minus_one_pt);
        cc->Rescale(x2[i]);
        // Check Do Rescale
        cc->EvalSubMutable(x2[i],y_pt);
        cc->EvalMultMutable(one_minus_x0,x2[i]);
        cc->Rescale(one_minus_x0);
        z[i] = cc->EvalAdd(one_minus_x0,ones_pt);
    }

}

void calculateTree(CryptoContext<DCRTPoly> &cc, vector<Ciphertext<DCRTPoly>> &res, vector<Ciphertext<DCRTPoly>> &z1, vector<Ciphertext<DCRTPoly>> &z2,
                   vector<Ciphertext<DCRTPoly>> &z3, vector<vector<double>> &leaves, usint testN) {

    vector<double> ones(testN,1);
    Plaintext ones_pt = cc->MakeCKKSPackedPlaintext(ones);

    Plaintext c1 = cc->MakeCKKSPackedPlaintext(leaves[0]);
    Plaintext c2 = cc->MakeCKKSPackedPlaintext(leaves[1]);
    Plaintext c3 = cc->MakeCKKSPackedPlaintext(leaves[2]);
    Plaintext c4 = cc->MakeCKKSPackedPlaintext(leaves[3]);
//    #pragma omp parallel for
    for (usint i = 0; i < res.size(); ++i) {

        //    f1 = z1 * z2 * c1;
        Ciphertext<DCRTPoly> temp;
        temp = cc->EvalMult(z2[i],c1); // depth 2
        cc->RescaleInPlace(temp); // depth 1
        res[i] = cc->EvalMult(temp,z1[i]); // depth 2

        //    f2 = z1 * (1 - z2) * c2;
        temp = cc->EvalSub(ones_pt,z2[i]); // depth 1
        temp = cc->EvalMult(temp,c2); // depth 2
        cc->RescaleInPlace(temp); // depth 1
        temp = cc->EvalMult(temp,z1[i]); // depth 2
        cc->EvalAddInPlace(res[i],temp); // depth 2

        //    f3 = (1 - z1) * z3 * c3;
        Ciphertext<DCRTPoly> temp2 = cc->EvalSub(ones_pt,z1[i]); // depth 1
        temp = cc->EvalMult(temp2,c3); // depth 2
        cc->RescaleInPlace(temp); // depth 1
        temp = cc->EvalMult(temp,z3[i]); // depth 2
        cc->EvalAddInPlace(res[i],temp); // depth 2

        //    f4 = (1 - z1) * (1 - z3) * c4;
        temp = cc->EvalSub(ones_pt,z3[i]); // depth 1
        temp = cc->EvalMult(temp,c4); // depth 2
        cc->RescaleInPlace(temp); // depth 1
        // temp2 = (1 - z1)
        temp = cc->EvalMult(temp,temp2); // depth 2
        cc->EvalAddInPlace(res[i],temp); // depth 2

        //    return f1 + f2 + f3 + f4;
        cc->RescaleInPlace(res[i]); // depth 1
    }
}

int main() {

    TimeVar t;
    TimeVar tAll;

    TIC(tAll);
//    usint nr_class = 11;
//    double keyGenTime(0.0);
//    double encryptionTime(0.0);
//    double computationTime(0.0);
//    double decryptionTime(0.0);
//    double endToEndTime(0.0);

    vector<double> trees_root;
    vector<double> trees_left;
    vector<double> trees_right;
    vector<vector<double>> trees_leaves;
    readTrees(trees_root,trees_left, trees_right, trees_leaves, "../xgboost-data/encodedTree.txt");

    vector<vector<complex<double>>> X0_root;
    vector<vector<complex<double>>> X0_left;
    vector<vector<complex<double>>> X0_right;
    readRowTestData(X0_root,X0_left, X0_right, "../xgboost-data/testDatax0_ordered_raw.txt");

    vector<vector<complex<double>>> X2_root;
    vector<vector<complex<double>>> X2_left;
    vector<vector<complex<double>>> X2_right;
    readRowTestData(X2_root,X2_left, X2_right, "../xgboost-data/testDatax2_ordered_raw.txt");

    vector<int> labels;
    readLabels(labels, "../xgboost-data/labels.txt");

    size_t testN = labels.size();  // how many sample for predicting
    cout << "Test Data Size: " << testN << endl;
    testN = 10;
    cout << "Tree and Test Data are ready to go...\n";

//// KEY GENERATION
    TIC(t);
    uint32_t multDepth = 3;

    // 0 means the library will choose it based on securityLevel
    uint32_t m = 8192;
//    uint32_t batchSize = m>>2;

    CryptoContext<DCRTPoly> cc =
            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(multDepth,
                                                                 54,
                                                                 0,
                                                                 HEStd_128_classic,
                                                                 m,
                                                                 APPROXAUTO,
                                                                 BV,
                                                                 4,
                                                                 2,
                                                                 54,
                                                                 0,
                                                                 SPARSE);

    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    cc->Enable(LEVELEDSHE);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeysGen(keyPair.secretKey);
    cc->EvalSumKeyGen(keyPair.secretKey);
//    keyGenTime = TOC(t);

    //// ENCODE and ENCRYPTION
    TIC(t);

    vector<Ciphertext<DCRTPoly>> X0_ro(testN);
    vector<Ciphertext<DCRTPoly>> X2_ro(testN);
    vector<Ciphertext<DCRTPoly>> X0_l(testN);
    vector<Ciphertext<DCRTPoly>> X2_l(testN);
    vector<Ciphertext<DCRTPoly>> X0_r(testN);
    vector<Ciphertext<DCRTPoly>> X2_r(testN);

//#pragma omp parallel for
    for (size_t i = 0; i < testN; i++) {
        Plaintext x0Temp = cc->MakeCKKSPackedPlaintext(X0_root[i]);
        Plaintext x2Temp = cc->MakeCKKSPackedPlaintext(X2_root[i]);
        X0_ro[i] = cc->Encrypt(keyPair.publicKey, x0Temp);
        X2_ro[i] = cc->Encrypt(keyPair.publicKey, x2Temp);

        x0Temp = cc->MakeCKKSPackedPlaintext(X0_left[i]);
        x2Temp = cc->MakeCKKSPackedPlaintext(X2_left[i]);
        X0_l[i] = cc->Encrypt(keyPair.publicKey, x0Temp);
        X2_l[i] = cc->Encrypt(keyPair.publicKey, x2Temp);

        x0Temp = cc->MakeCKKSPackedPlaintext(X0_right[i]);
        x2Temp = cc->MakeCKKSPackedPlaintext(X2_right[i]);
        X0_r[i] = cc->Encrypt(keyPair.publicKey, x0Temp);
        X2_r[i] = cc->Encrypt(keyPair.publicKey, x2Temp);
    }
    cout << "X is encrypted." << endl;

    //cout << "Correctness:" << evaluateAllSamples(trees, 128, 11, testDataX0, testDataX2, labels) << endl;

    //// XGBOOST COMPARE
    vector<Ciphertext<DCRTPoly>> z1(testN);
    vector<Ciphertext<DCRTPoly>> z2(testN);
    vector<Ciphertext<DCRTPoly>> z3(testN);

    compareCipher(cc,z1,X0_ro,X2_ro, trees_root,testN);
    compareCipher(cc,z3,X0_r,X2_r, trees_right,testN);
    compareCipher(cc,z2,X0_l,X2_l, trees_left,testN);

    vector<Ciphertext<DCRTPoly>> res(testN);

    calculateTree(cc,res,z1,z2,z3,trees_leaves,testN);

    vector<Ciphertext<DCRTPoly>> sum(testN);

//    for (usint i = 0; i < testN; ++i) {
//        sum[i] = cc->EvalSum(res[i],128);
//    }

    vector<vector<double>> results(testN);
    for (uint32_t i = 0; i < testN; ++i)
        results[i] = vector<double>(128*11);

    for (size_t s = 0; s < testN; ++s) {
        Plaintext result_dec_values;
        cc->Decrypt(keyPair.secretKey, res[s], &result_dec_values);
        result_dec_values->SetLength(128*11);
        for (size_t i = 0; i < 128*11; ++i) {
            auto value = result_dec_values->GetCKKSPackedValue();
            results[s][i] = value[i].real();
        }
    }
    vector<double> temp(11,0);
    for (usint i = 0; i < 11; ++i) {
        for (int j = 0; j < 128; ++j) {
            temp[i] += results[0][i*128+j];
        }
        cout << "Class " << i << " value: "<<temp[i] << endl;
    }

    cout << "Operation is done!! \n";






    return 0;
}
