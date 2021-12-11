#include <iostream>
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
        for (uint i = 0; i < (c1).size(); i++){
            c1[i] = c1[i] - c2[i];
            c2[i] = c2[i] - c4[i];
            c3[i] = c4[i] - c3[i];
        }
        trees_w.push_back(c1); // l1-l2
        trees_w.push_back(c2); // l2-l4
        trees_w.push_back(c3); // l4-l3
        trees_w.push_back(c4); // l4
    }
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


void checkAcc(vector<vector<double>> results ,vector<int> labels){
    int counter = 0;
    for(uint i = 0; i < results.size(); i++){
        if (labels[i] == distance(results[i].begin(),max_element(results[i].begin(), results[i].end()))){
            counter++;
        }
    }
    cout << counter << " " << results.size() << endl;

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
    #pragma omp parallel for
    for (usint i = 0; i <testN; i++){
        auto one_minus_x0 = cc->EvalSub(ones_pt,x0[i]);
        x2[i] = cc->EvalMult(x2[i],y_minus_one_pt);
        x2[i] = cc->Rescale(x2[i]);
        // Check Do Rescale
        x2[i]= cc->EvalSub(x2[i],y_pt);
        one_minus_x0 = cc->EvalMult(one_minus_x0,x2[i]);
        one_minus_x0= cc->Rescale(one_minus_x0);
        z[i] = cc->EvalAdd(one_minus_x0,ones_pt);
    }

}

void calculateTree(CryptoContext<DCRTPoly> &cc, vector<Ciphertext<DCRTPoly>> &res, vector<Ciphertext<DCRTPoly>> &z1, vector<Ciphertext<DCRTPoly>> &z2,
                   vector<Ciphertext<DCRTPoly>> &z3, Ciphertext<DCRTPoly> &leaves1, Ciphertext<DCRTPoly> &leaves2,Ciphertext<DCRTPoly> &leaves3,Ciphertext<DCRTPoly> &leaves4, usint testN) {

    vector<double> ones(128*11,1);
    Plaintext ones_pt = cc->MakeCKKSPackedPlaintext(ones);

#pragma omp parallel for
    for (usint i = 0; i < res.size(); ++i) {

        //    z2 * (l1-l2)+ (l2-l4);
        Ciphertext<DCRTPoly> temp;
        temp = cc->EvalMult(z2[i],leaves1); // depth 2
        temp = cc->Rescale(temp); // depth 1
        temp = cc->EvalAdd(temp,leaves2); // depth 1

        //    z1 * z2 * (l1-l2)+ (l2-l4);
        res[i] = cc->EvalMult(z1[i],temp); // depth 2

        //    (z1 - 1) * z3 * (l4-l3);
        Ciphertext<DCRTPoly> temp2 = cc->EvalSub(z1[i],ones_pt); // depth 1
        temp = cc->EvalMult(temp2,leaves3); // depth 2
        temp = cc->Rescale(temp); // depth 1
        temp = cc->EvalMult(temp,z3[i]); // depth 2
        res[i] = cc->EvalAdd(res[i],temp); // depth 2

        //    return f;
        res[i] = cc->Rescale(res[i]); // depth 1
        res[i] = cc->EvalAdd(res[i], leaves4); // depth 1

    }
}


int main() {
    TimeVar t;
    TimeVar tAll;

    TIC(tAll);
    double keyGenTime(0.0);
    double encryptionTime(0.0);
    double compCipher(0.0);
    double calcTree(0.0);
    double decryptionTime(0.0);
    double endToEndTime(0.0);

    vector<double> trees_root;
    vector<double> trees_left;
    vector<double> trees_right;
    vector<vector<double>> trees_leaves;
    readTrees(trees_root,trees_left, trees_right, trees_leaves, "../Data_XGB/encodedTree.txt");

    vector<vector<complex<double>>> X0_root;
    vector<vector<complex<double>>> X0_left;
    vector<vector<complex<double>>> X0_right;
    readRowTestData(X0_root,X0_left, X0_right, "../Data_XGB/demoTestData_x0.txt");

    vector<vector<complex<double>>> X2_root;
    vector<vector<complex<double>>> X2_left;
    vector<vector<complex<double>>> X2_right;
    readRowTestData(X2_root,X2_left, X2_right, "../Data_XGB/demoTestData_x2.txt");

    vector<int> labels;
    readLabels(labels, "../Data_XGB/labels.txt");

    size_t testN = labels.size();  // how many sample for predicting
    cout << "Test Data Size: " << testN << endl;
    cout << "Tree and Test Data are ready to go...\n";

    //// KEY GENERATION
    CryptoContext<DCRTPoly> cc = CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
                    5, 50, 1408, HEStd_128_classic);
    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    cc->Enable(LEVELEDSHE);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeysGen(keyPair.secretKey);
    keyGenTime = TOC(t);

    //// ENCODE and ENCRYPTION
    TIC(t);
    vector<Ciphertext<DCRTPoly>> X0_ro(testN);
    vector<Ciphertext<DCRTPoly>> X2_ro(testN);
    vector<Ciphertext<DCRTPoly>> X0_l(testN);
    vector<Ciphertext<DCRTPoly>> X2_l(testN);
    vector<Ciphertext<DCRTPoly>> X0_r(testN);
    vector<Ciphertext<DCRTPoly>> X2_r(testN);

    vector<Ciphertext<DCRTPoly>> leaves(4);

    #pragma omp parallel for
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
    encryptionTime = TOC(t);
    Plaintext temp = cc->MakeCKKSPackedPlaintext(trees_leaves[0]);
    auto leaves1 = cc->Encrypt(keyPair.publicKey, temp);
    temp = cc->MakeCKKSPackedPlaintext(trees_leaves[1]);
    auto leaves2 = cc->Encrypt(keyPair.publicKey, temp);
    temp = cc->MakeCKKSPackedPlaintext(trees_leaves[2]);
    auto leaves3 = cc->Encrypt(keyPair.publicKey, temp);
    temp = cc->MakeCKKSPackedPlaintext(trees_leaves[3]);
    auto leaves4 = cc->Encrypt(keyPair.publicKey, temp);
    //// XGBOOST COMPARE
    TIC(t);
    vector<Ciphertext<DCRTPoly>> z1(testN);
    vector<Ciphertext<DCRTPoly>> z2(testN);
    vector<Ciphertext<DCRTPoly>> z3(testN);


    compareCipher(cc, z1, X0_ro, X2_ro, trees_root, testN);
    compareCipher(cc, z3, X0_r, X2_r, trees_right, testN);
    compareCipher(cc, z2, X0_l, X2_l, trees_left, testN);
    compCipher = TOC(t);
    vector<Ciphertext<DCRTPoly>> res(testN);

    TIC(t);
    calculateTree(cc,res,z1,z2,z3,leaves1,leaves2,leaves3,leaves4,testN);

    calcTree = TOC(t);
    ////DECRYPTION
    TIC(t);
    vector<vector<double>> results(testN, vector<double>(128*11));

#pragma omp parallel for
    for (size_t s = 0; s < testN; ++s) {
        Plaintext result_dec_values;
        cc->Decrypt(keyPair.secretKey, res[s], &result_dec_values);
        result_dec_values->SetLength(128*11);
        for (size_t i = 0; i < 128*11; ++i) {
            auto value = result_dec_values->GetCKKSPackedValue();
            results[s][i] = value[i].real();
        }
    }
    decryptionTime = TOC(t);

    cout << "\nKey Generation Time: \t\t" << keyGenTime/1000 << " s" << endl;
    cout << "Encoding and Encryption Time: \t" << encryptionTime/1000 << " s" << endl;
    cout << "Compare Time: \t\t" << compCipher/1000 << " s" << endl;
    cout << "Calculation Time: \t\t" << calcTree/1000 << " s" << endl;
    cout << "Decryption & Decoding Time: \t" << decryptionTime/1000 << " s" << endl;

    endToEndTime = TOC(tAll);
    cout << "\nEnd-to-end Runtime: \t\t" << endToEndTime/1000 << " s" << endl;


    vector<vector<double>> sss(testN, vector<double>(11,0));
    for (usint n = 0; n < testN; ++n) {
        for (usint i = 0; i < 11; ++i) {
            for (int j = 0; j < 128; ++j) {
                sss[n][i] += results[n][i * 128 + j];
            }}
    }
    cout << "Operation is done!! \n";
    checkAcc(sss, labels);

    return 0;
}