/*
 *Created by sselcan
*/

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <utils/debug.h>
#include "palisade.h"

using namespace std;
using namespace lbcrypto;

/*
 * Original version of xgboost prediction module
 * data set {-1, 0, 1}
 *
 * */
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
            for(std::string::size_type i = 0; i < 128*11*3; i=i+3) {
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
                   vector<Ciphertext<DCRTPoly>> &z3, vector<vector<double>> &leaves, usint testN) {

    vector<double> ones(128*11,1);
    Plaintext ones_pt = cc->MakeCKKSPackedPlaintext(ones);

    Plaintext c1 = cc->MakeCKKSPackedPlaintext(leaves[0]);
    Plaintext c2 = cc->MakeCKKSPackedPlaintext(leaves[1]);
    Plaintext c3 = cc->MakeCKKSPackedPlaintext(leaves[2]);
    Plaintext c4 = cc->MakeCKKSPackedPlaintext(leaves[3]);
#pragma omp parallel for
    for (usint i = 0; i < res.size(); ++i) {

        //    f1 = z1 * z2 * c1;
        Ciphertext<DCRTPoly> temp;
        temp = cc->EvalMult(z2[i],c1); // depth 2
        temp = cc->Rescale(temp); // depth 1
        res[i] = cc->EvalMult(temp,z1[i]); // depth 2

        //    f2 = z1 * (1 - z2) * c2;
        temp = cc->EvalSub(ones_pt,z2[i]); // depth 1
        temp = cc->EvalMult(temp,c2); // depth 2
        temp = cc->Rescale(temp); // depth 1
        temp = cc->EvalMult(temp,z1[i]); // depth 2
        res[i] = cc->EvalAdd(res[i],temp); // depth 2

        //    f3 = (1 - z1) * z3 * c3;
        Ciphertext<DCRTPoly> temp2 = cc->EvalSub(ones_pt,z1[i]); // depth 1
        temp = cc->EvalMult(temp2,c3); // depth 2
        temp = cc->Rescale(temp); // depth 1
        temp = cc->EvalMult(temp,z3[i]); // depth 2
        res[i] = cc->EvalAdd(res[i],temp); // depth 2

        //    f4 = (1 - z1) * (1 - z3) * c4;
        temp = cc->EvalSub(ones_pt,z3[i]); // depth 1
        temp = cc->EvalMult(temp,c4); // depth 2
        temp = cc->Rescale(temp); // depth 1
        // temp2 = (1 - z1)
        temp = cc->EvalMult(temp,temp2); // depth 2
        res[i] = cc->EvalAdd(res[i],temp); // depth 2

        //    return f1 + f2 + f3 + f4;
        res[i] = cc->Rescale(res[i]); // depth 1
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

    //Read the encoded model
    vector<double> trees_root;
    vector<double> trees_left;
    vector<double> trees_right;
    vector<vector<double>> trees_leaves;
    readTrees(trees_root,trees_left, trees_right, trees_leaves, "../Data_XGB/encodedTree.txt");

    //Get the encoded test data
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
    TIC(t);
    CryptoContext<DCRTPoly> cc =
            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(4, 50, 1408, HEStd_128_classic);
    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);
    cc->Enable(LEVELEDSHE);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeysGen(keyPair.secretKey);
    cc->EvalSumKeyGen(keyPair.secretKey);
    vector<int> ilist = {1,2,4,8,16,32,64};
    cc->EvalAtIndexKeyGen(keyPair.secretKey, ilist);

    keyGenTime = TOC(t);
    //// ENCODE and ENCRYPTION
    TIC(t);
    vector<Ciphertext<DCRTPoly>> X0_ro(testN);
    vector<Ciphertext<DCRTPoly>> X2_ro(testN);
    vector<Ciphertext<DCRTPoly>> X0_l(testN);
    vector<Ciphertext<DCRTPoly>> X2_l(testN);
    vector<Ciphertext<DCRTPoly>> X0_r(testN);
    vector<Ciphertext<DCRTPoly>> X2_r(testN);

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
    //// XGBOOST COMPARE
    TIC(t);
    vector<Ciphertext<DCRTPoly>> z1(testN);
    vector<Ciphertext<DCRTPoly>> z2(testN);
    vector<Ciphertext<DCRTPoly>> z3(testN);


    compareCipher(cc, z1, X0_ro, X2_ro, trees_root, testN);
    compareCipher(cc, z3, X0_r, X2_r, trees_right, testN);
    compareCipher(cc, z2, X0_l, X2_l, trees_left, testN);
    compCipher = TOC(t);


    compCipher = TOC(t);
    TIC(t);

    vector<Ciphertext<DCRTPoly>> res(testN);

    calculateTree(cc,res,z1,z2,z3,trees_leaves,testN);
    calcTree = TOC(t);


    TIC(t);
    vector<Ciphertext<DCRTPoly>> sum(testN);
#pragma omp parallel for
    for (usint i = 0; i < testN; ++i) {
        for (usint j = 7; j > 0; j--){
            auto temp = cc->EvalAtIndex(res[i], pow(2,j-1));
            res[i] = cc->EvalAdd(res[i],temp);
        }
        sum[i] = res[i];
    }


    ////DECRYPTION
    vector<vector<double>> results(testN,vector<double>(11));
#pragma omp parallel for
    for (size_t s = 0; s < testN; ++s) {
        Plaintext result_dec_values;
        cc->Decrypt(keyPair.secretKey, res[s], &result_dec_values);
        result_dec_values->SetLength(128*11);
        auto value = result_dec_values->GetCKKSPackedValue();
        for (size_t i = 0; i < 11; ++i) {
            results[s][i] = value[i*128].real();
        }
    }
    decryptionTime = TOC(t);

    cout << "\nKey Generation Time: \t\t" << keyGenTime/1000 << " s" << endl;
    cout << "Encoding and Encryption Time: \t" << encryptionTime/1000 << " s" << endl;
    cout << "Compare Time: \t\t" << compCipher/1000 << " s" << endl;
    cout << "Calculate Tree Time Time: \t\t" << calcTree/1000 << " s" << endl;
    cout << "Decryption & Decoding Time: \t" << decryptionTime/1000 << " s" << endl;

    endToEndTime = TOC(tAll);
    cout << "\nEnd-to-end Runtime: \t\t" << endToEndTime/1000 << " s" << endl;


    cout << results[0].size() << endl;
    vector<double> sss(11,0);
    for (int i = 0; i < 11; ++i) {
        cout << "Class " << i << " value: "<<results[0][i] << endl;
    }

    cout << "Operation is done!! \n";
    checkAcc(results, labels);

    return 0;
}