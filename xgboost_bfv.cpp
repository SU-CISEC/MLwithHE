/*
Created by cerenyildirim
*/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include "seal/seal.h"
#include <ctime>
#include <climits>
#include <omp.h>
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace seal;

void readTrees(SEALContext &cc, vector<int64_t> & trees0, vector<int64_t> & trees1, vector<int64_t> & trees2, vector<vector<int64_t>> & trees_w ,string nameOfTheFile) {
    ifstream file(nameOfTheFile);
    BatchEncoder batch_encoder(cc);
    if (file.is_open())
    {
        int input;
        int totalNodes = 0;
        vector<int64_t> c1;
        vector<int64_t> c2;
        vector<int64_t> c3;
        vector<int64_t> c4;

        while (file >> input)
        {
            switch(totalNodes % 7) {
                case 0:
                    trees0[totalNodes / 7] = input;
                    break;
                case 1:
                    trees1[totalNodes / 7] = input;
                    break;
                case 2:
                    trees2[totalNodes / 7] = input;
                    break;
                case 3:
                    c1.push_back((input));
                    break;
                case 4:
                    c2.push_back(input);
                    break;
                case 5:
                    c3.push_back(input);
                    break;
                case 6:
                    c4.push_back(input);
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



void readRowTestData(SEALContext & context, vector<vector<int64_t>> & data0,vector<vector<int64_t>> & data1, vector<vector<int64_t>> & data2, string nameOfTheFile) {
    ifstream file(nameOfTheFile);
    BatchEncoder batch_encoder(context);
    if (file.is_open()) {
        string line;
        int64_t readerInt;

        while (getline(file, line)) {
            int index = 0;
            size_t slot_count = batch_encoder.slot_count();
            vector<int64_t> currentTestData0(slot_count, 0ULL);
            vector<int64_t> currentTestData1(slot_count, 0ULL);
            vector<int64_t> currentTestData2(slot_count, 0ULL);

            istringstream iss(line);
            for(std::string::size_type i = 0; i < line.size(); i=i+3) {
                iss >> readerInt;
                currentTestData0[index] = (readerInt);
                iss >> readerInt;
                currentTestData1[index] = (readerInt);
                iss >> readerInt;
                currentTestData2[index] = (readerInt);
                index++;
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

void compareCipher(KeyGenerator &keygen, Encryptor &encryptor, Evaluator &evaluator, BatchEncoder &batch_encoder, RelinKeys relin_keys, PublicKey public_key,vector<Ciphertext> &z, vector<Ciphertext> &x0, vector<Ciphertext> &x2, vector<int64_t> &y, int testN) {
// (1 - x0) * (x2 * (y - 1) - y) + 1

    size_t slot_count = batch_encoder.slot_count();
    vector<int64_t> ones(slot_count, 1ULL);
    vector<int64_t> y_minus_one(slot_count, 0ULL);

    for (uint i = 0; i <y.size(); i++)
        y_minus_one[i] = y[i] - ones[i];

    Plaintext ones_pt; batch_encoder.encode(ones, ones_pt);
    Ciphertext ones_ct; encryptor.encrypt(ones_pt, ones_ct);
    Plaintext y_minus_one_pt; batch_encoder.encode(y_minus_one, y_minus_one_pt);
    Plaintext y_pt; batch_encoder.encode(y, y_pt);
#pragma omp parallel for
    for (int i = 0; i <testN; i++){
        Ciphertext one_minus_x0; evaluator.sub(ones_ct, x0[i], one_minus_x0);
        //cout << "    + Noise budget in result: " << decryptor.invariant_noise_budget(x0[i]) << " bits" << endl;
        evaluator.multiply_plain_inplace(x2[i], y_minus_one_pt);
        //cout << "    + Noise budget in result: " << decryptor.invariant_noise_budget(x2[i]) << " bits" << endl;
        evaluator.relinearize_inplace(x2[i], relin_keys);
        // Check Do Rescale
        evaluator.sub_plain_inplace(x2[i],y_pt);
        evaluator.multiply_inplace(one_minus_x0,x2[i]);
        evaluator.relinearize_inplace(one_minus_x0, relin_keys);
        evaluator.add_plain(one_minus_x0, ones_pt, z[i]);
    }
}

void calculateTree(KeyGenerator &keygen, Encryptor &encryptor, Evaluator &evaluator, BatchEncoder &encoder, RelinKeys relin_keys, PublicKey public_key, vector<Ciphertext> &res, vector<Ciphertext> &z1, vector<Ciphertext> &z2,
                   vector<Ciphertext> &z3, vector<vector<int64_t>> &leaves, int testN) {

    size_t slot_count = encoder.slot_count();
    vector<int64_t> ones(slot_count, 1ULL);
    Plaintext ones_pt; encoder.encode(ones, ones_pt);
    Ciphertext ones_ct; encryptor.encrypt(ones_pt, ones_ct);

    Plaintext c1; encoder.encode(leaves[0], c1);
    Plaintext c2; encoder.encode(leaves[1], c2);
    Plaintext c3; encoder.encode(leaves[2], c3);
    Plaintext c4; encoder.encode(leaves[3], c4);
#pragma omp parallel for
    for (uint i = 0; i < res.size(); ++i) {

        //    z2 * (l1-l2)+ (l2-l4);
        Ciphertext temp;
        evaluator.multiply_plain(z2[i],c1, temp); // depth 2
        evaluator.relinearize_inplace(temp, relin_keys);
        evaluator.add_plain_inplace(temp,c2); // depth 1

        //    z1 * z2 * (l1-l2)+ (l2-l4);
        evaluator.multiply(temp,z1[i], res[i]); // depth 2

        //    (z1 - 1) * z3 * (l4-l3);
        Ciphertext temp2; evaluator.sub(z1[i], ones_ct, temp2); // depth 1
        evaluator.multiply_plain_inplace(temp2,c3); // depth 2
        evaluator.relinearize_inplace(temp2, relin_keys);
        evaluator.multiply_inplace(temp2, z3[i]); // depth 2
        evaluator.add_inplace(res[i], temp2); // depth 2


        //    return f1 + f2 + f3 + f4;
        evaluator.relinearize_inplace(res[i], relin_keys); // depth 1
        evaluator.add_plain_inplace(res[i], c4); // depth 1
    }
}

int main() {


    auto start = high_resolution_clock::now();

    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 16384 / 2;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));


    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    BatchEncoder batch_encoder(context);
    size_t slot_count = batch_encoder.slot_count();

    GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);

    auto keygent = high_resolution_clock::now();

    vector<int64_t> trees_root(slot_count, 0ULL);
    vector<int64_t> trees_left(slot_count, 0ULL);
    vector<int64_t> trees_right(slot_count, 0ULL);
    vector<vector<int64_t>> trees_leaves;
    readTrees(context, trees_root,trees_left, trees_right, trees_leaves, "../Data_XGB/encodedTree_BFV.txt");

    vector<vector<int64_t>> X0_root;
    vector<vector<int64_t>> X0_left;
    vector<vector<int64_t>> X0_right;
    readRowTestData(context, X0_root,X0_left, X0_right, "../Data_XGB/demoTestData_x0.txt");

    vector<vector<int64_t>> X2_root;
    vector<vector<int64_t>> X2_left;
    vector<vector<int64_t>> X2_right;
    readRowTestData(context, X2_root,X2_left, X2_right, "../Data_XGB/demoTestData_x2.txt");

    vector<int> labels;
    readLabels(labels, "../Data_XGB/labels.txt");

    size_t testN = labels.size();  // how many sample for predicting
    //testN = 1;
    cout << "Test Data Size: " << testN << endl;
    cout << "Tree and Test Data are ready to go...\n";
    auto readingt = high_resolution_clock::now();

    //// ENCODE and ENCRYPTION

    vector<Ciphertext> X0_ro(testN);
    vector<Ciphertext> X2_ro(testN);
    vector<Ciphertext> X0_l(testN);
    vector<Ciphertext> X2_l(testN);
    vector<Ciphertext> X0_r(testN);
    vector<Ciphertext> X2_r(testN);



#pragma omp parallel for
    for (size_t i = 0; i < testN; i++) {
        Plaintext x0Temp, x2Temp;
        batch_encoder.encode(X0_root[i], x0Temp);
        batch_encoder.encode(X2_root[i], x2Temp);
        encryptor.encrypt(x0Temp, X0_ro[i]);
        encryptor.encrypt(x2Temp, X2_ro[i]);

        batch_encoder.encode(X0_left[i], x0Temp);
        batch_encoder.encode(X2_left[i], x2Temp);
        encryptor.encrypt(x0Temp, X0_l[i]);
        encryptor.encrypt(x2Temp, X2_l[i]);

        batch_encoder.encode(X0_right[i], x0Temp);
        batch_encoder.encode(X2_right[i], x2Temp);
        encryptor.encrypt(x0Temp, X0_r[i]);
        encryptor.encrypt(x2Temp, X2_r[i]);
    }
    cout << "X is encrypted." << endl;
    //cout << "    + Noise budget in result: " << decryptor.invariant_noise_budget(X0_ro[0]) << " bits" << endl;
    auto enct = high_resolution_clock::now();

    //// XGBOOST COMPARE
    vector<Ciphertext> z1(testN);
    vector<Ciphertext> z2(testN);
    vector<Ciphertext> z3(testN);

    compareCipher(keygen, encryptor, evaluator, batch_encoder, relin_keys, public_key,z1,X0_ro,X2_ro, trees_root,testN);
    compareCipher(keygen, encryptor, evaluator, batch_encoder, relin_keys, public_key,z3,X0_r,X2_r, trees_right,testN);
    compareCipher(keygen, encryptor, evaluator, batch_encoder, relin_keys, public_key,z2,X0_l,X2_l, trees_left,testN);
    //cout << "    + Noise budget in result: " << decryptor.invariant_noise_budget(z1[0]) << " bits" << endl;

    auto compt = high_resolution_clock::now();

    vector<Ciphertext> res(testN);

    calculateTree(keygen, encryptor, evaluator, batch_encoder, relin_keys, public_key,res,z1,z2,z3,trees_leaves, testN);
    cout << "calc tree 1" << endl;
    //cout << "    + Noise budget in result: " << decryptor.invariant_noise_budget(res[0]) << " bits" << endl;

    auto calct = high_resolution_clock::now();
    vector<Ciphertext> sum(testN);
#pragma omp parallel for
    for (uint i = 0; i < testN; ++i) {
        vector<Ciphertext> temp_vec(11);
        Ciphertext temp;
        for (int k = 0; k < 7; k++) {
            evaluator.rotate_rows(res[i], pow(2, k), galois_keys, temp);
            evaluator.add_inplace(res[i], temp);
        }
        sum[i] = res[i];
    }

    vector<vector<int64_t>> results(testN, vector<int64_t>(11));
#pragma omp parallel for
    for (size_t s = 0; s < testN; ++s) {
        Plaintext result_dec_values;
        decryptor.decrypt(sum[s], result_dec_values);
        vector<int64_t> result;
        batch_encoder.decode(result_dec_values, result);
        results[s] = result;
    }
    auto dect = high_resolution_clock::now();

//Quick acc check
//    vector<double> sss(11,0);
//    int counter_all = 0;
//    int corrects = 0;
//    //#pragma omp parallel for
//    for (uint j = 0; j < testN; ++j) {
//        int index = 0;
//        int64_t max = INT_MIN;
//        for (uint i = 0; i < 11; ++i) {
//            if (max < results[j][i*128]) {
//                max = results[j][i*128];
//                index = i;
//            }
//        }
//        //cout << labels[j] << " " << index << endl;
//        if (labels[j] == index){
//            corrects++;
//        }
//        counter_all++;
//    }


    cout << "Operation is done!! \n";
    cout << "KeyGen time " << float(duration_cast<milliseconds>(keygent-start).count())/1000 << " seconds." << endl;
    cout << "reading time " << float(duration_cast<milliseconds>(readingt- keygent).count())/1000 << " seconds." << endl;
    cout << "Encode-Encrypt time " << float(duration_cast<milliseconds>(enct- readingt).count())/1000 << " seconds." << endl;
    cout << "Compare time " << float(duration_cast<milliseconds>(compt- enct).count())/1000 << " seconds." << endl;
    cout << "Calculate Tree time " << float(duration_cast<milliseconds>(calct- compt).count())/1000 << " seconds." << endl;
    cout << "Decyption time " << float(duration_cast<milliseconds>(dect- calct).count())/1000 << " seconds." << endl;
    cout << "End-to-End time " << float(duration_cast<milliseconds>(dect - start).count())/1000 - float(duration_cast<milliseconds>(readingt- keygent).count())/1000<< " seconds." << endl;
    cout << corrects << " " << counter_all << endl;
    return 0;
}
