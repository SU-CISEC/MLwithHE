//
// Created by Ferhat Yaman on 26.09.2020.
//

#include "utils.cpp"
using namespace std;

int main() {

    TimeVar t;
    TimeVar tAll;

    TIC(tAll);

    double keyGenTime(0.0);
    double encryptionTime(0.0);
    double computationTime(0.0);
    double decryptionTime(0.0);
    double endToEndTime(0.0);

    cout << "\n======SVM PALISADE Solution========\n" << std::endl;

    vector<double> yData;
    vector<vector<double>> sData;

    string X_Filename = "./data/XScaled.csv";
    size_t N = 200;
    size_t M = 46703;
    ReadXFile(sData, X_Filename,N,M);
    ReadYDataFile()

    return 0;
}