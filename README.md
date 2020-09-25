SVM and XGBoost Classification with HE 
=====================================
This project is created to predict the type of tumor given by mutations on genes location. This project started for a competition on [iDASH 2020](http://www.humangenomeprivacy.org/2020/competition-tasks.html).
The is implementing of classifying predictions on encrypted data which force us to use Homomorphic Encryption(HE). 
The project is written by 2 great HE library, SEAL and PALISADE.

This github repository includes the implementation of prediction function of Support Vector Machine(SVM) and XGBoost algorithms.
Those algorithms are inspired from [Scikit-Learn Library](https://scikit-learn.org/stable/index.html).

The repo includes the following files:
* svm-seal.cpp - research prototype for the SVM protocol on SEAL.
* svm-palisade.cpp - research prototype for the SVM protocol on PALISADE.
* xgboost-seal.cpp - research prototype for the XGBoost protocol on SEAL.
* xgboost-palisade.cpp - research prototype for the XGBoost protocol on PALISADE.
* data/random_sample.csv - an artificial random data set including 46327 features, 200 location (provided solely for demonstration purposes). 
* data/model_parameters.csv - pre-trained model parameters for every class()
 
How to Build and Run the Prototypes
=====================================

1. Install SEAL v3.5.8 from [Microsoft SEAL](https://github.com/microsoft/SEAL). Follow the instructions provided in [README.md](https://github.com/microsoft/SEAL/blob/master/README.md).
   Install PALISADE v1.9.1 from [PALISADE Development Repository](https://gitlab.com/palisade/palisade-development/-/tree/release-v1.9.1). Follow the instructions provided in [README.md](https://gitlab.com/palisade/palisade-development/-/blob/release-v1.9.1/README.md).
   

2. Clone this repository to a local directory and switch to this directory. 

3. Create a directory where the binaries will be built. The typical choice is a subfolder "build". In this case, run the following commands:
```
mkdir build
cd build
cmake ..
make
```
4. Run the following command to execute the SVM prototype:
```
./svm-seal
./svm-palisade 
```

or

Run the following command to execute the XGBoost prototype:
```
./xgboost-seal
./xgboost-palisade
```

5. The results will be written to the "data" folder. The following output files will be created for both prototypes:


Additional files with outputs of protocol-specific statistics will also be created.

Things To Do
=====================================

1. Parametric file reading for data and encryption of those 

2. Implementation of predict algorithm on plaintext data

3. Implement functions  with HE 