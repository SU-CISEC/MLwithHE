SVM Classification with HE 
=====================================
This project is created to predict the type of tumor given by mutations on genes location. This project is started for a competition on [iDASH 2020](http://www.humangenomeprivacy.org/2020/competition-tasks.html).
The is implementing of classifying predictions on encrypted data which force us to use Homomorphic Encryption(HE). 
The project is written by great HE library PALISADE.

This github repository includes the implementation of prediction function of Support Vector Machine(SVM) algorithm.
Those algorithms are inspired from [Scikit-Learn Library](https://scikit-learn.org/stable/index.html).

The repo includes the following files:
* svm-HE.cpp - research prototype for the SVM algorithm on PALISADE.
* data/data.test.txt - a generated test data set by script.sh, consists of 141 features.
* model/coef.csv - pre-trained model parameters for every class()
* model/rho.csv  - pre-trained model parameters for every class()
* script.sh - preprocessing script which is written in bash

**Test Data, Test Executables and Metric Scripts**

* data/SampleData - Contains sample data and Python Metric Scripts
* svm-HE-test.cpp - Runs on test Sample data and prepares different results other than data
* svm.cpp - Runs on sample data without HE as a debug tool for Homomorphic Encryption

How to Build and Run the Prototypes
=====================================

### Preprocessing

1. This script can be run with bash version 3 or higher. To learn your bash version you can run:
    ```
     $ bash --version
    ```
2. The first line of the script.sh file should be ```#!/usr/local/bin/bash``` or ```#!/bin/bash``` with respect to your bash PATH. The default path is ```#!/bin/bash```.

3. Set execute permission on your script:
    ```
     $ chmod +x script.sh
    ```
4. To run your script, enter:
    ```
     $ ./script.sh
    ```
5. Enter the path of name of the test files without _CNs.txt and _variants.txt extension. If names of the test files are: test_CNs.txt and test_variants.txt it should be
    ```
     $ /file_path/test
    ```


### Homomorphic Prediction

1. Install PALISADE v1.9.1 from [PALISADE Development Repository](https://gitlab.com/palisade/palisade-development/-/tree/release-v1.9.1). Follow the instructions provided in [README.md](https://gitlab.com/palisade/palisade-development/-/blob/release-v1.9.1/README.md).
   
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
    ./svm-HE 
    ```

5. The results will be written to the "data/results.csv" file.

 # Docker

1- To run container save this instruction as Dockerfile and follow the instructions below

```
FROM ferhatyaman/ml-with-he:latest

WORKDIR /home/user/MLwithHE/cmake-build-debug-docker

#Make sure your input data is inside the Challenge folder or change the name of folder
COPY ./Challenge ../data

CMD ["/bin/bash", "-c", "../data/script.sh ; ./svm-HE"]
```
2- Build dockerfile using command below
```
sudo docker build --tag svm-predictor .
```
3- Run Container from latest produced image and give file name as input accordingly
```
sudo docker run -it --name idash20-sabanci svm-predictor
#Give input as ../data/FILENAME Ex: ../data/Bladder_challenge  for Bladder_challenge_variants.txt
```
4- Copy results from container to localhost
```
sudo docker cp idash20-sabanci:/home/user/MLwithHE/data/results.csv .
```
XGBoost Classification with HE 
=====================================
This XGBoost Classification is a continuation of the project we started with SVM Classification with HE (see above).
This github repository includes the implementation of prediction function of XGBOOST algorithm. The algorithm is implemented using both PALISADE and SEAL libraries.

Following files added with XGBoost Support:
* xgboost.cpp - research prototype for the XGBoost algorithm on PALISADE.
* xgboost_v3.cpp - similar to version0, minor improvements on tree score calculation
* xgboost_encModel.cpp -  research prototype for the XGBoost algorithm on PALISADE. Both model and test data is encrypted
* xgboost_BFV.cpp - research prototype for the XGBoost algorithm on SEAL, BFV scheme is used
* Data_XGB - contains sample test data and encoded XGBoost model

For more detail see: [ML with HE: Privacy Preserving Machine Learning Inferences for Genome Studies](https://arxiv.org/abs/2110.11446)

### How to Build and Run the Prototypes

1. Install PALISADE v1.11.3 from [PALISADE Development Repository](https://gitlab.com/palisade/palisade-development/-/tree/release-v1.11.3). Follow the instructions provided in [README.md](https://gitlab.com/palisade/palisade-development/-/blob/release-v1.11.3/README.md).

2. Install SEAL v3.6.6 from [Microsoft SEAL Repository](https://github.com/microsoft/SEAL/releases/tag/v3.6.6). Follow the [instructions](https://github.com/microsoft/SEAL#building-microsoft-seal-manually)

3. Clone this repository to a local directory and switch to this directory.

5. Create a directory where the binaries will be built. The typical choice is a subfolder "build". In this case, run the following commands:

    ```
    mkdir build
    cd build
    cmake ..
    make
    ```

4. Run the following command to execute the SVM prototype:
    ```
    ./xgboost_bfv 
    ```

# Contributors
This project is supported by Sabanci University. Main contributors are below.
* Berke Dilekoğlu
* Ceren Yıldırım
* Ferhat Yaman
* Şeyma Mağara
* Dr. Öznur Taştan
* Dr. Kamer Kaya
* Dr. Erkay Savaş
