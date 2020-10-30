# To load image use command below
    #sudo docker load --input ml-with-he.tar
FROM ml-with-he:latest

WORKDIR /home/user/MLwithHE/cmake-build-debug-docker
#Make sure your input data is inside the challenge folder or change the name of folder
COPY ./challenge ../data

CMD ["/bin/bash", "-c", "../data/script.sh ; ./svm-HE"]

# Build dockerfile using command below
    #sudo docker build --tag svm-predictor .

# Run Container from latest produced image and give file name as input accordingly
    #sudo docker run -it --name idash20-sabanci svm-predictor
    #Give input as ../data/FILENAME Ex: ../data/Bladder_challenge  for Bladder_challenge_variants.txt

# Copy results from container to localhost
#sudo docker cp idash20-sabanci:/home/user/MLwithHE/data/results.csv .