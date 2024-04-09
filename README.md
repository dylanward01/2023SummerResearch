### This application serves to provide a way to find and analyze frequency bands in Univariate, Multivariate, and Functional data, in a graphic-based, user-friendly manner.

## Repository Contents
- The Dockerfiles folder is given just for transparancy's sake. It includes 2 dockerfiles- one utilized to build this repository as a docker image on a PC, and one utilized to build this repository as a docker image on a Mac. 
- Within the folder FrequencyBandAnalysis one can see:
  -   The file app.R, which contains the code for actually running the shiny application.
  -   The 3 other R files, and the 2 cpp files- all of which contain the actual functions and codes for the algorithms that are run.
  -   A folder entitled "Sample Data" which contains simulated data for each of the 3 types of data that this application works on. This folder is not required to be present for the application to run- it is simply there to give the user some sample data to explore the application, and the outputs that it gives.
  -   A folder that contains Vignettes, to better various facets of this application. 

## Running the Application 

Currently, this application can be run either through R, or through Docker. An explanation on how to run this application via any of those means is explained in this [vignette](https://github.com/dylanward01/FrequencyBandAnalysis/blob/main/FrequencyBandAnalysis/Vignettes/Running-The-Shiny-App.pdf), that is seen in this repository

## Within the Application

To see a detailed explanation of how to utilize the application once it's opened, and how to interpret its results, one can consult the [vignette](https://github.com/dylanward01/FrequencyBandAnalysis/blob/main/FrequencyBandAnalysis/Vignettes/Shiny-App-Usage.pdf) seen in this repository.

## Citations

Throughout this application, 3 different papers are mentioned and referred to. Their doi's, and the algorithm produced in each of those papers is mentioned below:
- (Univariate): doi.org/10.1080/01621459.2019.1671199
- (Multivariate): doi.org/10.48550/arXiv.2301.03664
- (Functional): doi.org/10.48550/arXiv.2102.01784
