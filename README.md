## Repository Contents
- The Dockerfiles folder is given just for transparancy's sake. It includes 2 dockerfiles- one utilized to build this repository as a docker image on a PC, and one utilized to build this repository as a docker image on a Mac. 
- Within the folder FrequencyBandAnalysis one can see:
  -   The file app.R, which contains the code for actually running the shiny application.
  -   The 3 other R files (all with functions in their names), and the 2 cpp files- all of which contain the actual functions and codes for the algorithms that are run.
  -   A folder entitled "Sample Data" which contains simulated data for each of the 3 types of data that this application works on. This folder is not required to be present for the application to run- it is simply there to give the user some sample data to explore the application, and the outputs that it gives.
  -   A folder that contains Vignettes, to better various facets of this application. 

## Running the Application in R

1. Download the FrequencyBandAnalysis folder seen above. 
2. Ensure that the packages listed in the app.R file are installed on your machine.
3. Proceed to run the app.R file, and wait for the application to open in a new window from your R interface.

## Running the Application through Docker

If one wants to use this application, but doesn't want to utilize R, they can run this application through Docker. Namely, Dockerhub serves to provide a repository where updated images for this application are located, as an image for this application is currently located at https://hub.docker.com/r/dylanward01/frequencybandanalysis. This method requires nothing to be downloaded from this repository, and only for docker to be installed and running on your machine. In order for the user to actually run the application from this, they would:
1. Open their command line, and run the command docker pull dylanward01/frequencybandanalysis:0.9.1
2. After the image finishes pulling, run the command docker run -d --rm -p 3838:3838 dylanward01/frequencybandanalysis:0.9.1 (or replacing the 3838:3838 with whatever port you wish to run this application through).
3. Once that command finishes, and the application has loaded (which will take a few minutes after the command finishes), open up your favorite web browser, and type in localhost:3838 (or again replacing the 3838 if you chose a different port number), and proceed to run through and enjoy the shiny application.

## Within the Application

To see a detailed explanation of how to utilize the application once it's opened, and how to interpret its results, one can consult the [vignette](https://github.com/dylanward01/FrequencyBandAnalysis/blob/main/FrequencyBandAnalysis/Vignettes/Shiny-App-Usage.pdf) seen in this repository.

## Citations

Throughout this application, 3 different papers are mentioned and referred to. Their doi's, and the algorithm produced in each of those papers is mentioned below:
- (Univariate): doi.org/10.1080/01621459.2019.1671199
- (Multivariate): doi.org/10.48550/arXiv.2301.03664
- (Functional): doi.org/10.48550/arXiv.2102.01784
