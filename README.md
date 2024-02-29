## Repository Contents

- Above, 2 things can be seen: A Dockerfile, and a folder entitled FrequencyBandAnalysis.
  - The Dockerfile is given to enable the user to create their own version of their own version of a docker container for this application. It is not needed to run the application through R, and is only needed if the user wishes to implement Docker. 
  - Within the folder FrequencyBandAnalysis one can see:
  -   The file app.R, which contains the code for actually running the shiny application.
  -   The 3 other R files (all with functions in their names), and the 2 cpp files- all of which contain the actual functions and codes for the algorithms that are run.
  -   A folder entitled "Papers", which contains the 3 papers that establish the methodology on which this application runs. This folder is not required to be present for the application to run- it is simply there to give the user easier access to the references for the algorithms and code that are implemented throughout this application.
  -   A folder entitled "Sample Data", which contains simulated data for each of the 3 types of data that this application works on. This folder is not requireed to be present for the application to run- it is simply there to give the user some sample data to explore the application, and the outputs that it gives.

## Running the Application in R

1.) Download the FrequencyBandAnalysis folder seen above- containing a file named app.R, 3 files containing functions in R, and 2 files containing functions in C++, among other things. Notably, the Dockerfile is not needed if you are going to run the app exclusively through R.
2.) Ensure that the packages listed in the app.R file are installed on your machine
3.) Proceed to run the app.R file, and wait for the application to open in a new window from your R interface.

## Within the Application

- Once the shiny application opens, the user will see 2 tabs- "File Upload" and "Simulation Setting".
-   The "File Upload" tab serves to give the user the ability to upload their own dataset, and run the corresponding algorithm on their uploaded data. Once the file is uploaded, some descriptors of the data will pop up on the left side of the application (some of which may prompt a selection), and a plot showing the data will appear in the application. Below the aforementioned descriptors, the user will be presented with a number of different parameters that they can provide a numerical input for. For these parameters, valid ranges for their values are explained in footnotes, denoted by *'s. After the necessary parameters have been given, the "Run" button can be pressed, and either the algorithm will run, or error messages that need to be resolved will pop up. As the algorithm runs, a message in the application will show- prompting the user to look to a given text file that will be created within the directory of the application, to view the progress of the currently running algorithm. Once the algorithm finishes, the message will disappear, and be replaced by a number of different plots that display the results from the algorithms. At the very bottom of the page. there is a button to click that can let you download the plots seen above into a 2-page pdf file.
    -  The dataset that is uploaded needs to be entirely numeric. There can be column names- provided that the header box is checkmarked prior to the file being uploaded. If the file is not read in as fully numeric, some message will be visible, which can give prompting on how to fix this issue.
- The "Simulation Setting" tab serves to give the user an introduction to the application, by using simulated data instead. The user can select the type of data that they want to look at (Univariate, Functional, or Multivariate), and then the type of that data they want to look at (Linear, White Noise, etc.). After those selections are made, the user can choose from 2 or 3 preset values for each of the parameters for the chosen algorithm. Again, these choices are laid out in footnotes denoted by *'s. After the necessary parameters have been given, the "Run" button can be pressed, and either the algorithm will run, or error messages that need to be resolved will pop up. As the algorithm runs, a message in the application will show- prompting the user to look to a given text file that will be created within the directory of the application, to view the progress of the currently running algorithm. Once the algorithm finishes, the message will disappear, and be replaced by a number of different plots that display the results from the algorithms. At the very bottom of the page. there is a button to click that can let you download the plots seen above into a 2-page pdf file.

## Docker

If one wants to use this application, but doesn't want to utilize R, they can run this application through Docker. There are 2 methods that the user could embark upon, to run the application this way:
  1. Firstly, they could use Dockerhub- as an image for this application is currently located at https://hub.docker.com/r/dylanward01/frequencybandanalysis. This method requires nothing to be downloaded from this repository, and only for docker to be installed and running on your machine. In order for the user to actually run the application from this, they would:
    1.1. Open their command line, and run the command docker pull dylanward01/frequencybandanalysis:5410
    1.2. After the image finishes pulling, run the command docker run -d --rm -p 3838:3838 dylanward01/frequencybandanalysis:5410 (or replacing the 3838:3838 with whatever port you wish to run this application through).
    1.3. Once that command finishes, and the application has loaded (which will take a few minutes after the command finishes), open up your favorite web browser, and type in localhost:3838 (or again replacing the 3838 if you chose a different port number), and proceed to run through and enjoy the shiny application.
  2. Secondly, they could create the docker container for this application themselves. This method requires both the contents of this repository to be downloaded, as well as for Docker to be installed and running on your machine. To do this, the yser would:
    2.1.  Open up their command line, and set the working directory to be the folder containing both the Dockerfile, and the subfolder entitled FrequencyBandAnalysis (ie, the folder you got from properly downloading this repository).
    2.2.  In their command line, run the command docker build .
    2.3.  After that command finishes, run the command docker images, to see the image ID of the image that was just created, and then run docker run -d --rm -p 3838:3838 (Insert entire Image ID here).
    2.4.  Once that command finishes, and the application has loaded (which will take a few minutes after the command finishes), open up your favorite web browser, and type in localhost:3838, and proceed to run through and enjoy the shiny application.
