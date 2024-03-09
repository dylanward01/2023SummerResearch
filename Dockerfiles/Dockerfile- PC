# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.3.1

RUN R -e "install.packages('devtools')"
RUN R -e "library(grid)"
RUN R -e "library(devtools);install_version('shiny', version = '1.7.4.1', dependencies= T);"
RUN R -e "library(devtools);install_version('momentchi2', version = '0.1.5', dependencies= T);"
RUN R -e "library(devtools);install_version('fields', version = '14.1', dependencies= T);"
RUN R -e "library(devtools);install_version('viridis', version = '0.6.4', dependencies= T);"
RUN R -e "library(devtools);install_version('signal', version = '0.7-7', dependencies= T);"
RUN R -e "library(devtools);install_version('fossil', version = '0.4.0', dependencies= T);"
RUN R -e "library(devtools);install_version('cowplot', version = '1.1.1', dependencies= T);"
RUN R -e "library(devtools);install_version('gridBase', version = '0.4-7', dependencies= T);"
RUN R -e "library(devtools);install_version('gridExtra', version = '2.3', dependencies= T);"
RUN R -e "library(devtools);install_version('ggplot2', version = '3.4.2', dependencies= T);"
RUN R -e "library(devtools);install_version('RcppArmadillo', version = '0.12.6.1.0', dependencies= T);"
RUN R -e "library(devtools);install_version('shinyjs', version = '2.1.0', dependencies= T);"
RUN R -e "library(devtools);install_version('fda', version = '6.1.4', dependencies= T);"
RUN R -e "library(devtools);install_version('plotly', version = '4.10.2', dependencies= T);"

COPY /FrequencyBandAnalysis ./app

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]