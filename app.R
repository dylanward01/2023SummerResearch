library(shiny)
library(momentchi2)
library(fields)
library(viridis)
library(signal)
library(fossil)
library(cowplot)
library(grid)
library(gridBase)
library(gridExtra)
library(ggplot2)
library(RcppArmadillo)
library(shinyjs)
library(fda)
library(plotly)


source("EBA_functions.R")
source("fEBA_Rfns.R")
Rcpp::sourceCpp("fEBA_072321.cpp")
df <- read.csv("xtsample.csv")
df <- df[[1]]

ui <- fluidPage(
  useShinyjs(),
  tags$h1("Empirical Frequency Band Analysis"),
  tags$a(tags$strong(tags$em("Source: Empirical Frequency Band Analysis of Nonstationary Time Series")), 
         href = "https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1671199"),
  tags$p(tags$em("Research reported in this publication was supported by the National Institute of General Medical Sciences 
         of the National Institutes of Health under Award Number R01GM140476. The content is solely the responsibility
         of the authors and does not necessarily represent the official views of the National Institutes of Health")),
  tags$p(tags$strong("Authors: Dylan Ward, Mohit Chhaparia, Kevin Gunalan Antony Michael Raj")), (tags$em(tags$u("Under the guidance of Professor Scott Bruce"))),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(condition = "input.tabselected==1",
                       radioButtons("type", "Time Series Type", 
                                    choices=c("Univariate","Functional", "Multivariate"), 
                                    selected="Univariate"), 
                       conditionalPanel(condition = "input.type == 'Univariate'", 
                       selectInput("Simsetting", "Simulation Setting",
                                   c("White Noise" = "W",
                                     "Linear" = "L",
                                     "Sinusoidal" = "S"), selected="S"),
                       selectInput(inputId="Time", label = "Choose total length of time series (T)", 
                                   choices = as.numeric(c(500,1000,5000,10000, 20000, 50000)), selected=1000), 
                       selectInput(inputId="Num", label = "Choose number of observations per approximately stationary block* (N)", 
                                   choices=100, selected=100),
                       selectInput(inputId="Tapers", label="Choose number of tapers to use in multitaper spectral estimator** (K)", 
                                  choices=10, selected=10),
                      selectInput(inputId = "Signi", label="Choose significance level", 
                                  choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05),
                       radioButtons(inputId = "TF", label = "Standardize", 
                                    c("True" = TRUE, "False" = FALSE), selected = FALSE),
                       actionButton("go", label = "Run"),
                      htmlOutput("Res"),
                      htmlOutput("Res1"),
                      ), 
                      conditionalPanel(condition = "input.type== 'Functional'", 
                                       radioButtons("Plot3D", "3D Plots", 
                                                    choices=c("Include","Exclude"), 
                                                    selected="Include"), 
                      selectInput("SimF1", "Simulation Setting",
                                  c("White Noise" = "W",
                                  "Linear" = "L",
                                  "Sinusoidal" = "S"), selected="S"),
                      selectInput(inputId = "TsF1", label="Choose total length of time series (T)", 
                                  choices = c( 500, 1000, 2000, 5000), selected=2000), 
                      selectInput(inputId = "RF1", label = "Choose number of points in functional domain (R)", 
                                  choices=seq(from=5, to=50, by=5), selected=5), 
                      selectInput(inputId = "NF1", label = "Choose number of observations per approximately stationary block* (N)",
                                  choices = 30, selected = 30), 
                      selectInput(inputId = "KF1", label = "Choose number of tapers to use in multitaper spectral estimator** (K)", 
                                  choices = 5, selected = 5),
                      selectInput(inputId = "RselF1", label = "Choose number of points in the functional domain to use in computing test statistics*** (Rsel)", 
                                  choices = c(5,10), selected = 5), 
                      selectInput(inputId = "AlphaF1", label="Choose significance level", 
                                  choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05),
                      radioButtons(inputId = "TF_F1", label = "Standardize", 
                                   c("True" = TRUE, "False" = FALSE), selected = FALSE),
                      actionButton("goF1", label = "Run"),
                      htmlOutput("F1_1"),
                      htmlOutput("F1_2"),
                      htmlOutput("F1_3")
                      ),
                      conditionalPanel(condition="input.type == 'Multivariate'", 
                                       htmlOutput("Check111"),
                      )),
      
      conditionalPanel(condition = "input.tabselected==2",
                       fileInput("file_csv", "Choose CSV File",
                                 multiple = TRUE,
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       checkboxInput("header", "Header", TRUE),
                       htmlOutput("UniVarDis"),
                       htmlOutput("NotUniVar"),
                       htmlOutput("BlankSpace"),
                       hidden(radioButtons(inputId = "Data_Checker", label = NULL, choices = c("Functional", "Multivariate"), selected=character(0))),
                       hidden(radioButtons("Plot3D_File", "3D Plots", 
                                    choices=c("Include","Exclude"), 
                                    selected="Include")),
                       htmlOutput("T_len"),
                       hidden(htmlOutput("Ts_Fxn_Dim")),
                       numericInput(inputId = "Num2", label = "Choose number of observations per approximately stationary block* (N)", 
                                   value = NULL, step = 1),
                       
                       numericInput(inputId = "Tapers2", label = "Choose number of tapers to use in multitaper spectral estimator** (K)", 
                                   value = NULL, step = 1),
                       
                       selectInput(inputId = "Signi2", label="Choose significance level", 
                                   choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05),
                       
                       radioButtons(inputId = "TF2", label = "Standardize", 
                                    c("True" = TRUE, "False" = FALSE), selected = FALSE),
                       hidden(numericInput(inputId = "Num_Fxna", label = "Choose number of observations per approximately stationary block* (N)", 
                                           value = NULL, step = 1)),
                       hidden(numericInput(inputId = "Tapers_Fxna", label = "Choose number of tapers to use in multitaper spectral estimator** (K)", 
                                           value = NULL, step = 1)),
                       hidden(selectInput(inputId = "Rsel_Fxna", label = "Choose number of points in the functional domain to use in computing test statistics*** (Rsel)", 
                                          choices = c(5,10), selected = 5)),
                       hidden(selectInput(inputId = "Signi_Fxna", label="Choose significance level", 
                                   choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05)),
                       
                       hidden(radioButtons(inputId = "TF_Fxna", label = "Standardize", 
                                    c("True" = TRUE, "False" = FALSE), selected = FALSE)),
                       hidden(actionButton("go_Fxna", label = "Run")),
                       actionButton("go2", label = "Run"), 
                       htmlOutput("res9"),
                       htmlOutput("res10"),
                       hidden(htmlOutput("Fxn_AA")),
                       hidden(htmlOutput("Fxn_BB")),
                       hidden(htmlOutput("Fxn_CC"))
                       ), 
     
                      
    ),
    mainPanel(
      tabsetPanel(type = "tabs", id = 'tabselected', selected = 1,
                  tabPanel("Simulation Setting", value = 1),
                  tabPanel("File Upload", value = 2)),
      conditionalPanel(condition = "input.tabselected==1",
                       conditionalPanel(condition = "input.type == 'Univariate'",
                       plotOutput("Image_Plot", height=1000, width=1000),
                       fluidRow(
                         column(width = 6, plotOutput("summ_out_uni", height = 500, width=500)), 
                         column(width = 6, plotOutput("summ_pval_uni", height = 500, width = 500)), 
                         #column(width = 4, dataTableOutput("summ_flat_uni"))
                       ),
                       downloadButton('downloadData','Download the Above Results') ,
                       br(),
                       br(),
                       ),
                       
                       conditionalPanel(condition = "input.type == 'Functional'",
                       plotlyOutput("Fxn_Plota", height=400, width=1000),
                       fluidRow(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                column(width=6,  htmlOutput("FxnPlotaDesc") ),
                                column(width=4, hidden(htmlOutput("test12121"))),
                                column(width=2, hidden(sliderInput(inputId = "x_F1", min=1, max=10, step=1, value=1,label=NULL, ticks = FALSE)))
                       ),
                       
                     
                       tags$head(tags$style(HTML('.irs-from, .irs-min, .irs-to, .irs-max, .irs-single {
                                                 visibility: hidden !important;
                                                 }' ))),
                       
                       plotOutput("Fxn_Plotb", height=600, width = 1000), 
                       fluidRow(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                column(width=6,  htmlOutput("FxnbPlotDesc") ),
                                column(width=4, htmlOutput("Blank2")),
                                column(width=2, hidden(sliderInput("plot1_FxnCheck",min=1,max=10,step=1,value=1,label=NULL, ticks = FALSE)))
                       ), 
                       tags$head(tags$style(HTML('.irs-from, .irs-min, .irs-to, .irs-max, .irs-single {
                                                 visibility: hidden !important;
                                                 }' ))),
                       conditionalPanel(condition = "input.Plot3D == 'Include'",
                       plotlyOutput("Plotly_Fxna", height=600, width=1000),
                       plotlyOutput("Plotly_Fxnb", height=600, width=1000), 
                       fluidRow(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                column(width=6,  htmlOutput("FxnPlot22Desc") ),
                                column(width=4, htmlOutput("Blank10100")),
                                column(width=2, hidden(sliderInput(inputId = "q11_F1", min=1, max=10, step=1, value=1, label=NULL, width="125%", ticks=FALSE)))
                       )
                       
                       )
      )),
      conditionalPanel(condition = "input.tabselected==2",
                       plotOutput("Image_Plota", height=400, width=1000),
                       plotOutput("Image_Plot2", height=600, width=1000), 
                       fluidRow(
                         column(width = 6, plotOutput("summ_out_uni_file", height = 500, width=500)), 
                         column(width = 6, plotOutput("summ_pval_uni_file", height = 500, width = 500)), 
                         #column(width = 4, dataTableOutput("summ_flat_uni"))
                       ),
                       downloadButton('downloadData1','Download the Above Results'),
                       br(),
                       br(),
                       conditionalPanel(condition = "input.Data_Checker== 'Functional'", 
                                        plotlyOutput("Test_Fxna_Plot1", width=1000), 
                                        fluidRow(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                          column(width=6,  htmlOutput("FxnPlotaDesc_AA") ),
                                          column(width=4, htmlOutput("test12121_AA")),
                                          column(width=2, sliderInput(inputId = "x_F1_AA", min=1, max=10, step=1, value=1,label=NULL, ticks = FALSE))
                                        ),
                                        tags$head(tags$style(HTML('.irs-from, .irs-min, .irs-to, .irs-max, .irs-single {
                                                 visibility: hidden !important;
                                                 }' ))),
                                        plotOutput("Fxn_Plotb_file", height=600, width=1000), 
                                        fluidRow(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                                 column(width=6,  htmlOutput("FxnbPlotDesc_file") ),
                                                 column(width=4, htmlOutput("Blank2_file")),
                                                 column(width=2, sliderInput("plot1_FxnCheck_file",min=1,max=10,step=1,value=1,label=NULL, ticks = FALSE))
                                        ),
                                        tags$head(tags$style(HTML('.irs-from, .irs-min, .irs-to, .irs-max, .irs-single {
                                                 visibility: hidden !important;
                                                 }' ))),
                                        conditionalPanel(condition = "input.Plot3D_File == 'Include'",
                                        plotlyOutput("Plotly_Fxna_file", height=600, width=1000),
                                        plotlyOutput("Plotly_Fxnb_file", height=600, width=1000), 
                                        fluidRow(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                                 column(width=6,  htmlOutput("FxnPlot22Desc_file") ),
                                                 column(width=4, htmlOutput("Blank10100_file")),
                                                 column(width=2, sliderInput(inputId = "q11_F1_file", min=1, max=10, step=1, value=1, label=NULL,  ticks=FALSE))
                                        ),
                                        
                                        )
                                        
                       )
      ), 
      )
  ))

server <- function(input,output, session) {
  
  output$BlankSpace <- renderText({
    paste("\n")
  })
  observe({
    RF1 <- as.numeric(input$RF1)
    if(RF1 == 5){
      updateSelectInput(session, "RselF1", choices=c(5), selected=5)
    } else {
      updateSelectInput(session, "RselF1", choices=c(5, 10), selected=5)
    }
  })
  observe({
    TF1 <- as.numeric(input$TsF1)
    new_vals <- floor(seq(from=30, to=TF1, length.out = 20))
    updateSelectInput(session, "NF1", choices=new_vals, selected=new_vals[2])
    
  })
  observe({
    NF1 <- as.numeric(input$NF1)
    if(is.na(NF1)){ 
    } else {
      new_vals <- floor(seq(from=1, to=floor(NF1/4 - 1), length.out = 5))
      updateSelectInput(session, "KF1", choices = new_vals, selected=new_vals[3])  
    }
    
  })
  observe({
    t2 <- as.numeric(input$TsF1)
    if(t2 == 500) {
      choices = c(62, 105)
    } else if(t2 == 1000) {
      choices = c(31, 100, 177)
    } else if(t2 == 5000) {
      choices = c(70, 292, 594 )
    } else if(t2 == 10000) {
      choices = c(100, 464, 1000)
    } else if(t2==2000) {
      choices = c(44, 158, 299)
    } else {
      choices=c(223, 1357, 3343)
    }
    updateSelectInput(session, "NF1", choices = choices, selected = choices[length(choices) - 1])
  })
  observe({
    t2 <- as.numeric(input$TsF1)
    sc_val <- (as.numeric(input$NF1))
    if(t2 == 500) {
      if(sc_val == 105){
        tap=10
      } else {
        tap=7
      }
    } else if(t2 == 1000) {
      if(sc_val == 31){
        tap=5
      } else if(sc_val == 100) {
        tap=10
      } else {
        tap=13
      }
    } else if(t2 == 5000) {
      if(sc_val == 70){
        tap=8
      } else if(sc_val == 292) {
        tap = c(17, 70)
      } else {
        tap=c(24, 120)
      }
    } else if(t2 == 10000) {
      if(sc_val == 100) {
        tap=10
      } else if(sc_val == 464){
        tap=c(21, 99)
      } else {
        tap=c(31, 177)
      }
    } else if(t2 == 2000){
      if(sc_val == 44) {
        tap=6
      } else if(sc_val == 158) {
        tap=c(12, 44)
      } else {
        tap=c(17, 71)
      }
    } else {
      if(sc_val == 223){
        tap=14
      } else if(sc_val == 1357){
        tap=c(36,223)
      } else {
        tap=c(57, 439)
      }
    }
    updateSelectInput(session, "KF1",choices=tap, selected = tap[1])
  })
  observe({
    t2 <- as.numeric(input$Time)
    if(t2 == 500) {
      choices = c(62, 105)
    } else if(t2 == 1000) {
      choices = c(31, 100, 177)
    } else if(t2 == 5000) {
      choices = c(70, 292, 594 )
    } else if(t2 == 10000) {
      choices = c(100, 464, 1000)
    } else if(t2==20000) {
      choices = c(141, 736, 1681)
    } else {
      choices=c(223, 1357, 3343)
    }
    updateSelectInput(session, "Num", choices = choices, selected = choices[length(choices) - 1])
  })
  observe({
    t2 <- as.numeric(input$Time)
    sc_val <- (as.numeric(input$Num))
    if(t2 == 500) {
      if(sc_val == 105){
        tap=10
      } else {
        tap=7
      }
    } else if(t2 == 1000) {
      if(sc_val == 31){
        tap=5
      } else if(sc_val == 100) {
        tap=10
        } else {
        tap=13
      }
    } else if(t2 == 5000) {
      if(sc_val == 70){
        tap=8
      } else if(sc_val == 292) {
        tap = c(17, 70)
      } else {
        tap=c(24, 120)
      }
    } else if(t2 == 10000) {
      if(sc_val == 100) {
        tap=10
      } else if(sc_val == 464){
        tap=c(21, 99)
      } else {
        tap=c(31, 177)
      }
    } else if(t2 == 20000){
      if(sc_val == 141) {
        tap=11
      } else if(sc_val == 736) {
        tap=c(27, 141)
      } else {
        tap=c(41, 262)
      }
    } else {
      if(sc_val == 223){
        tap=14
      } else if(sc_val == 1357){
        tap=c(36,223)
      } else {
        tap=c(57, 439)
      }
    }
    updateSelectInput(session, "Tapers",choices=tap, selected = tap[1])
  })
  output$Res <- renderText({
    paste(h6("*Choices are T", HTML(paste(tags$sup("1/2"))), ", T", HTML(paste(tags$sup("2/3"), 
                  ", or T", HTML(paste(tags$sup("3/4"), ", provided it satisfies 30", HTML("&le;"), 
                                    "N", HTML("&le;"), HTML(paste(tags$sup("T"))), "/", HTML(paste(tags$sub(2)))))))))
    })
  output$Res1 <- renderText({
    paste(h6("**Choices are N", HTML(paste(tags$sup("1/2"))), ", or N", HTML(paste(tags$sup("3/4"), 
                   ", provided it satisfies ", HTML(paste(tags$sup("floor(N/2)"))), 
                                     "/", HTML(paste(tags$sub(2))), " - 1> floor((K+1)(", HTML(paste(tags$sup("N"))), 
                                     "/", HTML(paste(tags$sub("N+1"))), "))"))))
  })
  output$res9 <- renderText({
    paste("*Please enter a dataframe")
  })
  output$F1_1 <- renderText({
    paste(h6("*Choices are T", HTML(paste(tags$sup("1/2"))), ", T", HTML(paste(tags$sup("2/3"))),
             ", or T",HTML(paste(tags$sup("3/4"))),", provided it satisfies 30 ", HTML("&le;"), "N", HTML("&le;"), "T"))
  })
  output$F1_2 <- renderText({
    paste(h6("**Choices are N", HTML(paste(tags$sup("1/2"))),", or N",HTML(paste(tags$sup("3/4"))),
             ", provided it satisfies 1 ", HTML("&le;"), "K", "<", "floor(N/4 - 1)"))
  })
  output$F1_3 <- renderText({
    paste(h6("***Choices are 5, or 10, provided it satisfies 1 ", HTML("&le;"), "Rsel", HTML("&le;"), "R"))
  })
  Freq_domain_vals <- eventReactive(input$RF1, {
    R_val <- as.numeric(input$RF1)
    domain <- numeric(0)
    if(is.na(R_val)){
      
    } else if(R_val == 5){
      domain <- c(1, 3, 5)
      nums <- seq(1:5)
      ifelse(domain %in% nums, domain, " ")
    } else if (R_val == 10){
      domain <- c(1, 3, 5, 7, 10)
    }else if (R_val == 15){
      domain <- c(1, 3, 5, 7, 10, 13, 15)
    }else if (R_val == 20){
      domain <- c(1, 3, 5, 7, 10, 13, 15, 17, 20)
    }else if (R_val == 25){
      domain <- c(1, 3, 5, 7, 10, 13, 15, 17, 20, 23, 25)
    }else if (R_val == 30){
      domain <- c(1, 3, 5, 7, 10, 13, 15, 17, 20, 23, 25, 27, 30)
    }else if (R_val == 35){
      domain <- c(1, 3, 5, 7, 10, 13, 15, 17, 20, 23, 25, 27, 30, 33, 35)
    }else if (R_val == 40){
      domain <- c(1, 3, 5, 7, 10, 13, 15, 17, 20, 23, 25, 27, 30, 33, 35, 37, 40)
    }else if (R_val == 45){
      domain <- c(1, 3, 5, 7, 10, 13, 15, 17, 20, 23, 25, 27, 30, 33, 35, 37, 40, 43, 45)
    }else if (R_val == 50){
      domain <- c(1, 3, 5, 7, 10, 13, 15, 17, 20, 23, 25, 27, 30, 33, 35, 37, 40, 43, 45, 47, 50)
    }
    (list("domain" = domain))
  })
  plot.listF1 <- eventReactive(input$goF1, ignoreNULL = TRUE, {
    ##get sim data
    nb=15; #number of basis functions used to generate white noise
    R=as.numeric(input$RF1); #number of points in functional domain
    Ts=as.numeric(input$TsF1); #length of time series
    seed=234; #seed for reproducibility
    B=floor(as.numeric(input$TsF1)/as.numeric(input$NF1)); #number of time blocks
    N=as.numeric(input$NF1); #number of observations per time block
    bw=floor((as.numeric(input$KF1) + 1) / (as.numeric(input$NF1) + 1)); #bandwidth for multitaper spectral estimator
    K=as.numeric(input$KF1); #number of tapers for multitaper spectral estimator
    std=as.numeric(input$TF_F1); #standardize variance for points in functional domain (TRUE) or not (FALSE)
    freq=seq(from=0,by=1/N,length.out=floor(N/2)+1); #Fourier frequencies
    Rsel=as.numeric(input$RselF1); #number of points in functional domain used for test statistics
    
    if (input$SimF1 == "W"){
      X=fws.sim(nb=nb,gsz=R,Ts=Ts,seed=seed);
      pse=fhat(X,N,K,Rsel,std);
       
       cmpnt="1-1"; #select component to view
       dimnames(pse) <- list(freq,apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-"),1:B);
       plot.x <- floor((1:B) * (Ts/B)); plot.y <- as.numeric(rownames(pse)); plot.z <- pse
       plot.cmp <- apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-")
       plot.main <- "Multitaper Autospectrum"; plot.data = X
     
      
    } else if (input$SimF1 == "L") {
      X=f3bL.sim(nb=nb,gsz=R,Ts=Ts,seed=seed);
      pse=fhat(X,N,K,Rsel,std);
       cmpnt="1-1"; #select component to view
       dimnames(pse) <- list(freq,apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-"),1:B);
       plot.x <- floor((1:B) * (Ts/B)); plot.y <-as.numeric( rownames(pse)); plot.z <- pse
       plot.cmp <- apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-")
       plot.main <- "Multitaper Autospectrum"; plot.data = X
       
      
    } else if (input$SimF1 == "S") {
      X=f3bS.sim(nb=nb,gsz=R,Ts=Ts,seed=seed);
       pse=fhat(X,N,K,Rsel,std);
       cmpnt="1-1"; #select component to view
       dimnames(pse) <- list(freq,apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-"),1:B);
       plot.x <- floor((1:B) * (Ts/B)); plot.y <- as.numeric(rownames(pse)); plot.z <- pse
       plot.cmp <- apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-")
       plot.main <- "Multitaper Autospectrum"; plot.data = X
       
    }
    output$Fxn_Plota <- renderPlotly({
         a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(X[1,])), y=X[1,])) + 
           xlab("Functional Domain") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
           scale_x_continuous(limits=c(0,1), expand=c(0,0))
         ggplotly(a)
       })
    conf <- numeric(length(pse))
    dim(conf) <- dim(pse)
    dimnames(conf) <- dimnames(pse)
    for(k in 1:dim(pse)[3]){
      for(j in 1:dim(pse)[2]){
        first_col <- ((j - 1) %% (as.numeric(input$RselF1))) + 1 
        second_col <- ((j - 1) %/% (as.numeric(input$RselF1))) + 1
        cmpt_1 <- paste(first_col, "-", first_col, sep="")
        cmpt_2 <- paste(second_col, "-", second_col, sep="")
        for(i in 1:dim(pse)[1]){
          if(first_col == second_col){
            conf[i,j,k] <- Re(pse[i,j,k])
          } else {
            conf[i,j,k] <- Re((Mod(pse[i,j,k])**2) / (pse[i,cmpt_1,k] * pse[i,cmpt_2,k]))
          }
        }
      }
    }
    plot.z <- conf
    
    output$FxnPlotaDesc <- renderText({
      paste(h4("Currently viewing timepoint "))
    })
    
    output$FxnbPlotDesc <- renderText({
      paste(h4("Currently viewing component "))
    })
    output$FxnPlot22Desc <- renderText({
      paste(h4("Currently viewing component "))
    })
    show("x_F1")
    show("q_F1")
    show("q11_F1")
    show("Plotly_Fxna")
    show("plot1_FxnCheck")
    show("test12121")
    show("Blank2")
    show("Blank10100")
    updateSliderInput(session, "x_F1", min = 1, max=as.numeric(input$TsF1), value=1, step=1)
    updateSelectInput(session, "q_F1", choices=plot.cmp, selected = plot.cmp[1])
    updateSliderInput(session, "q11_F1", min=1, max=length(plot.cmp), value=1)
    updateSliderInput(session, "plot1_FxnCheck", min=1, max=length(plot.cmp), value=1)
    output$Blank2 <- renderText({
      paste(h4(strong((paste("1-1")))))
    })
    output$Blank10100 <- renderText({
      paste(h4(strong((paste("1-1")))))
    })
    output$test12121 <- renderText({
      paste(h4(strong((paste("1")))))
    })
    output$Plotly_Fxna <- renderPlotly({
      plot_ly(x = ~seq(from=0, to=1, length.out = ncol(X)), 
              y = ~seq(from=1, to=nrow(X)), 
              z = ~X) %>% add_surface() %>% layout(
        title = "3D Representation of Simulated Data",
        scene = list(
        xaxis = list(title="Functional Domain"), 
        yaxis = list(title = "Timepoint"), 
        zaxis = list(title="Value")
      )) %>% colorbar(title = "Value", len=1)
    })
    show("Fxn_Row")
    plot.main_2 <- "Estimated Coherence"
    plot.main_3D <- "3D Representation of Autospectrogram"
    plot.main_3D.2 <- "3D Representation of Coherence"
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
        plot.main = plot.main, plot.cmp = plot.cmp, plot.data = plot.data, 
        plot.main_2 = plot.main_2, plot.main_3D = plot.main_3D, plot.main_3D.2 = plot.main_3D.2)
    
  });
  
  output$Fxn_Plotb <- renderPlot({
    image.plot(x=plot.listF1()[[1]],y=plot.listF1()[[2]],z=suppressWarnings(t(Re(plot.listF1()[[3]][,"1-1",]))), 
               axes = TRUE, col = inferno(256), 
               main = plot.listF1()[[4]],xlab='Time',ylab='Hz',xaxs="i"); 
  })
  
  output$Plotly_Fxnb <- renderPlotly({
    plot_ly(y=~plot.listF1()[[1]], x=~plot.listF1()[[2]], z=~t(Re(plot.listF1()[[3]][,"1-1",])))  %>%layout(title=plot.listF1()[[8]], 
                                                                                                            scene = list( 
           xaxis = list(title='Frequency',range = c(0.5, 0)), 
           yaxis = list(title="Timepoint"), 
           zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
  })
  observeEvent(input$plot1_FxnCheck, ignoreNULL = TRUE, {
    curr_num <- as.numeric(input$plot1_FxnCheck)
    if(is.na(plot.listF1()[[4]])){
      
    } else {
      curr_comp <- plot.listF1()[[5]][curr_num]
      output$Blank2 <- renderText({
        paste(h4(strong((paste(curr_comp)))))
      })
      if(strsplit(curr_comp, "-")[[1]][1] == strsplit(curr_comp, "-")[[1]][2]){
        output$Fxn_Plotb <- renderPlot({
          image.plot(x=plot.listF1()[[1]],y=plot.listF1()[[2]],z=suppressWarnings(t(Re(plot.listF1()[[3]][,curr_comp,]))), 
                     axes = TRUE, col = inferno(256), 
                     main = plot.listF1()[[4]],xlab='Time',ylab='Hz',xaxs="i")
        })
          
      } else {
        output$Fxn_Plotb <- renderPlot({
          image.plot(x=plot.listF1()[[1]],y=plot.listF1()[[2]],z=suppressWarnings(t(Re(plot.listF1()[[3]][,curr_comp,]))), 
                     axes = TRUE, col = inferno(256), 
                     main = plot.listF1()[[7]],xlab='Time',ylab='Hz',xaxs="i"); 
        })
          
      }
        
    
    }
    
  })
   
  observeEvent(input$q11_F1, ignoreNULL = TRUE, {
    curr_num <- as.numeric(input$q11_F1)
    if(is.na(plot.listF1()[[4]])){

    } else {
      curr_comp <- plot.listF1()[[5]][curr_num]
      output$Blank10100 <- renderText({
        paste(h4(strong((paste(curr_comp)))))
      })
      vals <- as.numeric(strsplit(curr_comp, "-")[[1]])
      if(vals[1] == vals[2]){
        output$Plotly_Fxnb <- renderPlotly({
          plot_ly(y=~plot.listF1()[[1]], x=~plot.listF1()[[2]], z=~t(Re(plot.listF1()[[3]][,curr_comp,])))  %>%layout(title=plot.listF1()[[8]],
                                                                                                                      scene = list(
                                                                                                                        xaxis = list(title='Frequency',range = c(0.5, 0)),
                                                                                                                        yaxis = list(title="Timepoint"),
                                                                                                                        zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
          
          
        }) 
      } else {
        output$Plotly_Fxnb <- renderPlotly({
          plot_ly(y=~plot.listF1()[[1]], x=~plot.listF1()[[2]], z=~t(Re(plot.listF1()[[3]][,curr_comp,])))  %>%layout(title=plot.listF1()[[9]],
                                                                                                                      scene = list(
                                                                                                                        xaxis = list(title='Frequency',range = c(0.5, 0)),
                                                                                                                        yaxis = list(title="Timepoint"),
                                                                                                                        zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
          
          
        })
      }
     
                                                                                                              

    }

  }) 
  # observeEvent(input$Fxn_Row, ignoreNULL = TRUE, {
  #   if(input$Fxn_Row == 'No'){
  #     hide("Plotly_Fxna.5")
  #     hide("x11_F1")
  #     hide("FxnPlot11Desc")
  #   } else {
  #     show("Plotly_Fxna.5")
  #     show("x11_F1")
  #     updateSelectInput(session, "x11_F1", choices = seq(from=1, to = dim(plot.listF1()[[6]])[1], by=1), selected=1)
  #     show("FxnPlot11Desc")
  #     output$Plotly_Fxna.5 <- renderPlot({
  #       ggplot() + geom_line(aes(x=seq(from=1, to=as.numeric(input$RF1), length.out=length(plot.listF1()[[6]][1,])), y=plot.listF1()[[6]][1,])) + 
  #         xlab("Functional Domain") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) 
  #     })
  #     output$FxnPlot11Desc <- renderText({
  #       paste(h4("Currently viewing timepoint "))
  #     })
  #     
  #   }
  # })
  # 
  plot.listaa <- eventReactive(input$file_csv, ignoreNULL = FALSE,  {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    le <- length(file[[1]])
    if(le == 0){
      
    } else {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    dataf <- read.csv(file$datapath, header = input$header)
    dataf <- dataf[[1]]
    list(dataf=dataf)}}
    )
  plot.listbb <- eventReactive(input$file_csv, ignoreNULL = FALSE,  {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    le <- length(file[[1]])
    if(le == 0){
      
    } else {
      file <- input$file_csv
      ext <- tools::file_ext(file$datapath)
      dataf <- read.csv(file$datapath, header = input$header)[,-1]
      updateNumericInput(session, "Num_Fxna", value  = floor(sqrt(dim(dataf)[1])))
      updateNumericInput(session, "Tapers_Fxna", value = floor(dim(dataf)[1] ** 0.25))
      if(dim(dataf)[2] == 5){
        updateSelectInput(session, "Rsel_Fxna", choices = 5, selected = 5)
      }
      main = c("check")
      updateSliderInput(session, "x_F1_AA", max=dim(dataf)[1], min=1, value=1, step=1)
      list(dataf=dataf, main = main)}}
  )
  observeEvent(input$file_csv, {
    if(is.na(plot.listbb()[[2]])){

    } else {
      output$Test_Fxna_Plot1<- renderPlotly({
        a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=dim(plot.listbb()[[1]])[2]), y=as.numeric(plot.listbb()[[1]][1,]))) +
          xlab("Functional Domain") + ylab("") + ggtitle("Observed Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) +
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
        ggplotly(a)
      })
      output$test12121_AA <- renderText({
        paste(h4(strong("1")))
      })
      hide("Blank2_file")
      hide("FxnbPlotDesc_file")
      hide("plot1_FxnCheck_file")
      hide("Fxn_Plotb_file")
      output$Fxn_Plotb_file <- renderPlot({
        
      })
      hide("FxnPlot22Desc_file")
      hide("Blank10100_file")
      hide("q11_F1_file")
      output$Plotly_Fxna_file <-renderPlotly({
        
      })
      output$Plotly_Fxnb_file <- renderPlotly({
          
        })
      
      

    }
  })
  
  observeEvent(input$go_Fxna, {
    show("Blank2_file")
    show("FxnbPlotDesc_file")
    show("plot1_FxnCheck_file")
    show("Fxn_Plotb_file")
    show("FxnPlot22Desc_file")
    show("Blank10100_file")
    show("q11_F1_file")
  })
  observeEvent(input$file_csv, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    #dataf <- read.csv(file$datapath, header = input$header)[,-1]
    dataf <- read.csv(file$datapath, header = input$header)
    dims <- dim(dataf)[2]
    if(dims == 1){
      hide("NotUniVar")
      show("UniVarDis")
      output$UniVarDis <- renderText({
        paste(strong("This is a Univariate Time Series"))
      })
      hide("Data_Checker")
      updateRadioButtons(session, "Data_Checker", selected=character(0))
      hide("Plot3D_File")
      hide("Num_Fxna")
      hide("Tapers_Fxna")
      hide("Ts_Fxn_Dim")
      hide("Rsel_Fxna")
      hide("Fxn_AA")
      hide("Fxn_BB")
      hide("Fxn_CC")
      hide("Signi_Fxna")
      hide("TF_Fxna")
      hide("go_Fxna")
      show("T_len")
      show("Num2")
      show("Tapers2")
      show("Signi2")
      show("TF2")
      show("go2")
      show("res9")
      show("res10")
      show("Image_Plota")
      show("Image_Plot2")
      show("downloadData1")
      show("downloadData1_A1")
      show("summ_out_uni_file")
      show("summ_pval_uni_file")
    } else {
      dataf <- dataf[,-1]
      hide("UniVarDis")
      show("NotUniVar")
      output$NotUniVar <- renderText({
        paste(strong("This time series has multiple components. Choose what type of 
                     time series this is. "))
      })
      output$Ts_Fxn_Dim <- renderText({
        paste(strong(paste("Dimensions of Time Series (T x R): ", dim(dataf)[1], "x", dim(dataf)[2] )))
      })
      if(dim(dataf)[2] < 10){
        output$Fxn_CC <- renderText({
          paste(h6("***Valid choices are 5 as we need to satisfy 1", HTML("&le;"), 
                   "Rsel", HTML("&le;")," R"))
        })  
      } else {
        output$Fxn_CC <- renderText({
          paste(h6("***Valid choices are 5, and 10 as we need to satisfy 1", HTML("&le;"), 
                   "Rsel", HTML("&le;")," R"))
        })
      }
      
      show("Data_Checker")
      #show("Plot3D_File")
      hide("T_len")
      hide("Num2")
      hide("Tapers2")
      hide("Signi2")
      hide("TF2")
      hide("go2")
      hide("res9")
      hide("res10")
      hide("Image_Plota")
      hide("Image_Plot2")
      hide("downloadData1")
      hide("downloadData1_A1")
      hide("summ_out_uni_file")
      hide("summ_pval_uni_file")
    }
  })
  observeEvent(input$Data_Checker, {
    if(input$Data_Checker == "Functional"){
      show("Plot3D_File")
      show("Num_Fxna")
      show("Tapers_Fxna")
      show("Rsel_Fxna")
      show("Signi_Fxna")
      show("TF_Fxna")
      show("go_Fxna")
      show("Ts_Fxn_Dim")
      
      show("Fxn_AA")
      output$Fxn_AA <- renderText({
        paste(h6("*Valid choices range from 30 to ", dim(plot.listbb()[[1]])[1], "as we need to satisfy 30", HTML("&le;"), 
                 "N", HTML("&le;"), "T"))
      })
      show("Fxn_BB")
      output$Fxn_BB <- renderText({
        paste(h6("**Valid choices range from 1 to ", floor(sqrt(dim(plot.listbb()[[1]])[1]) / 4 - 1), "as we need to satisfy 1", HTML("&le;"), 
                 "K < floor(N/4 - 1)"))
      })
      show("Fxn_CC")
      # if(dim(plot.listbb()[[1]])[2] == 5){
      #   output$Fxn_CC <- renderText({
      #     paste(h6("***Valid choices are 5 as we need to satisfy 1", HTML("&le;"), 
      #              "Rsel", HTML("&le;")," R"))
      #   })  
      # } else {
      #   output$Fxn_CC <- renderText({
      #     paste(h6("***Valid choices are 5, and 10 as we need to satisfy 1", HTML("&le;"), 
      #              "Rsel", HTML("&le;")," R"))
      #   })
      # }
      
    } else {
      hide("Num_Fxna")
      hide("Tapers_Fxna")
      hide("Rsel_Fxna")
      hide("Signi_Fxna")
      hide("TF_Fxna")
      hide("go_Fxna")
      hide("Ts_Fxn_Dim")
    }
  })
  observeEvent(input$file_csv, {
    
   output$T_len <- renderText({
     paste(strong("Total Length of Time Series (T): ", length(plot.listaa()[[1]])))
   })
   output$Blank <- renderText({
     paste("")
   })
    output$Image_Plota <- renderPlot({
      print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.listaa()[[1]])), y= plot.listaa()[[1]])) + xlab("Time") +
               ylab("") + ggtitle("Observed Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
               scale_x_continuous(limits=c(0,1), expand=c(0,0)))
      
    })
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    le <- length(file[[1]])
    if(le == 0){
      
    } else {
      file <- input$file_csv
      ext <- tools::file_ext(file$datapath)
      dataf <- read.csv(file$datapath, header = input$header)
      dataf <- dataf[[1]]
      
    updateNumericInput(session, "Num2", value=floor(sqrt(length(dataf))))
    updateNumericInput(session, "Tapers2", value=floor(0.15 * sqrt(length(dataf))))
    output$res9 <- renderText({
        paste(h6("*Valid choices range from 30 to ", floor(length(dataf)/2), "as we need to satisfy 30", HTML("&le;"), 
                                                                                                   "N", HTML("&le;"), HTML(paste(tags$sup("T"))), "/", HTML(paste(tags$sub(2)))))
    })
    
    output$res10 <- renderText({
          paste("**Valid choices range from 1 to ", floor(sqrt(length(dataf))*0.24))
    }) }
  })
  observeEvent(input$Num2, ignoreNULL = FALSE, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    le <- length(file[[1]])
    if(le == 0){
      output$res10 <- renderText({
          paste()
       })
    } else {
    T_B <- input$Num2
    K <- input$Tapers2
    if(!is.na(T_B) & !is.na(K)){
      f_part <- c(1, floor(T_B/2 + 1))
      diff <- abs(diff(f_part)) / 2
      for(i in 1:T_B){
        temp_k <- i
        bw <- floor((temp_k+1)*(T_B/(T_B+1))) + 1
        if(bw < diff) {
          
        } else {
          break
        }
      }
      bw <- floor((K+1)*(T_B/(T_B+1))) + 1
      f_part <- c(1, floor(T_B/2 + 1))
      flo_mea <- floor(mean(f_part))
      max_tap <- ((flo_mea - 1) * ((T_B+1)/T_B)) - 1
      mmm <- abs(diff(f_part))
      output$res10 <- renderText({
        paste(h6("**Valid choices range from 1 to ", (i-1), "as we need to satisfy ", HTML(paste(tags$sup("floor(N/2)"))), 
                 "/", HTML(paste(tags$sub(2))), " - 1> floor((K+1)(", HTML(paste(tags$sup("N"))), 
                 "/", HTML(paste(tags$sub("N+1"))), "))"))
      })
    } 
    }
    
  })
  observeEvent(input$Num_Fxna, {
    curr_num <- as.numeric(input$Num_Fxna)
    output$Fxn_BB <- renderText({
      paste(h6("*Valid choices range from 1 to ", floor((curr_num/ 4 - 1)) , "as we need to satisfy 1", HTML("&le;"), 
               "K < floor(N/4 - 1)"))
    })
  })
  observeEvent(input$x_F1_AA, ignoreNULL = FALSE, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    le <- length(file[[1]])
    if(le == 0){
      
    } else if(is.na(plot.listbb()[[2]])){
      
    } else {
      get_num <- as.numeric(input$x_F1_AA)
      
      output$test12121_AA <- renderText({
        paste(h4(strong((paste(get_num)))))
      })
      output$Test_Fxna_Plot1 <- renderPlotly({
        a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(plot.listbb()[[1]][get_num,])), y=as.numeric(plot.listbb()[[1]][get_num,]))) + 
          xlab("Functional Domain") + ylab("") + ggtitle("Observed Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
        ggplotly(a)
      })
    }
    
  })
  output$res10 <- renderText({
    
  })
  plot.list2 <- eventReactive(input$go2, ignoreNULL = FALSE, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    dataf <- read.csv(file$datapath, header = input$header)
    dataf <- as.vector(dataf[[1]], mode = "numeric")
    ebaoutfu <- eba.search(X=dataf,N= as.numeric(input$Num2),K=as.numeric(input$Tapers2),std=input$TF2,alpha=as.numeric(input$Signi2))
    plot.x = ebaoutfu$mtspec$t
    plot.y = ebaoutfu$mtspec$f
    plot.z = t(ebaoutfu$mtspec$mtspec)
    plot.main = "Multitaper Spectrogram"
    plot.h = ebaoutfu$part.final[c(-1,-length(ebaoutfu$part.final))]
    plot.data = dataf
    plot.log <- ebaoutfu$log
    plot.pvals <- ebaoutfu$pvals
    plot.flat <- ebaoutfu$flat
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
         plot.main = plot.main, plot.h = plot.h, plot.data = plot.data, 
         plot.log = plot.log, plot.pvals = plot.pvals, plot.flat = plot.flat)})
  
  show("FxnPlotaDesc_AA")
  output$FxnPlotaDesc_AA <- renderText({
    paste(h4("Currently viewing timepoint "))
  })
  plot.listFxn2 <- eventReactive(input$go_Fxna,  {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    dataf <- read.csv(file$datapath, header = input$header)[,-1]
    dataf <- as.matrix(dataf)
    nrow = nrow(dataf); ncol = ncol(dataf)
    nb=15; #number of basis functions used to generate white noise
    R=dim(dataf)[2]; #number of points in functional domain
    Ts=dim(dataf)[1]; #length of time series
    seed=234; #seed for reproducibility
    B=floor(dim(dataf)[1]/as.numeric(input$Num_Fxna)); #number of time blocks
    N=as.numeric(input$Num_Fxna); #number of observations per time block
    bw=floor((as.numeric(input$Tapers_Fxna) + 1) / (as.numeric(input$Num_Fxna) + 1)); #bandwidth for multitaper spectral estimator
    K=as.numeric(input$Tapers_Fxna); #number of tapers for multitaper spectral estimator
    std=as.numeric(input$TF_Fxna); #standardize variance for points in functional domain (TRUE) or not (FALSE)
    freq=seq(from=0,by=1/N,length.out=floor(N/2)+1); #Fourier frequencies
    Rsel=as.numeric(input$Rsel_Fxna); #number of points in functional domain used for test statistics
    pse=fhat(dataf,N,K,Rsel,std);
    cmpnt="1-1"; #select component to view
    dimnames(pse) <- list(freq,apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-"),1:B);
    plot.x <- floor((1:B) * (Ts/B)); plot.y <- as.numeric(rownames(pse)); plot.z <- pse
    plot.cmp <- apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-")
    plot.main <- "Multitaper Autospectrum"; plot.data = dataf
    
    
    conf <- numeric(length(pse))
    dim(conf) <- dim(pse)
    dimnames(conf) <- dimnames(pse)
    for(k in 1:dim(pse)[3]){
      for(j in 1:dim(pse)[2]){
        first_col <- ((j - 1) %% (as.numeric(Rsel))) + 1 
        second_col <- ((j - 1) %/% (as.numeric(Rsel))) + 1
        cmpt_1 <- paste(first_col, "-", first_col, sep="")
        cmpt_2 <- paste(second_col, "-", second_col, sep="")
        for(i in 1:dim(pse)[1]){
          if(first_col == second_col){
            conf[i,j,k] <- Re(pse[i,j,k])
          } else {
            conf[i,j,k] <- Re((Mod(pse[i,j,k])**2) / (pse[i,cmpt_1,k] * pse[i,cmpt_2,k]))
          }
        }
      }
    }
    plot.z <- conf
    
    plot.data = dataf
    
    
    output$FxnbPlotDesc_file <- renderText({
      paste(h4("Currently viewing component "))
    })
    output$FxnPlot22Desc_file <- renderText({
      paste(h4("Currently viewing component "))
    })
    
    show("x_F1_AA")
    #show("q_F1")
    show("q11_F1_file")
    show("Plotly_Fxna_file")
    show("Plotly_Fxnb_file")
    show("plot1_FxnCheck_file")
    show("test12121_AA")
    show("Blank2_file")
    #show("Blank10100")
    #updateSliderInput(session, "x_F1_AA", min = 1, max=Ts, value=1, step=1)
    #updateSelectInput(session, "q_F1", choices=plot.cmp, selected = plot.cmp[1])
    updateSliderInput(session, "q11_F1_file", min=1, max=length(plot.cmp), value=1)
    updateSliderInput(session, "plot1_FxnCheck_file", min=1, max=length(plot.cmp), value=1)
    output$Blank2_file <- renderText({
      paste(h4(strong((paste("1-1")))))
    })
    output$Blank10100_file <- renderText({
      paste(h4(strong((paste("1-1")))))
    })
    output$test12121_AA <- renderText({
      paste(h4(strong((paste("1")))))
    })
    output$Plotly_Fxna_file <- renderPlotly({
      plot_ly(x = ~seq(from=0, to = 1, length.out=ncol(dataf)),
              y = ~seq(from=1, to=nrow(dataf)), 
              z = ~dataf) %>% add_surface() %>% layout(
        title = "3D Representation of Simulated Data",
        scene = list(
          xaxis = list(title="Functional Domain"), 
           yaxis = list(title = "Timepoint"), 
           zaxis = list(title="Value")
         )) %>% colorbar(title = "Value", len=1)
     })
    # show("Fxn_Row")
    plot.main_2 <- "Estimated Coherence"
    #list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
    #     plot.main = plot.main, plot.cmp = plot.cmp, plot.data = plot.data, plot.main_2 = plot.main_2)
    plot.main_3D <- "3D Representation of Autospectrogram"
    plot.main_3D.2 <- "3D Representation of Coherence"
    
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
         plot.main = plot.main, plot.cmp = plot.cmp, plot.data = plot.data, 
         nrow = nrow, ncol = ncol, plot.main_2 = plot.main_2, 
         plot.main_3D = plot.main_3D, plot.main_3D.2 = plot.main_3D.2)})
  
  output$Fxn_Plotb_file <- renderPlot({
    image.plot(x=plot.listFxn2()[[1]],y=plot.listFxn2()[[2]],z=suppressWarnings(t(Re(plot.listFxn2()[[3]][,"1-1",]))), 
               axes = TRUE, col = inferno(256), 
               main = plot.listFxn2()[[4]],xlab='Time',ylab='Hz',xaxs="i"); 
  })
  output$Plotly_Fxnb_file <- renderPlotly({
    plot_ly(y=~plot.listFxn2()[[1]], x=~plot.listFxn2()[[2]], z=~t(Re(plot.listFxn2()[[3]][,"1-1",])))  %>%layout(title=plot.listFxn2()[[10]],
                                                                                                                      scene = list(
                                                                                                                        xaxis = list(title='Frequency',range = c(0.5, 0)),
                                                                                                                        yaxis = list(title="Timepoint"),
                                                                                                                        zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
  })
  observeEvent(plot.listFxn2()[[6]],{ 
    output$Fxn_Plotb_file <- renderPlot({
      image.plot(x=plot.listFxn2()[[1]],y=plot.listFxn2()[[2]],z=suppressWarnings(t(Re(plot.listFxn2()[[3]][,"1-1",]))), 
                 axes = TRUE, col = inferno(256), 
                 main = plot.listFxn2()[[4]],xlab='Time',ylab='Hz',xaxs="i"); 
    })
    output$Plotly_Fxnb_file <- renderPlotly({
      plot_ly(y=~plot.listFxn2()[[1]], x=~plot.listFxn2()[[2]], z=~t(Re(plot.listFxn2()[[3]][,"1-1",])))  %>%layout(title=plot.listFxn2()[[10]],
                                                                                                                    scene = list(
                                                                                                                      xaxis = list(title='Frequency',range = c(0.5, 0)),
                                                                                                                      yaxis = list(title="Timepoint"),
                                                                                                                      zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
    })
    })
  observeEvent(input$plot1_FxnCheck_file, ignoreNULL = TRUE, {
    curr_num <- as.numeric(input$plot1_FxnCheck_file)
    if(is.na(plot.listFxn2()[[4]])){
      
    } else {
      curr_comp <- plot.listFxn2()[[5]][curr_num]
      output$Blank2_file <- renderText({
        paste(h4(strong((paste(curr_comp)))))
      })
      if(strsplit(curr_comp, "-")[[1]][1] == strsplit(curr_comp, "-")[[1]][2]){
        output$Fxn_Plotb_file <- renderPlot({
          image.plot(x=plot.listFxn2()[[1]],y=plot.listFxn2()[[2]],z=suppressWarnings(t(Re(plot.listFxn2()[[3]][,curr_comp,]))),
                     axes = TRUE, col = inferno(256),
                     main = plot.listFxn2()[[4]],xlab='Time',ylab='Hz',xaxs="i")
        })
        
      } else {
        output$Fxn_Plotb_file <- renderPlot({
          image.plot(x=plot.listFxn2()[[1]],y=plot.listFxn2()[[2]],z=suppressWarnings(t(Re(plot.listFxn2()[[3]][,curr_comp,]))),
                     axes = TRUE, col = inferno(256),
                     main = plot.listFxn2()[[9]],xlab='Time',ylab='Hz',xaxs="i");
        })
      }
    }
  })
  observeEvent(input$q11_F1_file, ignoreNULL = TRUE, {
    curr_num <- as.numeric(input$q11_F1_file)
    if(is.na(plot.listFxn2()[[4]])){
      
    } else {
      curr_comp <- plot.listFxn2()[[5]][curr_num]
      output$Blank10100_file <- renderText({
        paste(h4(strong((paste(curr_comp)))))
      })
      vals <- as.numeric(strsplit(curr_comp, "-")[[1]])
      if(vals[1] == vals[2]){
        output$Plotly_Fxnb_file <- renderPlotly({
          plot_ly(y=~plot.listFxn2()[[1]], x=~plot.listFxn2()[[2]], z=~t(Re(plot.listFxn2()[[3]][,curr_comp,])))  %>%layout(title=plot.listFxn2()[[10]],
                                                                                                                            scene = list(
                                                                                                                              xaxis = list(title='Frequency',range = c(0.5, 0)),
                                                                                                                              yaxis = list(title="Timepoint"),
                                                                                                                              zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
          
          
        })
      } else {
        output$Plotly_Fxnb_file <- renderPlotly({
          plot_ly(y=~plot.listFxn2()[[1]], x=~plot.listFxn2()[[2]], z=~t(Re(plot.listFxn2()[[3]][,curr_comp,])))  %>%layout(title=plot.listFxn2()[[11]],
                                                                                                                            scene = list(
                                                                                                                              xaxis = list(title='Frequency',range = c(0.5, 0)),
                                                                                                                              yaxis = list(title="Timepoint"),
                                                                                                                              zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
          
          
        })
      }
      
      
      
    }
    
  }) 
  output$Image_Plot2 <- renderPlot({
    par(mar=c(4,4,12,12))
    vp.top <- viewport(height=unit(0.4, "npc"), width=unit(0.8, "npc"),
                       just=c( "bottom"), y=0.6, x=0.475)
    plot.new()
    image.plot(x=plot.list2()[[1]], y=plot.list2()[[2]], z=plot.list2()[[3]], 
               axes = TRUE, col = inferno(256), 
               xlab='Time',ylab='Hz',xaxs="i", 
               bigplot = c(0.075, .675, .125, .925), smallplot = c(.7, .75, .125, .925));title(plot.list2()[[4]], line=0.75); 
    abline(h=plot.list2()[[5]], col = "skyblue", lwd=3); 
    
    vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                      just=c("left", "top"), y=0.625, x=0.7)
    len <- length(plot.list2()[[5]])
    vals <- plot.list2()[[5]]
    if(len == 0){
      str <- "(0, 0.5),"
    } else if (len == 1) {
      str <- paste("(0, ", round(vals, 3), "), [", round(vals, 3), ", 0.5),", sep="")
    } else {
      str <- paste("(0", sep="")
      for(i in 1:len){
        str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
      }
      str <- paste(str, ",", "0.5),", sep="")
    }
    spp <- strsplit(str, "),")[[1]]
    for(a in 1:length(spp)){
      spp[a] <- paste(spp[a], ")", sep="")
    }
    pp <- data.frame("Predicted Frequency Bands" = spp)
    colnames(pp) <- c("Predicted \n Frequency Bands")
    grid.table(pp, vp=vp.br, rows=NULL)
    
    vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                     just=c("left", "top"), y=0.8, x=0.7)
    grid.polygon(x=c(0.29, 0.29,0.71, 0.71), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
    jj <- grid.legend(c("Predicted Partition Points"), gp=gpar(lty=1, lwd=3, col=c("skyblue")), vp=vp.r, 
                      draw=TRUE)
  })
  
  observeEvent(input$x_F1, ignoreNULL = FALSE, {
    curr_row <- as.numeric(input$x_F1)
    output$test12121 <- renderText({
      paste(h4(strong(paste(curr_row))))
    })
    if(is.na(curr_row)){
      
    } else {
      output$Fxn_Plota <- renderPlotly({
        a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(plot.listF1()[[6]][curr_row,])), y=plot.listF1()[[6]][curr_row,])) + 
          xlab("Functional Domain") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits = c(0,1), expand=c(0,0))
        ggplotly(a)
      })
    }
    
  })
  #####
  observeEvent(input$q_F1, ignoreNULL = TRUE, {
    curr_comp <- as.character(input$q_F1)
    if(is.na(curr_comp)){
      
    } else {
      if(is.na(plot.listF1()[[4]])){
        
      } else {
        if(as.numeric(strsplit(curr_comp, "-")[[1]])[1] == as.numeric(strsplit(curr_comp, "-")[[1]][2])){
          title = "Multitaper Autospectrum"
        } else {
          title = "Estimated Coherence"
        }
      output$Fxn_Plotb <- renderPlot({
        image.plot(x=plot.listF1()[[1]],y=plot.listF1()[[2]],z=suppressWarnings(t(Re(plot.listF1()[[3]][,curr_comp,]))), 
                   axes = TRUE, col = inferno(256), 
                   main = title ,xlab='Time',ylab='Hz',xaxs="i"); 
      
      })
      }
    }
    
  })
  
  
  plot.list <- eventReactive(input$go, ignoreNULL = FALSE, {
    set.seed(823819)
    X = eba.simdata(T=as.numeric(input$Time))
    
    if (input$Simsetting == "W"){
      ebaout.wn <- eba.search(X=X$wn,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
      plot.x = ebaout.wn$mtspec$t
      plot.y = ebaout.wn$mtspec$f
      plot.z = t(ebaout.wn$mtspec$mtspec)
      plot.main = "Multitaper Spectrogram for White Noise Setting"
      plot.h = as.numeric(ebaout.wn$part.final[c(-1,-length(ebaout.wn$part.final))])
      plot.data = X$wn  
      plot.log = ebaout.wn$log
      plot.pvals = ebaout.wn$pvals
      plot.flat = ebaout.wn$flat
      
      
    } else if (input$Simsetting == "L") {
      ebaout.bL <- eba.search(X=X$bL,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
      plot.x = ebaout.bL$mtspec$t
      plot.y = ebaout.bL$mtspec$f
      plot.z = t(ebaout.bL$mtspec$mtspec)
      plot.main = "Multitaper Spectrogram for Linear Setting"
      plot.h = as.numeric(ebaout.bL$part.final[c(-1,-length(ebaout.bL$part.final))])
      plot.data = X$bL
      plot.log = ebaout.bL$log
      plot.pvals = ebaout.bL$pvals
      plot.flat = ebaout.bL$flat
      
    } else if (input$Simsetting == "S") {
      ebaout.bS <- eba.search(X=X$bS,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
      plot.x = ebaout.bS$mtspec$t
      plot.y = ebaout.bS$mtspec$f
      plot.z = t(ebaout.bS$mtspec$mtspec)
      plot.main = "Multitaper Spectrogram for Sinusoidal Setting"
      plot.h = as.numeric(ebaout.bS$part.final[c(-1,-length(ebaout.bS$part.final))])
      plot.data = X$bS
      plot.log = ebaout.bS$log
      plot.pvals = ebaout.bS$pvals
      plot.flat = ebaout.bS$flat
      
    }
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
         plot.main = plot.main, plot.h = plot.h, plot.data = plot.data, 
         plot.log = plot.log, plot.pvals = plot.pvals, plot.flat = plot.flat)
    
  });
  
  
  output$Image_Plot <- renderPlot({
    par(mar=c(4,4,12,12))
    vp.top <- viewport(height=unit(0.4, "npc"), width=unit(0.8, "npc"),
                       just=c( "bottom"), y=0.6, x=0.475)
    plot.new()
    image.plot(x=plot.list()[[1]], y=plot.list()[[2]], z=plot.list()[[3]], 
               axes = TRUE, col = inferno(256), 
               xlab='Time',ylab='Hz',xaxs="i", 
               bigplot = c(.1, .55, .1, .5), smallplot = c(.6, .65, .1, .5));title(plot.list()[[4]], line=0.75); 
    abline(h=plot.list()[[5]], col = "skyblue", lwd=3); abline(h=c(0.15, 0.35), col="lawngreen", lwd=3)
    
    vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                      just=c("left", "top"), y=0.5, x=0.6)
    act <- c("(0, 0.15)", "[0.15, 0.35)", "[0.35, 0.5)")
    len <- length(plot.list()[[5]])
    vals <- plot.list()[[5]]
    if(len == 0){
      str <- "(0, 0.5),"
    } else if (len == 1) {
      str <- paste("(0, ", round(vals, 3), "), [", round(vals, 3), ", 0.5),", sep="")
    } else {
      str <- paste("(0", sep="")
      for(i in 1:len){
        str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
      }
      str <- paste(str, ",", "0.5),", sep="")
    }
    spp <- strsplit(str, "),")[[1]]
    for(a in 1:length(spp)){
      spp[a] <- paste(spp[a], ")", sep="")
    }
    max_len <- max(length(act), length(spp))
    if(length(act) == length(spp)){
      
    } else if(length(act) > length(spp)){
      sp_l <- length(spp) + 1
      for(i in sp_l: length(act)){
        spp[i] <- ""
      }
    } else {
      ac_l <- length(act) + 1
      for(i in ac_l: length(spp)){
        act[i] <- ""
      }
    }
    pp <- data.frame("Actual Frequency Bands" = act, "Predicted Frequency Bands" = spp)
    colnames(pp) <- c("Actual \n Frequency Bands", "Predicted \n Frequency Bands")
    grid.table(pp, vp=vp.br, rows=NULL)
    
    vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                     just=c("left", "top"), y=0.65, x=0.6)
    grid.polygon(x=c(0.29, 0.29,0.71, 0.71), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
    jj <- grid.legend(c("Predicted Partition Points", "Actual Partition Points"), gp=gpar(lty=1, lwd=3, col=c("skyblue", "lawngreen")), vp=vp.r, 
                      draw=TRUE)
    
    print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.list()[[6]])), y= plot.list()[[6]])) + xlab("Time") +
             ylab("") + ggtitle("Simulated Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
             scale_x_continuous(limits=c(0,1), expand=c(0,0)), vp=vp.top)
    
    
     
  }); 
  output$summ_out_uni <- renderPlot({
    pvals <- round(plot.list()[[7]][,4], 5)
    pval.th <- round(plot.list()[[7]][,5], 5)
    Sig <- character(length(pvals))
    for(i in 1:length(Sig)){
      if(pvals[i] < pval.th[i]){
        Sig[i] <- "TRUE"
      } else {
        Sig[i] <- "FALSE"
      }
    }
    pp <- data.frame("Frequency" = round(plot.list()[[7]][,2], 3), "P-Value" = round(plot.list()[[7]][,4], 5), 
                     "P-Value\nThreshold" = round(plot.list()[[7]][,5], 5), "Significance" = as.character(Sig))
    colnames(pp) <- c("Frequency", "P-Value", "P-Value \n Threshold", "Significant")
    table <- tableGrob(pp, rows=NULL)
    title <- textGrob(expression(bold("Summary of Partition \n      Point Tests")))
    blank1 <- textGrob("")
    
    len <- length(plot.list()[[5]])
    vals <- plot.list()[[5]]
    if(len == 0){
      str <- "(0, 0.5),"
    } else if (len == 1) {
      str <- paste("(0, ", round(vals, 3), "), [", round(vals, 3), ", 0.5),", sep="")
    } else {
      str <- paste("(0", sep="")
      for(i in 1:len){
        str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
      }
      str <- paste(str, ",", "0.5),", sep="")
    }
    spp <- strsplit(str, "),")[[1]]
    for(a in 1:length(spp)){
      spp[a] <- paste(spp[a], ")", sep="")
    }
    pvals <- plot.list()[[9]][,2]
    Res <- character(length(pvals))
    Sig2 <- numeric(length(pvals))
    for(i in 1:length(Res)){
      if(pvals[i] < 0.05){
        Res[i] = "Segment has \n nonflat spectrum"
        Sig2[i] = "TRUE"
      } else {
        Res[i] = "Segment has \n flat spectrum"
        Sig2[i] = "FALSE"
      }
    }
    blank1 <- textGrob(""); blank2 <- textGrob("")
    new_tab <- data.frame("Frequency Bands" = spp, "P-Values" = round(as.numeric(pvals), 5),"Significant" = Sig2 ,"Results" = Res)
    colnames(new_tab) <- c("Frequency \n Bands", "P-Value", "Significant", "Results")
    test1 <- tableGrob(new_tab, rows = NULL); 
    title2 <- textGrob(expression(bold("Summary of Testing for Flat \n Spectrum in Each Segment")))
    #grid.arrange(test1)
    grid.arrange(title, table, title2, test1,blank2, heights = c(0.75,0.75,0.85,0.75, 1) ,nrow = 5)
  })
  output$summ_out_uni_file <- renderPlot({
    pvals <- round(plot.list2()[[7]][,4], 5)
    pval.th <- round(plot.list2()[[7]][,5], 5)
    Sig <- character(length(pvals))
    for(i in 1:length(Sig)){
      if(pvals[i] < pval.th[i]){
        Sig[i] <- "TRUE"
      } else {
        Sig[i] <- "FALSE"
      }
    }
    pp <- data.frame("Frequency" = round(plot.list2()[[7]][,2], 3), "P-Value" = round(plot.list2()[[7]][,4], 5), 
                     "P-Value\nThreshold" = round(plot.list2()[[7]][,5], 5), "Significance" = as.character(Sig))
    colnames(pp) <- c("Frequency", "P-Value", "P-Value \n Threshold", "Significant")
    table <- tableGrob(pp, rows=NULL)
    title <- textGrob(expression(bold("Summary of Partition \n      Point Tests")))
    blank1 <- textGrob("")
    
    len <- length(plot.list2()[[5]])
    vals <- plot.list2()[[5]]
    if(len == 0){
      str <- "(0, 0.5),"
    } else if (len == 1) {
      str <- paste("(0, ", round(vals, 3), "), [", round(vals, 3), ", 0.5),", sep="")
    } else {
      str <- paste("(0", sep="")
      for(i in 1:len){
        str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
      }
      str <- paste(str, ",", "0.5),", sep="")
    }
    spp <- strsplit(str, "),")[[1]]
    for(a in 1:length(spp)){
      spp[a] <- paste(spp[a], ")", sep="")
    }
    pvals <- plot.list2()[[9]][,2]
    Res <- character(length(pvals))
    Sig2 <- numeric(length(pvals))
    for(i in 1:length(Res)){
      if(pvals[i] < 0.05){
        Res[i] = "Segment has \n nonflat spectrum"
        Sig2[i] = "TRUE"
      } else {
        Res[i] = "Segment has \n flat spectrum"
        Sig2[i] = "FALSE"
      }
    }
    blank1 <- textGrob(""); blank2 <- textGrob("")
    new_tab <- data.frame("Frequency Bands" = spp, "P-Values" = round(as.numeric(pvals), 5),"Significant" = Sig2 ,"Results" = Res)
    colnames(new_tab) <- c("Frequency \n Bands", "P-Value", "Significant", "Results")
    test1 <- tableGrob(new_tab, rows = NULL); 
    title2 <- textGrob(expression(bold("Summary of Testing for Flat \n Spectrum in Each Segment")))
    #grid.arrange(test1)
    grid.arrange(title, table, title2, test1,blank2, heights = c(0.75,0.75,0.85,0.75, 1) ,nrow = 5)
  })
  output$summ_pval_uni <- renderPlot({
    #pvals <- plot.list()[[8]]
    
    #freqs <- names(pvals)
    #vals <- unname(pvals)
    
    ggplot() + geom_point(aes(x = as.numeric(plot.list()[[8]][,1]), y = as.numeric(plot.list()[[8]][,2]))) + xlim(c(0,0.5)) + ylim(c(0,1)) + 
      xlab("Frequency") + ylab("P-Value") + ggtitle("P-Values for Testing Partition Points") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
      geom_vline(xintercept = plot.list()[[5]], linetype = "dashed") + scale_x_continuous(expand=c(0,0), limits=c(0,0.5)) + scale_y_continuous(expand = c(0,0), limits=c(0,1))
  })
  output$summ_pval_uni_file <- renderPlot({
    
    
    ggplot() + geom_point(aes(x = as.numeric(plot.list2()[[8]][,1]), y = as.numeric(plot.list2()[[8]][,2]))) + xlim(c(0,0.5)) + ylim(c(0,1)) + 
      xlab("Frequency") + ylab("P-Value") + ggtitle("P-Values for Testing Partition Points") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
      geom_vline(xintercept = plot.list2()[[5]], linetype = "dashed") + scale_x_continuous(expand=c(0,0), limits=c(0,0.5)) + scale_y_continuous(expand = c(0,0), limits=c(0,1))
  })
  output$downloadData <- downloadHandler(
    filename = function(){
      paste("Simulated_Output_Results","pdf",sep = ".") 
    },
    content = function(file){
      pdf(file, paper = "USr", width = 1100, height=600, onefile = TRUE)
      par(mar=c(4,4,12,12))
      vp.top <- viewport(height=unit(0.4, "npc"), width=unit(0.8, "npc"),
                         just=c( "bottom"), y=0.6, x=0.475)
      #plot.new()
      image.plot(x=plot.list()[[1]], y=plot.list()[[2]], z=plot.list()[[3]], 
                 axes = TRUE, col = inferno(256), 
                 xlab='Time',ylab='Hz',xaxs="i", 
                 bigplot = c(.1, .55, .1, .5), smallplot = c(.6, .65, .1, .5));title(plot.list()[[4]], line=0.75); 
      abline(h=plot.list()[[5]], col = "skyblue", lwd=3); abline(h=c(0.15, 0.35), col="lawngreen", lwd=3)
      
      vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.43, "npc"), 
                        just=c("left", "top"), y=0.5, x=0.65)
      act <- c("(0, 0.15)", "[0.15, 0.35)", "[0.35, 0.5)")
      len <- length(plot.list()[[5]])
      vals <- plot.list()[[5]]
      if(len == 0){
        str <- "(0, 0.5),"
      } else if (len == 1) {
        str <- paste("(0, ", round(vals, 3), "),[", round(vals, 3), ", 0.5),", sep="")
      } else {
        str <- paste("(0", sep="")
        for(i in 1:len){
          str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
        }
        str <- paste(str, ",", "0.5),", sep="")
      }
      spp <- strsplit(str, "),")[[1]]
      for(a in 1:length(spp)){
        spp[a] <- paste(spp[a], ")", sep="")
      }
      max_len <- max(length(act), length(spp))
      if(length(act) == length(spp)){
        
      } else if(length(act) > length(spp)){
        sp_l <- length(spp) + 1
        for(i in sp_l: length(act)){
          spp[i] <- ""
        }
      } else {
        ac_l <- length(act) + 1
        for(i in ac_l: length(spp)){
          act[i] <- ""
        }
      }
      pp <- data.frame("Actual Frequency Bands" = act, "Predicted Frequency Bands" = spp)
      colnames(pp) <- c("Actual \n Frequency Bands", "Predicted \n Frequency Bands")
      
        vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.43, "npc"), 
                          just=c("left", "top"), y=0.5, x=0.63)
        grid.table(pp, vp=vp.br, rows=NULL)
        vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.5, "npc"), 
                         just=c("left", "top"), y=0.65, x=0.58)
        grid.polygon(x=c(0.29, 0.29,0.71, 0.71), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
        
       jj <- grid.legend(c("Predicted Partition Points", "Actual Partition Points"), gp=gpar(lty=1, lwd=3, col=c("skyblue", "lawngreen")), vp=vp.r, 
                        draw=TRUE)
      
      print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.list()[[6]])), y= plot.list()[[6]])) + xlab("Time") +
               ylab("") + ggtitle("Simulated Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
               scale_x_continuous(limits=c(0,1), expand=c(0,0)), vp=vp.top)
      #plot.new()
      vp.r <- viewport(height=unit(1, "npc"), width=unit(0.5, "npc"), 
                      y=0.5, x=0.75)
      vp.l <- viewport(height=unit(1, "npc"), width=unit(0.5, "npc"), 
                       y=0.5, x=0.25)
      pvals <- round(plot.list()[[7]][,4], 5)
      pval.th <- round(plot.list()[[7]][,5], 5)
      Sig <- character(length(pvals))
      for(i in 1:length(Sig)){
        if(pvals[i] < pval.th[i]){
          Sig[i] <- "TRUE"
        } else {
          Sig[i] <- "FALSE"
        }
      }
      pp <- data.frame("Frequency" = round(plot.list()[[7]][,2], 3), "P-Value" = round(plot.list()[[7]][,4], 5), 
                       "P-Value\nThreshold" = round(plot.list()[[7]][,5], 5), "Significance" = as.character(Sig))
      colnames(pp) <- c("Frequency", "P-Value", "P-Value \n Threshold", "Significant")
      table <- tableGrob(pp, rows=NULL)
      title <- textGrob(expression(bold("Summary of Partition \n      Point Tests")))
      blank1 <- textGrob("")
      
      len <- length(plot.list()[[5]])
      vals <- plot.list()[[5]]
      if(len == 0){
        str <- "(0, 0.5),"
      } else if (len == 1) {
        str <- paste("(0, ", round(vals, 3), "), [", round(vals, 3), ", 0.5),", sep="")
      } else {
        str <- paste("(0", sep="")
        for(i in 1:len){
          str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
        }
        str <- paste(str, ",", "0.5),", sep="")
      }
      spp <- strsplit(str, "),")[[1]]
      for(a in 1:length(spp)){
        spp[a] <- paste(spp[a], ")", sep="")
      }
      pvals <- plot.list()[[9]][,2]
      Res <- character(length(pvals))
      Sig2 <- numeric(length(pvals))
      for(i in 1:length(Res)){
        if(pvals[i] < 0.05){
          Res[i] = "Segment has \n nonflat spectrum"
          Sig2[i] = "TRUE"
        } else {
          Res[i] = "Segment has \n flat spectrum"
          Sig2[i] = "FALSE"
        }
      }
      blank1 <- textGrob(""); blank2 <- textGrob("")
      new_tab <- data.frame("Frequency Bands" = spp, "P-Values" = round(as.numeric(pvals), 5),"Significant" = Sig2 ,"Results" = Res)
      colnames(new_tab) <- c("Frequency \n Bands", "P-Value", "Significant", "Results")
      test1 <- tableGrob(new_tab, rows = NULL); 
      title2 <- textGrob(expression(bold("Summary of Testing for Flat \n Spectrum in Each Segment")))
      #grid.arrange(test1)
      grid.arrange(title, table, title2, test1,blank2, heights = c(0.75,0.75,0.85,0.75, 1) ,nrow = 5, vp = vp.l)
      print(ggplot() + geom_point(aes(x = as.numeric(plot.list()[[8]][,1]), y = as.numeric(plot.list()[[8]][,2]))) + xlim(c(0,0.5)) + ylim(c(0,1)) + 
              xlab("Frequency") + ylab("P-Value") + ggtitle("P-Values for Testing Partition Points") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
              geom_vline(xintercept = plot.list()[[5]], linetype = "dashed") , vp=vp.r)
      
      dev.off()
    }
  )
  output$downloadData1 <- downloadHandler(
    filename = function(){
      paste("Observed_Output_Results","pdf",sep = ".") 
    },
    content = function(file){
      pdf(file, paper = "USr", width = 1100, height=600, onefile = TRUE)
      par(mar=c(4,4,12,12))
      vp.top <- viewport(height=unit(0.4, "npc"), width=unit(0.8, "npc"),
                         just=c( "bottom"), y=0.6, x=0.475)
      image.plot(x=plot.list2()[[1]], y=plot.list2()[[2]], z=plot.list2()[[3]], 
                 axes = TRUE, col = inferno(256), 
                 xlab='Time',ylab='Hz',xaxs="i", 
                 bigplot = c(.125, .575, .125, .525), smallplot = c(.6, .65, .1, .5));title(plot.list2()[[4]], line=0.75); 
      abline(h=plot.list2()[[5]], col = "skyblue", lwd=3); 
      
      len <- length(plot.list2()[[5]])
      vals <- plot.list2()[[5]]
      if(len == 0){
        str <- "(0, 0.5),"
      } else if (len == 1) {
        str <- paste("(0, ", round(vals, 3), "), [", round(vals, 3), ", 0.5),", sep="")
      } else {
        str <- paste("(0", sep="")
        for(i in 1:len){
          str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
        }
        str <- paste(str, ",", "0.5),", sep="")
      }
      spp <- strsplit(str, "),")[[1]]
      for(a in 1:length(spp)){
        spp[a] <- paste(spp[a], ")", sep="")
      }
      pp <- data.frame("Predicted Frequency Bands" = spp)
      colnames(pp) <- c("Predicted \n Frequency Bands")
      
        vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.43, "npc"), 
                          just=c("left", "top"), y=0.5, x=0.63)
        grid.table(pp, vp=vp.br, rows=NULL)
        vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.5, "npc"), 
                         just=c("left", "top"), y=0.65, x=0.58)
        grid.polygon(x=c(0.29, 0.29,0.71, 0.71), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
        
    
      
      jj <- grid.legend(c("Predicted Partition Points"), gp=gpar(lty=1, lwd=3, col=c("skyblue")), vp=vp.r, 
                        draw=TRUE)
      
      print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.list2()[[6]])), y= plot.list2()[[6]])) + xlab("Time") +
               ylab("") + ggtitle("Observed Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
               scale_x_continuous(limits=c(0,1), expand=c(0,0)), vp=vp.top)
      vp.r <- viewport(height=unit(1, "npc"), width = unit(0.5, "npc"), 
                       x=0.75, y=0.5)
      vp.l <- viewport(height=unit(1, "npc"), width = unit(0.5, "npc"), 
                       x=0.25, y=0.5)
      pvals <- round(plot.list2()[[7]][,4], 5)
      pval.th <- round(plot.list2()[[7]][,5], 5)
      Sig <- character(length(pvals))
      for(i in 1:length(Sig)){
        if(pvals[i] < pval.th[i]){
          Sig[i] <- "TRUE"
        } else {
          Sig[i] <- "FALSE"
        }
      }
      pp <- data.frame("Frequency" = round(plot.list2()[[7]][,2], 3), "P-Value" = round(plot.list2()[[7]][,4], 5), 
                       "P-Value\nThreshold" = round(plot.list2()[[7]][,5], 5), "Significance" = as.character(Sig))
     
      colnames(pp) <- c("Frequency", "P-Value", "P-Value \n Threshold", "Significant")
      table <- tableGrob(pp, rows=NULL)
      title <- textGrob(expression(bold("Summary of Partition \n      Point Tests")))
      blank1 <- textGrob("")
      #grid.arrange(title, table, heights = c(2, 1, 5))
      
      len <- length(plot.list2()[[5]])
      vals <- plot.list2()[[5]]
      if(len == 0){
        str <- "(0, 0.5),"
      } else if (len == 1) {
        str <- paste("(0, ", round(vals, 3), "), [", round(vals, 3), ", 0.5),", sep="")
      } else {
        str <- paste("(0", sep="")
        for(i in 1:len){
          str <- paste(str, ", ",round(vals[i], 3),"),[", round(vals[i], 3), sep="")
        }
        str <- paste(str, ",", "0.5),", sep="")
      }
      spp <- strsplit(str, "),")[[1]]
      for(a in 1:length(spp)){
        spp[a] <- paste(spp[a], ")", sep="")
      }
      pvals <- plot.list2()[[9]][,2]
      Res <- character(length(pvals))
      Sig2 <- numeric(length(pvals))
      for(i in 1:length(Res)){
        if(pvals[i] < 0.05){
          Res[i] = "Segment has \n nonflat spectrum"
          Sig2[i] = "TRUE"
        } else {
          Res[i] = "Segment has \n flat spectrum"
          Sig2[i] = "FALSE"
        }
      }
      blank1 <- textGrob(""); blank2 <- textGrob("")
      new_tab <- data.frame("Frequency Bands" = spp, "P-Values" = round(as.numeric(pvals), 5),"Significant" = Sig2 ,"Results" = Res)
      colnames(new_tab) <- c("Frequency \n Bands", "P-Value", "Significant", "Results")
      test1 <- tableGrob(new_tab, rows = NULL); 
      title2 <- textGrob(expression(bold("Summary of Testing for Flat \n Spectrum in Each Segment")), gp=gpar(fontface = "bold"))
      grid.arrange(title, table, title2, test1,blank2, heights = c(0.75,0.75,0.85,0.75, 1) ,nrow = 5,
                   vp = vp.l)
      print (ggplot() + geom_point(aes(x = as.numeric(plot.list2()[[8]][,1]), y = as.numeric(plot.list2()[[8]][,2]))) + xlim(c(0,0.5)) + ylim(c(0,1)) + 
               xlab("Frequency") + ylab("P-Value") + scale_x_continuous(expand=c(0,0), limits=c(0,0.5)) + scale_y_continuous(expand = c(0,0), limits=c(0,1)) + ggtitle("P-Values for Testing Partition Points") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
               geom_vline(xintercept = plot.list2()[[5]], linetype = "dashed"), vp = vp.r)
      dev.off()
    }
  )
  
  
}

shinyApp(ui = ui, server = server)