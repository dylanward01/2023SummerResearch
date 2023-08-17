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
library(fEBAcpp)
library(shinyjs)
library(fda)
library(plotly)

source("EBA_functions.R")
source("fEBA_Rfns.R")
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
  tags$p(tags$strong("Authors: Dylan Ward, Kevin Gunalan Antony Michael Raj")), (tags$em(tags$u("Under the guidance of Professor Scott Bruce"))),
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
                                   choices=as.numeric(31), selected=31),
                       selectInput(inputId="Tapers", label="Choose number of tapers to use in multitaper spectral estimator** (K)", 
                                  choices=as.numeric(5), selected=5),
                      selectInput(inputId = "Signi", label="Choose significance level", 
                                  choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05),
                       radioButtons(inputId = "TF", label = "Standardize", 
                                    c("True" = TRUE, "False" = FALSE), selected = FALSE),
                       actionButton("go", label = "Run"),
                      htmlOutput("Res"),
                      htmlOutput("Res1"),
                      ), 
                      conditionalPanel(condition = "input.type== 'Functional'", 
                      selectInput("SimF1", "Simulation Setting",
                                  c("White Noise" = "W",
                                  "Linear" = "L",
                                  "Sinusoidal" = "S"), selected="S"),
                      selectInput(inputId = "TsF1", label="Choose total length of time series (T)", 
                                  choices = c(200, 500, 1000, 2000, 5000), selected=2000), 
                      selectInput(inputId = "RF1", label = "Choose number of points in functional domain (R)", 
                                  choices=seq(from=5, to=50, by=5), selected=5), 
                      selectInput(inputId = "NF1", label = "Choose number of observations per approximately stationary block* (N)",
                                  choices = 30, selected = 30), 
                      selectInput(inputId = "KF1", label = "Choose number of tapers to use in multitaper spectral estimator** (K)", 
                                  choices = NULL, selected = NULL),
                      selectInput(inputId = "RselF1", label = "Choose number of points in the functional domain to use in computing test statistics*** (Rsel)", 
                                  choices = NULL, selected = NULL), 
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
                       hidden(radioButtons(inputId = "Data_Checker", label = "Choose", choices = c("Functional", "Multivariate"), selected="Functional")),
                       htmlOutput("T_len"),
                       numericInput(inputId = "Num2", label = "Choose number of observations per approximately stationary block* (N)", 
                                   value = NULL, step = 1),
                       
                       numericInput(inputId = "Tapers2", label = "Choose number of tapers to use in multitaper spectral estimator** (K)", 
                                   value = NULL, step = 1),
                       
                       selectInput(inputId = "Signi2", label="Choose significance level", 
                                   choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05),
                       
                       radioButtons(inputId = "TF2", label = "Standardize", 
                                    c("True" = TRUE, "False" = FALSE), selected = FALSE),
                       actionButton("go2", label = "Run"), 
                       htmlOutput("res9"),
                       htmlOutput("res10"),
                       ), 
     
                      
    ),
    mainPanel(
      tabsetPanel(type = "tabs", id = 'tabselected', selected = 1,
                  tabPanel("Simulation Setting", value = 1),
                  tabPanel("File Upload", value = 2)),
      conditionalPanel(condition = "input.tabselected==1",
                       conditionalPanel(condition = "input.type == 'Univariate'",
                       plotOutput("Image_Plot", height=1000),
                       radioButtons(inputId = "downloadType", label = "Select download type", choices = list("png","pdf")),
                       downloadButton('downloadData','Download the plot') ),
                       conditionalPanel(condition = "input.type == 'Functional'",
                       plotlyOutput("Fxn_Plota", height=400),
                       splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                         cellWidths = c( "20%", "20%", "30%"),htmlOutput("FxnPlotaDesc"), 
                                   hidden(selectInput(inputId = "x_F1", label=NULL, choices = NULL, selected = NULL)), 
                                   htmlOutput("test12121")),
                     
                       plotOutput("Fxn_Plotb", height=600), 
                       splitLayout(cellWidths = c("20%", "25%", "25%", "30%"), htmlOutput("Blank101"),
                                   htmlOutput("FxnbPlotDesc"), 
                                   hidden(selectInput(inputId = "q_F1", label=NULL, choices=NULL, selected = NULL)),
                                   htmlOutput("Blank2")), 
                       splitLayout(cellWidths = c("80%", "20%"),
                                   (plotlyOutput("Plotly_Fxna", height=600)), 
                                   verticalLayout(plotOutput("Blank12121", height = 150), 
                                   hidden(radioButtons("Fxn_Row", label="View a Single Row?", choices=c("Yes", "No" ), 
                                                       selected = "No"))
                                   )
                                   ),
                       hidden(plotOutput("Plotly_Fxna.5", height = 400)),
                       splitLayout(tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                   cellWidths = c( "20%", "20%", "30%"),htmlOutput("FxnPlot11Desc"), 
                                   hidden(selectInput(inputId = "x11_F1", label=NULL, choices = NULL, selected = NULL)), 
                                   htmlOutput("test00")),
                       plotlyOutput("Plotly_Fxnb", height=600), 
                       splitLayout(cellWidths = c("20%", "25%", "25%", "30%"), htmlOutput("Blank10100"),
                                   htmlOutput("FxnPlot22Desc"), 
                                   hidden(selectInput(inputId = "q11_F1", label=NULL, choices=NULL, selected = NULL)),
                                   htmlOutput("Blank2023"))
      )),
      conditionalPanel(condition = "input.tabselected==2",
                       plotOutput("Image_Plota", height=400),
                       plotOutput("Image_Plot2", height=600), 
                       downloadButton('downloadData1','Download the plot')
      ), 
      )
  ))

server <- function(input,output, session) {
  
  output$BlankSpace <- renderText({
    paste("\n")
  })
  observe({
    TF1 <- as.numeric(input$TsF1)
    new_vals <- floor(seq(from=30, to=TF1/2, length.out = 10))
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
    R_val <- as.numeric(input$RF1)
    new_vals <- seq(from=1, to=R_val, by=1)
    updateSelectInput(session, "RselF1", choices = new_vals, selected = new_vals[length(new_vals)-1])
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
    paste(h6("*Choices must satisfy 30 ", HTML("&le;"), "N", HTML("&le;"), "T/2"))
  })
  output$F1_2 <- renderText({
    paste(h6("**Choices must satisfy 1 ", HTML("&le;"), "K", "<", "floor(N/4 - 1)"))
  })
  output$F1_3 <- renderText({
    paste(h6("***Choices must satisfy 1 ", HTML("&le;"), "Rsel", HTML("&le;"), "R"))
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
      output$Fxn_Plota <- renderPlotly({
        a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(X[1,])), y=X[1,])) + 
          xlab("Time") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
        ggplotly(a)
      })
      pse=fhat(X,N,K,Rsel,std);
       
       cmpnt="1-1"; #select component to view
       dimnames(pse) <- list(freq,apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-"),1:B);
       plot.x <- floor((1:B) * (Ts/B)); plot.y <- as.numeric(rownames(pse)); plot.z <- pse
       plot.cmp <- apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-")
       plot.main <- "Multitaper Autospectrum"; plot.data = X
     
      
    } else if (input$SimF1 == "L") {
      X=f3bL.sim(nb=nb,gsz=R,Ts=Ts,seed=seed);
      output$Fxn_Plota <- renderPlotly({
        a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(X[1,])), y=X[1,])) + 
          xlab("Time") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
        ggplotly(a)
      })
      pse=fhat(X,N,K,Rsel,std);
       cmpnt="1-1"; #select component to view
       dimnames(pse) <- list(freq,apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-"),1:B);
       plot.x <- floor((1:B) * (Ts/B)); plot.y <-as.numeric( rownames(pse)); plot.z <- pse
       plot.cmp <- apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-")
       plot.main <- "Multitaper Autospectrum"; plot.data = X
       
      
    } else if (input$SimF1 == "S") {
      X=f3bS.sim(nb=nb,gsz=R,Ts=Ts,seed=seed);
      output$Fxn_Plota <- renderPlotly({
        a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(X[1,])), y=X[1,])) + 
          xlab("Time") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
        ggplotly(a)
      })
       pse=fhat(X,N,K,Rsel,std);
       cmpnt="1-1"; #select component to view
       dimnames(pse) <- list(freq,apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-"),1:B);
       plot.x <- floor((1:B) * (Ts/B)); plot.y <- as.numeric(rownames(pse)); plot.z <- pse
       plot.cmp <- apply(expand.grid(1:Rsel,1:Rsel),1,paste,collapse = "-")
       plot.main <- "Multitaper Autospectrum"; plot.data = X
       
    }
    
    
    output$FxnPlotaDesc <- renderText({
      paste(h4("Currently viewing row "))
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
    updateSelectInput(session, "x_F1", choices = seq(from=1, to = dim(X)[1], by=1), selected=1)
    updateSelectInput(session, "q_F1", choices=plot.cmp, selected = plot.cmp[1])
    updateSelectInput(session, "q11_F1", choices=plot.cmp, selected = plot.cmp[1])
    output$Plotly_Fxna <- renderPlotly({
      plot_ly(z = ~X) %>% add_surface() %>% layout(scene = list(
        xaxis = list(title="Functional Domain",
                     range = c(0, as.numeric(input$RF1) + 0.5),
                     ticktext = (seq(1, as.numeric(input$RF1))), 
                     tickvals = (seq(0, as.numeric(input$RF1) - 1)), 
                     tickmode = "array"), 
        yaxis = list(title = "Row"), 
        zaxis = list(title="Value")
      )) %>% colorbar(title = "Value", len=1)
    })
    show("Fxn_Row")
    
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
        plot.main = plot.main, plot.cmp = plot.cmp, plot.data = plot.data)
    
  });
  
  output$Fxn_Plotb <- renderPlot({
    image.plot(x=plot.listF1()[[1]],y=plot.listF1()[[2]],z=suppressWarnings(t(Re(plot.listF1()[[3]][,"1-1",]))), 
               axes = TRUE, col = inferno(256), 
               main = plot.listF1()[[4]],xlab='Time',ylab='Hz',xaxs="i"); 
  })
  
  output$Plotly_Fxnb <- renderPlotly({
    plot_ly(y=~plot.listF1()[[1]], x=~plot.listF1()[[2]], z=~t(Re(plot.listF1()[[3]][,"1-1",])))  %>%layout(scene = list( 
           xaxis = list(title='Frequency',range = c(0, 0.5)), 
           yaxis = list(title="Time"), 
           zaxis = list(title="Value"))) %>% add_surface() %>% colorbar(title="Value", len=1)
  })
  observeEvent(input$Fxn_Row, ignoreNULL = TRUE, {
    if(input$Fxn_Row == 'No'){
      hide("Plotly_Fxna.5")
    } else {
      show("Plotly_Fxna.5")
      output$Plotly_Fxna.5 <- renderPlot({
        ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(plot.listF1()[[6]][1,])), y=plot.listF1()[[6]][1,])) + 
          xlab("Time") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
      })
      show("x11_F1")
      output$FxnPlot11Desc <- renderText({
        paste(h4("Currently viewing row "))
      })
      updateSelectInput(session, "x11_F1", choices = seq(from=1, to = dim(plot.listF1()[[6]])[1], by=1), selected=1)
    }
  })
 
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
  observeEvent(input$file_csv, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    dataf <- read.csv(file$datapath, header = input$header)
    dims <- dim(dataf)[2]
    if(dims == 1){
      hide("NotUniVar")
      show("UniVarDis")
      output$UniVarDis <- renderText({
        paste(strong("This is a Univariate Time Series"))
      })
      hide("Data_Checker")
      show("T_len")
      show("Num2")
      show("Tapers2")
      show("Signi2")
      show("TF2")
      show("go2")
      show("res9")
      show("res10")
      show("Image_Plota")
      show("Image_Plotb")
      show("downloadData1")
    } else {
      hide("UniVarDis")
      show("NotUniVar")
      output$NotUniVar <- renderText({
        paste(strong("Choose what type of time series this is"))
      })
      show("Data_Checker")
      hide("T_len")
      hide("Num2")
      hide("Tapers2")
      hide("Signi2")
      hide("TF2")
      hide("go2")
      hide("res9")
      hide("res10")
      hide("Image_Plota")
      hide("Image_Plotb")
      hide("downloadData1")
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
      mmm <- abs(diff(f.part))
      output$res10 <- renderText({
        paste(h6("**Valid choices range from 1 to ", (i-1), "as we need to satisfy ", HTML(paste(tags$sup("floor(N/2)"))), 
                 "/", HTML(paste(tags$sub(2))), " - 1> floor((K+1)(", HTML(paste(tags$sup("N"))), 
                 "/", HTML(paste(tags$sub("N+1"))), "))"))
      })
    } 
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
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
         plot.main = plot.main, plot.h = plot.h, plot.data = plot.data)})
  
  
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
    if(is.na(curr_row)){
      
    } else {
      output$Fxn_Plota <- renderPlotly({
        a <- ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(plot.listF1()[[6]][curr_row,])), y=plot.listF1()[[6]][curr_row,])) + 
          xlab("Time") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
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
      output$Fxn_Plotb <- renderPlot({
        image.plot(x=plot.listF1()[[1]],y=plot.listF1()[[2]],z=suppressWarnings(t(Re(plot.listF1()[[3]][,curr_comp,]))), 
                   axes = TRUE, col = inferno(256), 
                   main = plot.listF1()[[4]],xlab='Time',ylab='Hz',xaxs="i"); 
      
      })
      }
    }
    
  })
  
  observeEvent(input$x11_F1, ignoreNULL = FALSE, {
    curr_row <- as.numeric(input$x11_F1)
    if(is.na(curr_row)){
      
    } else {
      output$Plotly_Fxna.5 <- renderPlot({
        ggplot() + geom_line(aes(x=seq(from=0, to=1, length.out=length(plot.listF1()[[6]][curr_row,])), y=plot.listF1()[[6]][curr_row,])) + 
          xlab("Time") + ylab("") + ggtitle("Simulated Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
          scale_x_continuous(limits=c(0,1), expand=c(0,0))
      })
    }
    
  })
  
  observeEvent(input$q11_F1, ignoreNULL = TRUE, {
    curr_comp <- as.character(input$q11_F1)
    if(is.na(curr_comp)){
      
    } else {
      if(is.na(plot.listF1()[[4]])){
        
      } else {
        output$Plotly_Fxnb <- renderPlotly({
          plot_ly(y=~plot.listF1()[[1]], x=~plot.listF1()[[2]], z=~t(Re(plot.listF1()[[3]][,curr_comp,])))  %>%layout(scene = list( 
            xaxis = list(title='Frequency',range = c(0, 0.5)), 
            yaxis = list(title="Time"), 
            zaxis = list(title="Value"), 
            camera = list(eye = list(x=0, y=0, z=0)))) %>% add_surface() %>% colorbar(title="Value", len=1)
          
        })
      }
    }
    
  })
  plot.list <- eventReactive(input$go, ignoreNULL = FALSE, {
    X = eba.simdata(T=as.numeric(input$Time))
    
    if (input$Simsetting == "W"){
      ebaout.wn <- eba.search(X=X$wn,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
      plot.x = ebaout.wn$mtspec$t
      plot.y = ebaout.wn$mtspec$f
      plot.z = t(ebaout.wn$mtspec$mtspec)
      plot.main = "Multitaper Spectrogram for White Noise Setting"
      plot.h = as.numeric(ebaout.wn$part.final[c(-1,-length(ebaout.wn$part.final))])
      plot.data = X$wn  
      
      
      
    } else if (input$Simsetting == "L") {
      ebaout.bL <- eba.search(X=X$bL,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
      plot.x = ebaout.bL$mtspec$t
      plot.y = ebaout.bL$mtspec$f
      plot.z = t(ebaout.bL$mtspec$mtspec)
      plot.main = "Multitaper Spectrogram for Linear Setting"
      plot.h = as.numeric(ebaout.bL$part.final[c(-1,-length(ebaout.bL$part.final))])
      plot.data = X$bL
      
      
    } else if (input$Simsetting == "S") {
      ebaout.bS <- eba.search(X=X$bS,N= as.numeric(input$Num),K=as.numeric(input$Tapers),std=input$TF,alpha=as.numeric(input$Signi))
      plot.x = ebaout.bS$mtspec$t
      plot.y = ebaout.bS$mtspec$f
      plot.z = t(ebaout.bS$mtspec$mtspec)
      plot.main = "Multitaper Spectrogram for Sinusoidal Setting"
      plot.h = as.numeric(ebaout.bS$part.final[c(-1,-length(ebaout.bS$part.final))])
      plot.data = X$bS
      
    }
    
    list(plot.x = plot.x, plot.y = plot.y, plot.z = plot.z, 
         plot.main = plot.main, plot.h = plot.h, plot.data = plot.data)
    
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
  
  output$downloadData <- downloadHandler(
    filename = function(){
      paste("Simulated_Output",input$downloadType,sep = ".") 
    },
    content = function(file){
      if(input$downloadType == "png") png(file, width = 1000, height=600)
      else pdf(file, paper = "USr", width = 1100, height=600, onefile = FALSE)
      par(mar=c(4,4,12,12))
      vp.top <- viewport(height=unit(0.4, "npc"), width=unit(0.8, "npc"),
                         just=c( "bottom"), y=0.6, x=0.475)
      plot.new()
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
      #colnames(pp) <- stringr::str_replace_all(colnames(df), "\\n", "<br>")
      grid.table(pp, vp=vp.br, rows=NULL)
      
      vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.43, "npc"), 
                       just=c("left", "top"), y=0.65, x=0.65)
      grid.polygon(x=c(0.29, 0.29,0.71, 0.71), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
      jj <- grid.legend(c("Predicted Partition Points", "Actual Partition Points"), gp=gpar(lty=1, lwd=3, col=c("skyblue", "lawngreen")), vp=vp.r, 
                        draw=TRUE)
      
      print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.list()[[6]])), y= plot.list()[[6]])) + xlab("Time") +
               ylab("") + ggtitle("Simulated Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
               scale_x_continuous(limits=c(0,1), expand=c(0,0)), vp=vp.top)
      
      
      dev.off()
    }
  )
  output$downloadData1 <- downloadHandler(
    filename = function(){
      paste("Observed_Output",input$downloadType,sep = ".") 
    },
    content = function(file){
      if(input$downloadType == "png") png(file, width = 1000, height=600)
      else pdf(file, paper = "USr", width = 1100, height=600, onefile = FALSE)
      par(mar=c(4,4,12,12))
      vp.top <- viewport(height=unit(0.4, "npc"), width=unit(0.8, "npc"),
                         just=c( "bottom"), y=0.6, x=0.475)
      plot.new()
      image.plot(x=plot.list2()[[1]], y=plot.list2()[[2]], z=plot.list2()[[3]], 
                 axes = TRUE, col = inferno(256), 
                 xlab='Time',ylab='Hz',xaxs="i", 
                 bigplot = c(.125, .575, .125, .525), smallplot = c(.6, .65, .1, .5));title(plot.list2()[[4]], line=0.75); 
      abline(h=plot.list2()[[5]], col = "skyblue", lwd=3); 
      
      vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                        just=c("left", "top"), y=0.5, x=0.6)
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
                       just=c("left", "top"), y=0.65, x=0.6)
      grid.polygon(x=c(0.29, 0.29,0.71, 0.71), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
      jj <- grid.legend(c("Predicted Partition Points"), gp=gpar(lty=1, lwd=3, col=c("skyblue")), vp=vp.r, 
                        draw=TRUE)
      
      print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.list2()[[6]])), y= plot.list2()[[6]])) + xlab("Time") +
               ylab("") + ggtitle("Observed Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
               scale_x_continuous(limits=c(0,1), expand=c(0,0)), vp=vp.top)
      
      
      dev.off()
    }
  )
  
  
}

shinyApp(ui = ui, server = server)