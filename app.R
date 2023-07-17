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

source("EBA_functions.R")
df <- read.csv("xtsample.csv")
df <- df[[1]]

ui <- fluidPage(
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
                       radioButtons("type", "Method Type", 
                                    choices=c("Type A","Type B"), 
                                    selected="Type A"), 
                       conditionalPanel(condition = "input.type == 'Type A'", 
                       selectInput("Simsetting", "Simulation Setting",
                                   c("White Noise" = "W",
                                     "Linear" = "L",
                                     "Sinusoidal" = "S"), selected="S"),
                       selectInput(inputId="Time", label = "Choose total length of time series", 
                                   choices = as.numeric(c(500,1000,5000,10000, 20000, 50000)), selected=1000), 
                       #sliderInput(inputId = "Time", label = "Choose total length of time series", 
                        #           value = 50000, min = 1000, max = 50000, step = 1000),
                       
                       selectInput(inputId="Num", label = "Choose number of observations per approximately stationary block*", 
                                   choices=as.numeric(31), selected=31),
                       #sliderInput(inputId = "Num", label = "Choose number of observations per approximately stationary block", 
                        #           value = 500, min = 50, max = 1000, step = 50),
                       
                       #sliderInput(inputId = "Tapers", label = "Choose number of tapers to use in multitaper spectral estimator", 
                      #             value = 15, min = 0, max = 50, step = 1),
                       
                      selectInput(inputId="Tapers", label="Choose number of tapers to use in multitaper spectral estimator**", 
                                  choices=as.numeric(5), selected=5),
                      
                      selectInput(inputId = "Signi", label="Choose significance level", 
                                  choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05),
                       #radioButtons(inputId = "Signi", label = "Choose significance level", 
                       #             c("0.05" = 0.05, "0.01" = 0.01), selected = "0.05"),
                       
                       radioButtons(inputId = "TF", label = "Standardize", 
                                    c("True" = TRUE, "False" = FALSE), selected = FALSE),
                       actionButton("go", label = "Run")),
                      #h6(paste("*choices are", expression(T^0.25)))
                      #tags$h6(paste("*choices are T", HTML(paste(tags$sup("0.25"))), " T" ))
                      htmlOutput("Res"),
                      htmlOutput("Res1"),
                      ),
      
      conditionalPanel(condition = "input.tabselected==2",
                       fileInput("file_csv", "Choose CSV File",
                                 multiple = TRUE,
                                 accept = c("text/csv",
                                            "text/comma-separated-values,text/plain",
                                            ".csv")),
                       checkboxInput("header", "Header", TRUE),
                       numericInput(inputId = "Num2", label = "Choose number of observations per approximately stationary block*", 
                                   value = NULL, step = 1),
                       
                       numericInput(inputId = "Tapers2", label = "Choose number of tapers to use in multitaper spectral estimator**", 
                                   value = NULL, step = 1),
                       
                       selectInput(inputId = "Signi2", label="Choose significance level", 
                                   choices=as.numeric(seq(0.01,0.05, by=0.01)), selected = 0.05),
                       
                       radioButtons(inputId = "TF2", label = "Standardize", 
                                    c("True" = TRUE, "False" = FALSE), selected = FALSE),
                       actionButton("go2", label = "Run"), 
                       textOutput("res9"),
                       textOutput("res10"),
                       ), 
     
                      
    ),
    mainPanel(
      tabsetPanel(type = "tabs", id = 'tabselected', selected = 1,
                  tabPanel("Simulation Setting", value = 1),
                  tabPanel("File Upload", value = 2)),
      conditionalPanel(condition = "input.tabselected==1",
                       plotOutput("Image_Plot", height=1000),
                       radioButtons(inputId = "downloadType", label = "Select download type", choices = list("png","pdf")),
                       downloadButton('downloadData','Download the plot') , 
      ),
      conditionalPanel(condition = "input.tabselected==2",
                       plotOutput("Image_Plot2", height=1000)
      ), 
      )
  ))

server <- function(input,output, session) {
  observe({
    t2 <- as.numeric(input$Time)
    if(t2 == 500) {
      choices = 105
    } else if(t2 == 1000) {
      choices = c(31, 177)
    } else if(t2 == 5000) {
      choices = c(70, 594 )
    } else if(t2 == 10000) {
      choices = c(100, 1000)
    } else if(t2==20000) {
      choices = c(141, 1681)
    } else {
      choices=c(223, 3343)
    }
    updateSelectInput(session, "Num", choices = choices, selected = choices[1])
  })
  observe({
    t2 <- as.numeric(input$Time)
    sc_val <- (as.numeric(input$Num))
    if(t2 == 500) {
      tap=10
    } else if(t2 == 1000) {
      if(sc_val == 31){
        tap=5
      } else {
        tap=13
      }
    } else if(t2 == 5000) {
      if(sc_val == 70){
        tap=8
      } else {
        tap=c(24, 120)
      }
    } else if(t2 == 10000) {
      if(sc_val == 100) {
        tap=10
      } else {
        tap=c(31, 177)
      }
    } else if(t2 == 20000){
      if(sc_val == 141) {
        tap=11
      } else {
        tap=c(41, 262)
      }
    } else {
      if(sc_val == 223){
        tap=14
      } else {
        tap=c(57, 439)
      }
    }
    updateSelectInput(session, "Tapers",choices=tap, selected = tap[1])
  })
  output$Res <- renderText({
    paste(h6("*Choices are whichever of T", HTML(paste(tags$sup("0.25"))), ", T", HTML(paste(tags$sup("0.5"), 
                  ", T", HTML(paste(tags$sup("0.75"), "that are valid under the FRESH Procedure"))))))
    })
  output$Res1 <- renderText({
    paste(h6("**Choices are whichever of N", HTML(paste(tags$sup("0.5"))), ", N", HTML(paste(tags$sup("0.75"), 
                   ", N", HTML(paste(tags$sup("0.833"), "that are valid under the FRESH Procedure"))))))
  })
  output$res9 <- renderText({
    paste("*Please enter a dataframe")
  })
  observeEvent(input$file_csv, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    dataf <- read.csv(file$datapath, header = input$header)
    dataf <- dataf[[1]]
    updateNumericInput(session, "Num2", value=floor(sqrt(length(dataf))))
    updateNumericInput(session, "Tapers2", value=floor(0.15 * sqrt(length(dataf))))
    output$res9 <- renderText({
        paste("*Valid choices range from 30 to ", floor(length(dataf)/2))
    })
    output$res10 <- renderText({
          paste("**Valid choices range from 1 to ", floor(sqrt(length(dataf))*0.24))
    })
  })
  observeEvent(input$Num2, {
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
    bw <- floor((K+1)*(T_B/(T_B+1)))
    f_part <- c(1, floor(T_B/2 + 1))
    flo_mea <- floor(mean(f_part))
    max_tap <- ((flo_mea - 1) * ((T_B+1)/T_B)) - 1
    #if(max_val == bw){
    #  max_val <- max_val - 1
    #}
    output$res10 <- renderText({
      paste("**Valid choices range from 1 to ", floor(max_tap))
    })
    }
    
  })
  #observeEvent(inpit$Tapers2, {
  #  T_B <- input$Num2
  #  K <- input$Tapers2
  #  bw <- floor((K+1)*(T_B/(T_B+1)))
  #  max_val <- floor(T_B*0.24)
  #  if(max_val <= bw){
  #    max_val <- bw - 1
  #  }
  #  
  #})
  #observe({
  #  file <- input$file_csv
  #  ext <- tools::file_ext(file$datapath)
  #  le <- length(file[[1]])
  #  if(le == 0){
  #    output$res9 <- renderText({
  #      paste("*Please input a Dataframe")
  #    })
  #  } else {
  #  dataf <- read.csv(file$datapath, header = input$header)
  #  dataf <- dataf[[1]]
  #  #updateNumericInput(session, "Num2", value=floor(sqrt(length(dataf))))
  #  #updateNumericInput(session, "Tapers2", value=min(floor(2 * 0.15 * floor(sqrt(length(dataf))) - 1),
  #  #                                                 floor(length(dataf) ** 0.375)))
  # 
  #  output$res9 <- renderText({
  #    paste("*Valid choices range from 30 to ", length(dataf)/2)
  #  })
  #  gg <- input$Num2
  #  if(gg < 30) {
  #    updateNumericInput(session, "Num2", value=30)
  #  } 
  #  if (gg > (length(dataf)/2)) {
  #    updateNumericInput(session, "Num2", value=length(dataf)/2)
  #  }
  #  
  # }
  #})
  #observe({
  #  file <- input$file_csv
  #  ext <- tools::file_ext(file$datapath)
  #  le <- length(file[[1]])
  #  if(le == 0){
  #    output$res10 <- renderText({
  #      paste()
  #    })
  #  } else {
  #    dataf <- read.csv(file$datapath, header = input$header)
  #    dataf <- dataf[[1]]
  #    T_B <- input$Num2
  #    output$res10 <- renderText({
  #      paste("**Valid choices range from 1 to ", floor(T_B*0.24))
  #    })
  #    #(session, "Tapers2", value=floor(T_B ** (0.75)))
  #    #output$res10 <- renderText({
  #    #  paste("**Valid choices range from 1 to ", floor(2*T_B*.15-1))
  #    #})
  #    t2 <- input$Tapers2
  #    if(t2 < 1){
  #      updateNumericInput(session, "Tapers2", value = 1)
  #    }
  #    if(t2 > floor(T_B * 0.24)){
  #      updateNumericInput(session, "Tapers2", value=floor(T_B*0.24) - 1)
  #    }
  #     #if(t2 > floor(2*T_B*.15-1)) {
  #      #updateNumericInput(session, "Tapers2", value= (floor((2*T_B*.15-1) / 2)))
  #    #}
  #  }
  # 
  #})
  
  #observe({
  #  #res9
  #  if(exists("dataf")){
  #    updateTextInput(session, "res9", value=paste("Valid choices are any values between 30 and ", length(dataf)/2))
  #  } else {
  #    updateTextInput(session, "res9", value="Nooooo")
  #  }
  #})
  
  output$res10 <- renderText({
    
  })
  plot.list2 <- eventReactive(input$go2, ignoreNULL = FALSE, {
    file <- input$file_csv
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    dataf <- read.csv(file$datapath, header = input$header)
    dataf <- dataf[[1]]
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
               bigplot = c(.1, .55, .1, .5), smallplot = c(.6, .65, .1, .5));title(plot.list2()[[4]], line=0.75); 
    abline(h=plot.list2()[[5]], col = "skyblue", lwd=3); 
    
    vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                      just=c("left", "top"), y=0.5, x=0.6)
    len <- length(plot.list2()[[5]])
    vals <- plot.list2()[[5]]
    if(len == 0){
      str <- "(0, 0.5),"
    } else if (len == 1) {
      str <- paste("(0, ", round(vals, 3), "), (", round(vals, 3), ", 0.5),", sep="")
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
    pp <- data.frame("Predicted Partition Bands" = spp)
    colnames(pp) <- c("Predicted \n Partition Bands")
    grid.table(pp, vp=vp.br, rows=NULL)
    
    vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                     just=c("left", "top"), y=0.65, x=0.6)
    grid.polygon(x=c(0.31, 0.31,0.69, 0.69), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
    jj <- grid.legend(c("Predicted Partition"), gp=gpar(lty=1, col=c("skyblue")), vp=vp.r, 
                      draw=TRUE)
    
    print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.list2()[[6]])), y= plot.list2()[[6]])) + xlab("Time") +
             ylab("") + ggtitle("Inputted Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
             scale_x_continuous(limits=c(0,1), expand=c(0,0)), vp=vp.top)
    
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
    abline(h=plot.list()[[5]], col = "skyblue", lwd=3); abline(h=c(0.15, 0.35), col="red", lwd=3)
    
    vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                      just=c("left", "top"), y=0.5, x=0.6)
    act <- c("(0, 0.15)", "[0.15, 0.35)", "[0.35, 0.5)")
    len <- length(plot.list()[[5]])
    vals <- plot.list()[[5]]
    if(len == 0){
      str <- "(0, 0.5),"
    } else if (len == 1) {
      str <- paste("(0, ", round(vals, 3), "), (", round(vals, 3), ", 0.5),", sep="")
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
    pp <- data.frame("Actual Partition Bands" = act, "Predicted Partition Bands" = spp)
    colnames(pp) <- c("Actual \n Partition Bands", "Predicted \n Partition Bands")
    grid.table(pp, vp=vp.br, rows=NULL)
    
    vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.4, "npc"), 
                     just=c("left", "top"), y=0.65, x=0.6)
    grid.polygon(x=c(0.31, 0.31,0.69, 0.69), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
    jj <- grid.legend(c("Predicted Partition", "Actual Partition"), gp=gpar(lty=1, col=c("skyblue", "red")), vp=vp.r, 
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
      abline(h=plot.list()[[5]], col = "skyblue", lwd=3); abline(h=c(0.15, 0.35), col="red", lwd=3)
      
      vp.br <- viewport(height=unit(0.5, "npc"), width=unit(0.43, "npc"), 
                        just=c("left", "top"), y=0.5, x=0.65)
      act <- c("(0, 0.15)", "[0.15, 0.35)", "[0.35, 0.5)")
      len <- length(plot.list()[[5]])
      vals <- plot.list()[[5]]
      if(len == 0){
        str <- "(0, 0.5),"
      } else if (len == 1) {
        str <- paste("(0, ", round(vals, 3), "), (", round(vals, 3), ", 0.5),", sep="")
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
      pp <- data.frame("Actual Partition Bands" = act, "Predicted Partition Bands" = spp)
      colnames(pp) <- c("Actual \n Partition Bands", "Predicted \n Partition Bands")
      #colnames(pp) <- stringr::str_replace_all(colnames(df), "\\n", "<br>")
      grid.table(pp, vp=vp.br, rows=NULL)
      
      vp.r <- viewport(height=unit(0.5, "npc"), width=unit(0.43, "npc"), 
                       just=c("left", "top"), y=0.65, x=0.65)
      grid.polygon(x=c(0.31, 0.31,0.69, 0.69), y=c(0.6,0.4, 0.4,0.6 ), vp=vp.r)
      jj <- grid.legend(c("Predicted Partition", "Actual Partition"), gp=gpar(lty=1, col=c("skyblue", "red")), vp=vp.r, 
                        draw=TRUE)
      
      print( ggplot() + geom_line(aes(x=seq(0,1,length.out = length(plot.list()[[6]])), y= plot.list()[[6]])) + xlab("Time") +
               ylab("") + ggtitle("Simulated Time Series Data") + theme(plot.title = element_text(face="bold", hjust=0.5)) + 
               scale_x_continuous(limits=c(0,1), expand=c(0,0)), vp=vp.top)
      
      
      dev.off()
    }
  )
  
  
}

shinyApp(ui = ui, server = server)