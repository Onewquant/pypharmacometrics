library(shiny)

ui = fluidPage(
  titlePanel("Vancomycin TDM (Under Test)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1","Choose CSV file",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
      br(),
      numericInput(inputId="TIME", label="Next Dosing Time", value=0, min=0),
      numericInput(inputId="AMT", label="Amount", value=1000, min=0),
      numericInput(inputId="RATE", label="Rate", value=2000, min=0),
      numericInput(inputId="II", label="Interdose Interval", value=12, min=0),
      numericInput(inputId="ADDL", label="Additional Counts", value=0, min=0),
      actionButton("goButton", "Go!"),
      br(),
      br(),
      uiOutput("xlim_slider")
      ),
    mainPanel(
      tableOutput("contents"),
      tags$hr(),
      plotOutput("plot")
    )
  )
)

server = function(input, output, session){
  
  # getDATAi = reactive({
  #   req(input$file1)
  #   DATAi = expandDATA(convDT(read.csv(input$file1$datapath, na.strings = c("",".","NA"), as.is =TRUE)))
  #   updateNumericInput(session, "TIME", value=ceiling(max(DATAi[,"TIME"])))
  #   DATAi[,c("ID","TIME","AMT","RATE","DV","MDV","SEX","AGE","BWT","SCR","CLCR")]
  # })

  getDATAi = reactive({
    req(input$file1)
    
    # raw_data = read.csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/sample_vanco_loading_maintenance.csv", na.strings = c("", ".", "NA"), as.is = TRUE)
    raw_data <- read.csv(input$file1$datapath, na.strings = c("", ".", "NA"), as.is = TRUE)
    data_prepped <- convDT(raw_data)           # TIME 컬럼 붙여줌 (numerical time으로)

    # 데이터에 ADDL 컬럼이 있고, 값이 실제로 > 0인 경우만 확장 적용
    if ("ADDL" %in% names(data_prepped) && any(!is.na(data_prepped$ADDL) & data_prepped$ADDL > 0)) {
      data_prepped <- expandDATA(data_prepped)
    }

    data_prepped <- data_prepped[order(data_prepped$ID, data_prepped$TIME), ]                          # ID, TIME 으로 정렬
    updateNumericInput(session, "TIME", value = ceiling(max(data_prepped$TIME)))                       # shiny session에서 INPUT 에 최대 시간보다 크도록 설정
    return(data_prepped[, c("ID","TIME","AMT","RATE","DV","MDV","SEX","AGE","BWT","SCR","CLCR")])      # prep data 반환
  })
  
  
  # DATAi = (data_prepped[, c("ID","TIME","AMT","RATE","DV","MDV","SEX","AGE","BWT","SCR","CLCR")])
  
  getEBE = reactive({
    DATAi = getDATAi()
    EBE(PredVanco, DATAi, TH, OM, SG)
  })
  
  # rEBE = EBE(PredVanco, DATAi, TH, OM, SG)
  
  getPI = reactive({
    DATAi = getDATAi()
    rEBE = getEBE()
    input$goButton
    isolate({
      calcTDM(PredVanco, DATAi, TH, SG, rEBE, input$TIME, input$AMT, input$RATE, input$II, input$ADDL) # AC: Additional Counts
    })
  })
  
  # PI = calcTDM(PredVanco, DATAi, TH, SG, rEBE, 50, 1000, 1000, 12, 10) # AC: Additional Counts
  
  output$contents = renderTable({
    DATAi = getDATAi()
    return(DATAi)
  })
  
  output$xlim_slider = renderUI({
    DATAi = getDATAi()
    minX = min(DATAi[,"TIME"])
    maxX = max(DATAi[,"TIME"]) + input$II*input$ADDL
    sliderInput("range","Time Range", min=minX, max=ceiling(maxX), value = c(minX,maxX))
  })
  
  # minX = min(DATAi[,"TIME"])
  # maxX = max(DATAi[,"TIME"]) + 12*10
  
  output$plot = renderPlot({
    PI = getPI()
    minX = min(PI$x)
    maxX = max(PI$x)
    if (input$range[1] != minX | input$range[2] != maxX) {
      xlm = input$range
    } else{
      xlm = c(minX, maxX)
    }
    
    plot(0, 0, type="n", xlab="Time", ylab="Concentration +/- 2SD", xlim=xlm, ylim=c(min(PI$ypiLL), max(PI$ypiUL)))
    points(PI$x[!is.na(PI$y)], PI$y[!is.na(PI$y)], pch=16)
    lines(PI$x, PI$y2, lty=1)
    lines(PI$x, PI$yciLL, lty=2, col="red")
    lines(PI$x, PI$yciUL, lty=2, col="red")
    lines(PI$x, PI$ypiLL, lty=3, col="blue")
    lines(PI$x, PI$ypiUL, lty=3, col="blue")
    abline(h=c(5,15,25,35), lty=2)
    })
}

source("C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/tdm_calculation_engine.R")
TH = c(3.8135955291021233, 39.889510090195238, 44.981835351176571, 2.0055189192561507)         # CL, V1, V2, Q?
OM = matrix(c( 0.10855133849022583,       -1.52445093837639736E-002,  -0.19698189309256298,       0.11914555131547180,
               -1.52445093837639736E-002, 3.01016276351715687E-003,   2.31285256338822770E-002,   -6.61597586201947800E-003,
               -0.19698189309256298,      2.31285256338822770E-002,   1.0420667710781930,         -0.19223488058085114,
               0.11914555131547180,       -6.61597586201947800E-003,  -0.19223488058085114,       0.416
               ), nrow = 4, byrow = TRUE)
SG = matrix(c(0.14019731106615912^2, 0.00000000,
              0.00000000,            1.8662549226759475^2), nrow=2)

shinyApp(ui, server)

