library(shiny)

ui = fluidPage(
  titlePanel("Vancomycin TDM (Under Test)"),
  sidebarlayout(
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
      tags$hr()
      plotOutput("plot")
    )
  )
)

server = function(input, output, session){
  
  getDATAi = reactive({
    req(input$file1)
    DATAi = expandDATA(convDT(read.csv(input$file1$datapath, na.strings = c("",".","NA"), as.is =TRUE)))
    updateNumericInput(session, "TIME", value=ceiling(max(DATAi[,"TIME"])))
    DATAi[,c("ID","TIME","AMT","RATE","DV","MDV","SEX","AGE","BWT","SCR","CLCR")]
  })
  
  getEBE = reactive({
    DATAi = getDATAi()
    EBE(PredVanco, DATAi, TH, OM, SG)
  })
  
  getPI = reactive({
    DATAi = getDATAi()
    rEBE = getEBE()
  })
  
  
}