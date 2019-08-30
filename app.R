#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)
library(Rcpp)
sourceCpp("subsetsum.cpp")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Peptide Mass Checker"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         textInput("pept",
                     "Peptide String:"), 
         textInput("cap","Total Cap Mass: ",value=41.0265), 
         bsTooltip("cap", "Default assumes N-terminal acetyl group and C-terminal amide", placement = "top", trigger = "hover"),
         actionButton("trunc","List LCMS Output for Truncated Side Products"), 
         textInput("mass","Assumed Peptide Mass:"), 
         textInput("tol","Tolerance: ",value=1.5), 
         bsTooltip("tol","Search for side products with masses within tolerance from assumed mass"), 
         numericInput("ignore","# Initial Amino Acids Assumed Not to Fail to Couple:",value=1,min=1), 
         actionButton("search","Search For Mass in All Potential Side Products"), 
         bsTooltip("ignore", "Assumes synthesis starts from RHS of peptide", placement = "bottom", trigger = "hover"), 
         
         h6("Additional Options [Not Yet Functional]"), 
         checkboxInput("adj", "Protecting Group Mass Adjustments"),
         conditionalPanel(
           condition = "input.adj == true",
           textInput("pos", "Protecting Group Positions (Comma-Separated List)",placeholder = "e.g. 1,4,7"),
           textInput("madj","Protecting Group Masses (Comma-Separated List) ",placeholder= "e.g. 12.1,13.2,1.8"), 
           bsTooltip("pos", "Position 1 is first amino acid on LHS of peptide", placement = "bottom", trigger = "hover") 
         )
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         uiOutput("download_space"), 
         tableOutput("tableview")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  verifyComponents <- function(pept) { 
    verify = TRUE
    acids = c('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')
    pep.vec <- strsplit(pept,split="")[[1]]
    verify  <- sum(pep.vec %in% acids)== length(pep.vec) 
    return(verify) 
  }
   
  peptide <- reactive ({ 
      validate(
        need(input$pept != "", "Error: Please input a peptide")
      )
      validate(
        need(verifyComponents(input$pept), "Error: Peptide must contain only upper-case letters that each correspond to an amino acid")
      )
      input$pept
    })
  
  tolerance <- reactive({
    as.numeric(input$tol)
  })
  
  mass <- reactive({ 
    as.numeric(input$mass) 
  })
  
  ignore <- reactive({
    input$ignore 
  })
  
  cap <- reactive({
    as.numeric(input$cap) 
  })
  
  # adjustment_list <- reactive ( { 
  #   pept <- peptide()  
  #   adjustments <- rep(0,nchar(pept)) 
  #   adjustments 
  # })
    
   truncTable <- eventReactive(input$trunc, {
     pept <- peptide() 
     df <- listTruncated(pept)
     data.frame(df) 
    }
    )
  
   observeEvent( input$trunc, {
     output$tableview <- renderTable(truncTable(),
                          caption="Possible LCMS Output From Truncated Side-Products [(mass+k)/k]",
                          caption.placement = getOption("xtable.caption.placement", "top"))
     output$download_space <- renderUI({
       downloadButton("downloadtrunc","Download as CSV")
     })
     
   },
    once=TRUE)
   
   searchTable <- eventReactive(input$search,{
     pept <- peptide() 
     mass <- mass() 
     tol <- tolerance() 
     ignore <- ignore() 
     cap <- cap() 
     df <- searchPeptides(pept,mass,tol,ignore,cap)
     data.frame(df) 
    })
   
   observeEvent(input$search, { 
     output$tableview <- renderTable(searchTable(),
                                     caption="Synthesis Side Products with Matching Masses",
                                     caption.placement = getOption("xtable.caption.placement", "top"))
      output$download_space <- renderUI({
       downloadButton("downloadsearch","Download as CSV") 
      }) 
   })
   
   # output$downloadsearch <-downloadHandler(
   #   filename=function() { 
   #     
   #  }
   # )
   output$downloadtrunc <- downloadHandler(  
      filename = function() { 
        paste(peptide(), "_truncated", ".csv", sep="")
      }, 
      content = function(file) { 
        write.csv(truncTable(),file,row.names=FALSE)
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

