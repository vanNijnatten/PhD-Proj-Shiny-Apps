source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "dependencies.r"))
source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "helpers.r"))
source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "data.r"))

source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "pages", "main", "maincontroller.r"))
source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "pages", "main", "mainview.r"))

source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "pages", "dataselect", "dataselectcontroller.r"))
source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "pages", "dataselect", "dataselectview.r"))

source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "pages", "tableshow", "tableshowcontroller.r"))
source(verbose=TRUE, chdir=TRUE, file.path(here::here(), "scr", "pages", "tableshow", "tableshowview.r"))


shiny::shinyApp(
  navbarPage(">>",
    mainview#,
    #dataselectvview,
    #datashowview,
    #navbarMenu("More",
    #  tabPanel("Table",
    #    DT::dataTableOutput("table")
    #  ),
    #  tabPanel("About",
    #    fluidRow(
    #      includeMarkdown("about.md")
    #    )
    #  )
    #)
  ),

  function(input, output, session) {
    maincontroller(input, output, session)
    datacontroller(input, output, session)
    showdatacontroller(input, output, session)

    observeEvent(input$stop, {
      stopApp(cat("Stopped the app."))
    })
  }
)
