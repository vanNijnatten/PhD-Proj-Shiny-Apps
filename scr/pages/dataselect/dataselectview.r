dataview <- shiny::tabPanel("Data",
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::radioButtons("dataset", "Dataset",
        c("GLUCOLD", "SHERLOCk", "NORM", "Stop Smoking")
      )
    ),
    shiny::mainPanel(
      h1("plot")
    )
  )
)
