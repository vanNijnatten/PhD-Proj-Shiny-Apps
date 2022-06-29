mainview <- shiny::tabPanel("Data",
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::h1("sidebar panel"),
      shiny::radioButtons("dataset", "Dataset",
        c("GLUCOLD", "")
      )
    ),
    shiny::mainPanel(
      h1("plot")
    )
  )
)

#shiny::fluidPage(
#  shiny::titlePanel("Main"),
#  shiny::sidebarLayout(
#    shiny::sidebarPanel(
#      shiny::h1("sidebar panel")
#    ),
#
#    shiny::mainPanel(
#      shiny::h1("main panel")
#      #shiny::tabsetPanel(
#      #  shiny::tabPanel(
#      #    "Volcano Plots",
#      #    shiny::plotOutput("volcano")
#      #  )
#      #)
#    )
#  )
#)
