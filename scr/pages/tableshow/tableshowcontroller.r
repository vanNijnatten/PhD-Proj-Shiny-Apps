showdatacontroller <- function(input, output, session) {
  #widget <- DT::datatable(
  #  dat = data.frame(x=1:5, y=6:10),
  #  selection = "none",
  #  rownames = FALSE,
  #  filter = list(position = "top"),
  #  escape = FALSE,
  #  extensions = "KeyTable",
  #  options = list(
  #    keys = TRUE,
  #    search = list(regex = TRUE),
  #    columnDefs = list(
  #      list(
  #        orderSequence = c("desc", "asc"),
  #        targets = "_all"
  #      ),
  #      list(
  #        className = "dt-center",
  #        targets = "_all"
  #      )
  #    ),
  #    autoWidth = TRUE,
  #    processing = FALSE,
  #    pageLength = 10,
  #    lengthMenu = list(c(5, 10, 25, 50, -1), c("5", "10", "25", "50", "All"))
  #  )
  #)
  #
  #output$tbl <- DT::renderDataTable(widget)
}
