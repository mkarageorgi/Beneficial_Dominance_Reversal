StatTable <- function(model, directory, display_title, save_title, pred_labels, dv_labels) {
  
  library('webshot')
  library('sjPlot')
  
  # Set default value for pred_labels if none provided
  if (missing(pred_labels) || length(pred_labels) == 0) {
    pred_labels <- NULL
  }
  
  # Set default value for dv_labels if none provided
  if (missing(dv_labels) || length(dv_labels) == 0) {
    dv_labels <- NULL 
  }
  
  save_dir <- file.path(directory)
  
  # Format the save title to lowercase and replace spaces with underscores
  formatted_save_title <- gsub(" ", "_", tolower(save_title))
  file_name <- file.path(save_dir, sprintf("%s", formatted_save_title))
  
  # Generate table
  table_result <- sjPlot::tab_model(
    model,
    pred.labels = pred_labels,
    title = display_title,  # Use the display title
    show.aic = TRUE,
    dv.labels = dv_labels,
    p.style = "numeric_stars",
    show.ci = FALSE,
    show.se = TRUE,
    show.stat = TRUE,
    use.viewer = FALSE,  # Disable viewer for better compatibility with scrolling
    CSS = list(
      css.depvarhead = 'font-weight: bold;',
      css.fontsize = 'font-size: 690px;',
      css.table = 'border-collapse: collapse; width: 100%; font-family: Arial;',
      css.tr = 'margin-bottom: 5px;',
      css.td = 'text-align: center;'  # Center-align the text within table cells
    )
  )
  
  # If the table is long, make it scrollable
  if (length(pred_labels) > 10) {
    scrollable_table <- paste0(
      '<div style="overflow-x:auto; max-height:400px; overflow-y:scroll;">',
      table_result$knitr,
      '</div>'
    )
    cat(scrollable_table, file = paste0(file_name, '_scroll.html'))  # Save with .html extension
  }
  
  # Save the main table with .html extension
  cat(table_result$knitr, file = paste0(file_name, '.html'))
  
  # Capture a screenshot of the HTML file
  webshot::webshot(
    paste0(file_name, '.html'),  # Ensure .html extension is included
    file.path(save_dir, sprintf("%s.pdf", formatted_save_title))
  )
}
