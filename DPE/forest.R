
#' Forest plot
#'
#' A data frame will be used for the basic layout of the forest plot.
#' Graphical parameters can be set using the \code{\link{forest_theme}} function.
#'
#' @param data Data to be displayed in the forest plot
#' @param est Point estimation. Can be a list for multiple columns
#' and/or multiple groups. If the length of the list is larger than
#' then length of \code{ci_column}, then the values reused for each column
#' and considered as different groups.
#' @param lower Lower bound of the confidence interval, same as \code{est}.
#' @param upper Upper bound of the confidence interval, same as \code{est}.
#' @param sizes Size of the point estimation box, can be a unit, vector or a list.
#' @param ref_line X-axis coordinates of zero line, default is 1. Provide an atomic
#'  vector if different reference line for each \code{ci_column} is desired.
#' @param vert_line Numerical vector, add additional vertical line at given value.
#' Provide a list of numerical vector element if different vertical line for each
#'  \code{ci_column} is desired.
#' @param ci_cols Color of the point estimation box and confidence intervals, can be a unit or vector .
#' @param ci_column Column number of the data the CI will be displayed.
#' @param xlog If TRUE, x-axis tick marks assume values are exponential, e.g.
#' for logistic regression (OR), survival estimates (HR), Poisson regression etc.
#' Provide a logical vector if different conversion for each \code{ci_column} is
#' desired.
#' @param is_summary A logical vector indicating if the value is a summary value,
#' which will have a diamond shape for the estimate. Can not be used with multiple
#' group `forestplot`.
#' @param xlim Limits for the x axis as a vector of length 2, i.e. c(low, high). It
#' will take the minimum and maximum of the lower and upper value if not provided.
#' This will apply to all CI columns if provided, and will be calculated automatically
#' for each column if not provided. This should be a list with the same length of
#' \code{ci_column} if different \code{xlim} for different column is desired.
#' @param ticks_at Set X-axis tick-marks point. This will apply to all CI columns if
#' provided, and will be calculated automatically for each column if not provided.
#' This should be a list if different \code{ticks_at} for different column is desired.
#' @param ticks_digits Number of digits for the x-axis, default is 1. This should be
#' a numerical vector if different rounding will be applied to different column.
#' @param arrow_lab Labels for the arrows, string vector of length two (left and
#' right). The theme of arrow will inherit from the x-axis. This should be a list
#' if different arrow labels for each column is desired.
#' @param xlab X-axis labels, it will be put under the x-axis. An atomic vector should 
#' be provided if different \code{xlab} for different column is desired.
#' @param footnote Footnote for the forest plot, will be aligned at left bottom
#' of the plot. Please adjust the line length with line break to avoid the overlap
#' with the arrow and/or x-axis.
#' @param title The text for the title.
#' @param nudge_y Horizontal adjustment to nudge groups by, must be within 0 to 1.
#' @param theme Theme of the forest plot, see \code{\link{forest_theme}} for
#' details.
#'
#' @return A \code{\link[gtable]{gtable}} object.
#' @example inst/examples/forestplot-example.R
#' @export
#'
#'

forest <- function(data,
                   est,
                   lower,
                   upper,
                   sizes = 0.4,
                   ref_line = ifelse(xlog, 1, 0),
                   vert_line = NULL,
                   ci_cols = "black",
                   ci_column,
                   xlog = FALSE,
                   is_summary = NULL,
                   xlim = NULL,
                   ticks_at = NULL,
                   ticks_digits = 0,
                   arrow_lab = NULL,
                   xlab = NULL,
                   footnote = NULL,
                   title = NULL,
                   nudge_y = 0,
                   theme = NULL){

  check_errors(data = data, est = est, lower = lower, upper = upper, sizes = sizes, 
               ref_line = ref_line, vert_line = vert_line, ci_cols = ci_cols, ci_column = ci_column,
               xlog = xlog, is_summary = is_summary, xlim = xlim, ticks_at = ticks_at,
               ticks_digits = ticks_digits, arrow_lab = arrow_lab, xlab = xlab,
               title = title)

  # Point sizes
  if(length(sizes) == 1 & !inherits(sizes, "list"))
    sizes <- rep(sizes, nrow(data))
  #ci line colors
  
  if(length(ci_cols) == 1 & !inherits(ci_cols, "list"))
    ci_cols <- rep(ci_cols, nrow(data))
  
  # Set theme
  if(is.null(theme)){
    theme <- forest_theme()
  }

  # For multiple ci_column
  if(length(ref_line) == 1)
    ref_line <- rep(ref_line, length(ci_column))

  if(!is.null(vert_line) && !inherits(vert_line, "list"))
    vert_line <- rep(list(vert_line), length(ci_column))

  if(length(xlog) == 1)
    xlog <- rep(xlog, length(ci_column))

  if(!is.null(xlim) && !inherits(xlim, "list"))
    xlim <- rep(list(xlim), length(ci_column))

  if(!is.null(ticks_at) && !inherits(ticks_at, "list"))
    ticks_at <- rep(list(ticks_at), length(ci_column))

  if(length(ci_column) != length(ticks_digits))
    ticks_digits <- rep(ticks_digits, length(ci_column))

  if(!is.null(arrow_lab) && !inherits(arrow_lab, "list"))
    arrow_lab <- rep(list(arrow_lab), length(ci_column))

  if(length(xlab) == 1)
    xlab <- rep(xlab, length(ci_column))

  if(inherits(est, "list") | inherits(lower, "list") | inherits(upper, "list")){

    if(!inherits(est, "list") | !inherits(lower, "list") | !inherits(upper, "list"))
      stop("Estimate, lower and upper must be a list if plan to plot CI for multiple columns.")

    # Calculate group number
    group_num <- length(est)/length(ci_column)
    ci_col_list <- rep(ci_column, group_num)

    # Replicate sizes to align with est and CIs
    if(!inherits(sizes, "list"))
      sizes <- rep(list(sizes), length(ci_col_list))

    theme <- make_group_theme(theme = theme, group_num = group_num)

    # Get color and pch
    color_list <- rep(theme$ci$col, each = length(ci_column))
    pch_list <- rep(theme$ci$pch, each = length(ci_column))
    lty_list <- rep(theme$ci$lty, each = length(ci_column))
    lwd_list <- rep(theme$ci$lwd, each = length(ci_column))

    # Check nudge_y
    if(nudge_y >= 1 || nudge_y < 0)
      stop("`nudge_y` must be within 0 to 1.")

    # Check nudge_y
    if(group_num > 1 & nudge_y == 0)
      nudge_y <- 0.1

    # Create nudge_y vector
    if(group_num > 1){
      if((group_num %% 2) == 0){
        rep_tm <- cumsum(c(nudge_y/2, rep(nudge_y, group_num)))
        nudge_y <- c(rep_tm[1:(group_num/2)], -rep_tm[1:(group_num/2)])
      }else{
        rep_tm <- cumsum(c(0, rep(nudge_y, group_num)))
        nudge_y <- c(0, rep_tm[2:(group_num/2)], -rep_tm[2:(group_num/2)])
      }
    }

    nudge_y <- rep(nudge_y, each = length(ci_column))

  }else{
    if(length(unique(c(length(est), length(lower), length(upper), nrow(data)))) != 1)
      stop("Estimate, lower and upper should be same length as data row number.")

    #if(!is.null(sizes) & length(sizes) != 1){
    #  sizes <- sqrt(1/sizes)
    #  sizes <- sizes/max(sizes, na.rm = TRUE)
    #}

    ci_col_list <- ci_column
    color_list <- theme$ci$col
    pch_list <- theme$ci$pch
    lty_list <- theme$ci$lty
    lwd_list <- theme$ci$lwd

    group_num <- 1

    est <- list(est)
    lower <- list(lower)
    upper <- list(upper)
    sizes <- list(sizes)

  }

  if(group_num > 1 || is.null(is_summary))
    is_summary <- rep(FALSE, nrow(data))

  # Positions of values in ci_column
  gp_list <- rep_len(1:(length(lower)/group_num), length(lower))

  # Check exponential
  for(i in seq_along(xlog)){
    if(xlog[i]){
      sel_num <- gp_list == i
      if (any(unlist(est[sel_num]) <= 0, na.rm = TRUE) ||
            any(unlist(lower[sel_num]) <= 0, na.rm = TRUE) ||
            any(unlist(upper[sel_num]) <= 0, na.rm = TRUE) ||
            (any(ref_line[i] <= 0)) ||
            (!is.null(vert_line) && any(unlist(vert_line[[i]]) <= 0, na.rm = TRUE)) ||
            (!is.null(xlim) && any(unlist(xlim[[i]]) < 0))) {
        stop("est, lower, upper, ref_line, vert_line and xlim should be larger than 0, if `xlog=TRUE`.")
      }
    }
  }

  # Set xlim to minimum and maximum value of the CI
  xlim <- lapply(seq_along(ci_column), function(i){
    sel_num <- gp_list == i
    make_xlim(xlim = xlim[[i]],
              lower = lower[sel_num],
              upper = upper[sel_num],
              ref_line = ref_line[i],
              ticks_at = ticks_at[[i]],
              is_exp = xlog[i])
  })

  # Set X-axis breaks if missing
  ticks_at <- lapply(seq_along(xlim), function(i){
    make_ticks(at = ticks_at[[i]],
               xlim = xlim[[i]],
               refline = ref_line[i],
               is_exp = xlog[i])
  })

  # Calculate heights
  col_height <- apply(data, 1, function(x){
                        max(convertHeight(stringHeight(x), "mm", valueOnly = TRUE))
                      })
  col_height <- unit(col_height, "mm")

  # Add increase heights for multiple groups
  if(group_num > 1){
    heights <- group_num*0.7*col_height + theme$tab_theme$core$padding[2]
    # Convert data to plot
    gt <- tableGrob(data,
                    theme = theme$tab_theme,
                    heights = heights,
                    rows = NULL)
  }else{
    gt <- tableGrob(data,
                    theme = theme$tab_theme,
                    rows = NULL)
  }

  # Column index
  col_indx <- rep_len(1:length(ci_column), length(ci_col_list))
  #browser()
  # Draw CI
  for(col_num in seq_along(ci_col_list)){

    if(xlog[col_indx[col_num]]){
      est[[col_num]] <- log(est[[col_num]])
      lower[[col_num]] <- log(lower[[col_num]])
      upper[[col_num]] <- log(upper[[col_num]])
    }

    for(i in 1:nrow(data)){
      if(is.na(est[[col_num]][i]))
        next

      if(is_summary[i]){
        draw_ci <- make_summary(est = est[[col_num]][i],
                                lower = lower[[col_num]][i],
                                upper = upper[[col_num]][i],
                                size = sizes[[col_num]][i],
                                xlim = xlim[[col_indx[col_num]]],
                                gp = theme$summary)
      }else {
        draw_ci <- makeci(est = est[[col_num]][i],
                          lower = lower[[col_num]][i],
                          upper = upper[[col_num]][i],
                          size = sizes[[col_num]][i],
                          xlim = xlim[[col_indx[col_num]]],
                          pch = pch_list[col_num],
                          gp = gpar(lty = lty_list[col_num],
                                    lwd = lwd_list[col_num],
                                    col = ci_cols[i]),
                          t_height = theme$ci$t_height,
                          nudge_y = nudge_y[col_num])
      }

      gt <- gtable_add_grob(gt, draw_ci,
                            t = i + 1,
                            l = ci_col_list[col_num],
                            b = i + 1,
                            r = ci_col_list[col_num],
                            clip = "off",
                            name = paste0("ci-", i, "-", col_num))
    }
  }

  tot_row <- nrow(gt)

  # Prepare X axis
  x_axis <- lapply(seq_along(xlim), function(i){
    make_xaxis(at = ticks_at[[i]],
               gp = theme$xaxis,
               ticks_digits = ticks_digits[i],
               x0 = ref_line[i],
               xlim = xlim[[i]],
               xlab = xlab[i],
               is_exp = xlog[i])
  })

  x_axht <- sapply(x_axis, function(x){
    convertHeight(sum(grobHeight(x$children)) + unit(.5, "lines"),
                  unitTo = "mm",
                  valueOnly = TRUE)
  })

  gt <- gtable_add_rows(gt, heights = unit(max(x_axht), "mm"))

  # Add footnote
  if(!is.null(footnote)){
    footnote_grob <- textGrob(label = footnote,
                              gp = theme$footnote,
                              x = 0,
                              just = c("left", "top"),
                              check.overlap = TRUE,
                              name = "footnote")

    gt <- gtable_add_grob(gt,
                          footnote_grob,
                          t = tot_row + 1,
                          l = 1,
                          b = tot_row + 1, r = min(ci_column),
                          clip = "off",
                          name = "footnote")
  }

  # Prepare arrow object and row to put it
  if(!is.null(arrow_lab)){
    arrow_grob <- lapply(seq_along(xlim), function(i){
      make_arrow(x0 = ref_line[i],
                 arrow_lab = arrow_lab[[i]],
                 gp = theme$xaxis,
                 is_exp = xlog[i],
                 xlim = xlim[[i]])
    })

    lb_ht <- sapply(arrow_grob, function(x){
      convertHeight(max(grobHeight(x$children)) + unit(.5, "lines"),
                    unitTo = "mm",
                    valueOnly = TRUE)
    })

    gt <- gtable_add_rows(gt, heights = unit(max(lb_ht), "mm"))

  }

  for(j in ci_column){
    idx <- which(ci_column == j)
    # Add reference line
    gt <- gtable_add_grob(gt,
                          vert_line(x = ref_line[idx],
                                    gp = theme$refline,
                                    xlim = xlim[[idx]],
                                    is_exp = xlog[idx]),
                          t = 2,
                          l = j,
                          b = tot_row, r = j,
                          clip = "off",
                          name = paste0("reference.line-", j))

    # Add the X-axis
    gt <- gtable_add_grob(gt, x_axis[[idx]],
                          t = tot_row + 1,
                          l = j,
                          b = tot_row + 1, r = j,
                          clip = "off",
                          name = paste0("xaxis-", j))

    # Add vertical line
    if(!is.null(vert_line))
      gt <- gtable_add_grob(gt,
                            vert_line(x = vert_line[[idx]],
                                      gp = theme$vertline,
                                      xlim = xlim[[idx]],
                                      is_exp = xlog[idx]),
                            t = 2,
                            l = j,
                            b = tot_row, r = j,
                            clip = "off",
                            name = paste0("vert.line-", j))

    # Add arrow
    if(!is.null(arrow_lab))
      gt <- gtable_add_grob(gt, arrow_grob[[idx]],
                            t = nrow(gt), l = j,
                            b = nrow(gt), r = j,
                            clip = "off",
                            name = paste0("arrow-", j))

  }

  # Add legend
  if(group_num > 1){

    by_row <- if(!theme$legend$position %in% c("top", "bottom") || is.null(theme$legend$position)) TRUE else FALSE

    legend <- theme$legend
    legend$pch <- theme$ci$pch
    legend$color <- theme$ci$col
    legend$lty <- theme$ci$lty


    leg_grob <- do.call(legend_grob, legend)

    if(by_row){
      gt <- gtable_add_cols(gt, widths = max(grobWidth(leg_grob$children)) + unit(.5, "lines"))
      gt <- gtable_add_grob(gt, leg_grob,
                            t = 2, l = ncol(gt),
                            b = nrow(gt)-1, r = ncol(gt),
                            clip = "off",
                            name = "legend")
    }else{
      add_pos <- ifelse(legend$position == "top", 0, -1)
      gt <- gtable_add_rows(gt, heights = max(grobHeight(leg_grob$children)) + unit(.5, "lines"), pos = add_pos)
      gt <- gtable_add_grob(gt, leg_grob,
                            t = if(add_pos == 0) 1 else nrow(gt), l = 1,
                            b = if(add_pos == 0) 1 else nrow(gt), r = ncol(gt),
                            clip = "off",
                            name = "legend")
    }
  }

  if(!is.null(title)){
    max_height <- max(convertHeight(stringHeight(title), "mm", valueOnly = TRUE))
    gt <- gtable_add_rows(gt, unit(max_height, "mm") + unit(2, "mm"), pos = 0)
    title_x <- switch(theme$title$just,
                      right = unit(1, "npc"),
                      left  = unit(0, "npc"),
                      center = unit(.5, "npc"))
    title_gb <- textGrob(label = title,
                         gp = theme$title$gp,
                         x = title_x,
                         just = theme$title$just,
                         check.overlap = TRUE,
                         name = "plot.title")

    gt <- gtable_add_grob(gt, title_gb,
                          t = 1,
                          b = 1,
                          l = 1,
                          r = ncol(gt),
                          clip = "off",
                          name = "plot.title")
  }

  # Add padding
  gt <- gtable_add_padding(gt, unit(0, "mm"))

  class(gt) <- union("forestplot", class(gt))

  return(gt)

}


#' Draw plot
#'
#' Print or draw forestplot.
#'
#' @param x forestplot to display
#' @param autofit If true, the plot will be autofit.
#' @param ... other arguments not used by this method
#' @return Invisibly returns the original forestplot.
#' @rdname print.forestplot
#' @method print forestplot
#' @export
print.forestplot <- function(x, autofit = FALSE, ...){

  if(autofit){
    # Auto fit the page
    x$widths <- unit(rep(1/ncol(x), ncol(x)), "npc")
    x$heights <- unit(rep(1/nrow(x), nrow(x)), "npc")
  }

  grid.newpage()
  grid.draw(x)

  invisible(x)
}

#' @method plot forestplot
#' @rdname print.forestplot
#' @export
plot.forestplot <- print.forestplot




#' Add text to forest plot
#'
#' This function can be used to add text to forest plot. The text can be span to
#' multiple rows and columns. The height of the row will be changed accordingly
#' if the text is added to only one row. The width of the text may exceeds the
#' columns provided if the text is too long.
#'
#' @param plot A forest plot object.
#' @param text A character or expression vector, see \code{\link[grid]{textGrob}}.
#' @param row Row to add the text, this will be ignored if the \code{part} is
#' "header".
#' @param col A numeric value or vector indicating the columns the text will be
#' added. The text will span over the column if a vector is given.
#' @param part Part to add text, body (default) or header.
#' @param just The justification of the text, \code{"center"} (default),
#' \code{"left"} or \code{"right"}.
#' @param gp An object of class \code{"gpar"}, this is the graphical parameter
#'  settings of the text. See \code{\link[grid]{gpar}}.
#' @param padding Padding of the text, default is \code{unit(1, "mm")}
#'
#' @return A \code{\link[gtable]{gtable}} object.
#' 
#' @export
#'
add_text <- function(plot,
                     text,
                     row = NULL,
                     col = NULL,
                     part = c("body", "header"),
                     just = c("center", "left", "right"),
                     gp = gpar(),
                     padding = unit(3,"mm")){
  
  if(!inherits(plot, "forestplot"))
    stop("plot must be a forestplot object.")
  
  if(!inherits(gp, "gpar"))
    stop("gp must be a gpar object.")
  
  if(!is.unit(padding))
    padding <- unit(padding, "mm")
  
  part <- match.arg(part)
  just <- match.arg(just)
  
  # Row must be provided for the body
  if(part == "body" & is.null(row))
    stop("Row must be defined if the text is interting to body.")
  
  # Align text
  tx_x <- switch(just,
                 right = unit(1, "npc") - padding,
                 left  = padding,
                 center = unit(.5, "npc"))
  
  l <- plot$layout
  
  # Header
  if(part == "header"){
    if(!is.null(row))
      row <- row + min(l$b[which(l$name == "colhead-fg")]) - 1
    else
      row <- max(l$b[which(l$name == "colhead-fg")])
  }else{
    row <- max(l$b[which(l$name == "colhead-fg")]) + row
  }
  
  # Span to whole plot if col is missing
  if(is.null(col))
    col <- 2:max(l$r)
  else
    col <- 1 + col # Add 1 to account for padding of the plot
  
  txt_grob <- textGrob(label = text,
                       gp = gp,
                       x = tx_x,
                       just = just,
                       check.overlap = TRUE,
                       name = "custom-text.add")
  
  # Change height of the text is added to one row only,
  # and the height of the text is larger than then row
  if(length(row) == 1 ){
    txt_height <- convertHeight(grobHeight(txt_grob), "mm", valueOnly = TRUE)
    row_height <- convertHeight(plot$heights[row], "mm", valueOnly = TRUE)
    
    if(txt_height > row_height)
      plot$heights[row] <- grobHeight(txt_grob)
    
  }
  
  plot <- gtable_add_grob(plot, txt_grob,
                          t = min(row),
                          b = max(row),
                          l = min(col),
                          r = max(col),
                          clip = "off",
                          name = "text.add")
  
  return(plot)
  
}




#' Add underline to cells
#'
#' This function can be used to add underline to cells.
#'
#' @param plot A forest plot object.
#' @param row A numeric value or vector indicating row number to add underline.
#' This is corresponding to the data row number. Remember to account for any
#' text inserted. This will be ignored if the \code{part} is  "header" and the
#' underline will be drawn under the header column.
#' @param col A numeric value or vector indicating the columns to add underline.
#' @param part The underline will be added to \code{"body"} (default) or
#' \code{"header"}.
#' @param gp An object of class \code{"gpar"}, graphical parameter to be passed
#' to \code{\link[grid]{segmentsGrob}}.
#'
#' @return A \code{\link[gtable]{gtable}} object.
#' 
#' @export
#'
add_underline <- function(plot,
                          row = NULL,
                          col = NULL,
                          part = c("body", "header"),
                          gp = gpar(lwd = 2.0)){
  
  if(!inherits(plot, "forestplot"))
    stop("plot must be a forestplot object.")
  
  if(!inherits(gp, "gpar"))
    stop("gp must be a gpar object.")
  
  part <- match.arg(part)
  
  if(part == "body" & is.null(row))
    stop("row must be provided for the body.")
  
  l <- plot$layout
  
  # Header
  if(part == "header"){
    if(is.null(row))
      row <- max(l$b[which(l$name == "colhead-fg")])
    else
      row <- row + min(l$b[which(l$name == "colhead-fg")]) - 1
  }else{
    row <- max(l$b[which(l$name == "colhead-fg")]) + row
  }
  
  # Span to whole plot if col is missing
  if(is.null(col))
    col <- 2:max(l$r)
  else
    col <- 1 + col # Add 1 to account for padding of the plot
  
  for(i in seq_along(row)){
    
    seg_gb <- segmentsGrob(x0 = unit(0,"npc"),
                           y0 = unit(0,"npc"),
                           x1 = unit(1,"npc"),
                           y1 = unit(0,"npc"),
                           gp = gp,
                           name = paste("underline", row[i], sep = "-"))
    
    plot <- gtable_add_grob(plot, seg_gb,
                            t = row[i],
                            b = row[i],
                            l = min(col),
                            r = max(col),
                            name = "underline")
  }
  
  return(plot)
  
  
}




#' Checking error for forest plot
#'
#' @inheritParams forest
#'
#' @keywords internal
#'
check_errors <- function(data,
                         est,
                         lower,
                         upper,
                         sizes,
                         ref_line,
                         vert_line,
                         ci_cols,
                         ci_column,
                         xlog ,
                         is_summary,
                         xlim,
                         ticks_at,
                         ticks_digits,
                         title,
                         arrow_lab,
                         xlab){
  
  if(!is.numeric(ci_column))
    stop("ci_column must be numeric atomic vector.")
  
  if(!is.null(title) && length(title) != 1)
    stop("title must be of length 1.")
  
  # Check length
  if(length(unique(c(length(est), length(lower), length(upper)))) != 1)
    stop("Estimate, lower and upper should have the same length.")
  
  # Check length for the summary
  if(!is.null(is_summary) && length(is_summary) != nrow(data))
    stop("is_summary should have same legnth as data rownumber.")
  
  # Check ref_line
  if(!is.numeric(ref_line) || !length(ref_line) %in% c(1, length(ci_column)))
    stop("ref_line should be of length 1 or the same length as ci_column.")
  
  # Check the xlog
  if(!is.logical(xlog) || !length(xlog) %in% c(1, length(ci_column)))
    stop("xlog must be logical and of length 1 or the same length as ci_column.")
  
  # Check the xlab
  if(!is.null(xlab) && !length(xlab) %in% c(1, length(ci_column)))
    stop("xlab must be of length 1 or the same length as ci_column.")
  
  # Check tick_digits
  if(!is.numeric(ticks_digits) || !length(ticks_digits) %in% c(1, length(ci_column)))
    stop("ticks_digits must be numeric of length 1 or same length as ci_column.")
  
  # If only one CI column
  if(length(ci_column) == 1){
    
    # Check vertical line
    if(!is.null(vert_line) && !is.numeric(vert_line))
      stop("vert_line must be a numeric vector.")
    
    # Check arrow
    if(!is.null(arrow_lab) & length(arrow_lab) != 2)
      stop("Arrow label must of length 2.")
    
    # Check xlim
    if(!is.null(xlim) && (!is.numeric(xlim) || length(xlim) != 2 || xlim[1] >= xlim[2]))
      stop("xlim must be numeric and of length 2, with first element less than the second.")
    
    # Check the break
    if(!is.null(ticks_at) && !is.numeric(ticks_at))
      stop("ticks_at must be numeric.")
    
    if(!is.null(ticks_at) && !is.null(xlim)){
      if(max(ticks_at) > max(xlim) || min(ticks_at) < min(xlim))
        warning("ticks_at is outside the xlim.")
    }
    
  }else{
    
    # Check vertical line
    if(!is.null(vert_line)){
      if(inherits(vert_line, "list")){
        if(length(vert_line) != length(ci_column))
          stop("vert_line must have the same length as ci_column.")
        cl <- sapply(vert_line, is.numeric)
        if(any(!cl))
          stop("vert_line must be all numeric.")
      }else {
        if(!is.numeric(vert_line))
          stop("vert_line must be a numeric vector.")
      }
    }
    
    # Check arrow
    if(!is.null(arrow_lab)){
      if(inherits(arrow_lab, "list")){
        if(length(arrow_lab) != length(ci_column))
          stop("arrow_lab must have the same length as ci_column.")
        cl <- sapply(arrow_lab, length) == 2
        if(any(!cl))
          stop("Elements in the arrow_lab must of length 2.")
      }else {
        if(!is.null(arrow_lab) & length(arrow_lab) != 2)
          stop("Arrow label must of length 2.")
      }
    }
    
    # Check xlim
    if(!is.null(xlim)){
      if(inherits(xlim, "list")){
        if(length(xlim) != length(ci_column))
          stop("xlim must have the same length as ci_column.")
        tst <- sapply(xlim, function(x){
          !is.numeric(x) || length(x) != 2 || x[1] >= x[2]
        })
        if(any(tst))
          stop("Elements in the xlim must be numeric and of length 2, with first element less than the second.")
        
      }else {
        if(!is.numeric(xlim) || length(xlim) != 2 || xlim[1] >= xlim[2])
          stop("xlim must be numeric and of length 2, with first element less than the second.")
      }
    }
    
    # Check the break
    if(!is.null(ticks_at)){
      if(inherits(ticks_at, "list")){
        if(length(ticks_at) != length(ci_column))
          stop("ticks_at must have the same length as ci_column.")
        
        cl <- sapply(ticks_at, is.numeric)
        if(any(!cl))
          stop("Elements in the ticks_at must be numeric.")
        
      }else {
        if(!is.numeric(ticks_at))
          stop("ticks_at must be numeric.")
      }
    }
    
  }
  
}


#' Edit forest plot
#'
#' This function is used to edit the graphical parameter of text and background
#' of the forest plot.
#'
#' @param plot A forest plot object.
#' @param row A numeric value or vector indicating row number to edit in the
#' dataset. Will edit the whole row if left blank for the body. This will be
#' ignored if the \code{part} is "header".
#' @param col A numeric value or vector indicating column to edit in the dataset.
#'  Will edit the whole column if left blank.
#' @param part Part to edit, body (default) or header.
#' @param which Which element to edit, text or background of the cell.
#' @param gp Pass \code{gpar} parameters, see \code{\link[grid]{gpar}}. It should
#' be passed as \code{gpar(col = "red")}.
#' 
#' @return A \code{\link[gtable]{gtable}} object.
#'
#' @export
#'
edit_plot <- function(plot,
                      row = NULL,
                      col = NULL,
                      part = c("body", "header"),
                      which = c("text", "background"),
                      gp){
  
  if(!inherits(plot, "forestplot"))
    stop("plot must be a forestplot object.")
  
  if(!inherits(gp, "gpar"))
    stop("gp must be a gpar object.")
  
  part <- match.arg(part)
  which <- match.arg(which)
  
  part <- switch(part,
                 body = "core",
                 header = "colhead")
  which <- switch(which,
                  text = "fg",
                  background = "bg")
  
  name_to_edit <- paste(part, which, sep = "-")
  
  l <- plot$layout
  
  # If body
  if(part == "core"){
    # Add number of header to the row
    if(!is.null(row))
      row <- row + max(l$b[which(l$name == "colhead-fg")])
    # Apply to whole body if missing row
    else
      row <- unique(l$b[which(l$name == "core-fg")])
  }else {
    # If header, add header part
    if(!is.null(row))
      row <- row + min(l$b[which(l$name == "colhead-fg")]) - 1
    else
      row <- min(l$b[which(l$name == "colhead-fg")])
  }
  
  if(!is.null(col))
    col <- col + 1
  else
    col <- 2:max(l$r)
  
  edit_cell(plot = plot, row = row, col = col, name = name_to_edit, gp = gp)
  
}

# Edit cell
edit_cell <- function(plot, row, col, name="core-fg", ...){
  l <- plot$layout
  ids <- which(l$t %in% row & l$l %in% col & l$name==name)
  for (id in ids){
    newgrob <- editGrob(plot$grobs[id][[1]], ...)
    plot$grobs[id][[1]] <- newgrob
  }
  plot
}

#' @title Create Forest Plot
#' @description This package uses gtable and gridExtra to overlay forest plots.
#' @name forestploter
#' @docType package
#' @author Alimu Dayimu \email{alimdayim@@hotmail.com}
#' @import grid
#' @keywords packagelibrary
#' @seealso \code{\link{grid}},\code{\link{gridExtra}}
#' @importFrom gtable gtable_add_grob gtable_add_rows gtable_add_padding gtable_add_cols
#' @importFrom gridExtra tableGrob ttheme_minimal
#' 
function()
  NULL


# Add vertical line
vert_line <- function(x, gp = grid::gpar(), xlim, is_exp = FALSE){
  
  if(is_exp)
    x <- log(x)
  
  segmentsGrob(x0 = unit(x,"native"),
               x1 = unit(x,"native"),
               y0 = unit(0.01,"npc"),
               y1 = unit(.99,"npc"),
               gp = gp,
               vp = viewport(xscale = xlim))
}


# Get grob corners
getCorners <- function(x) {
  list(xl=grobX(x, 180), xr=grobX(x, 0),
       yb=grobY(x, 270), yt=grobY(x, 90))
}

# Cehck if same length
same_len <- function(...){
  lst <- list(...)
  len <- vapply(lst, length, FUN.VALUE = 1L)
  length(unique(len)) == 1
}


#' Insert text to forest plot
#'
#' This function can be used to insert text to forest plot. Remember to adjust
#' for the row number if you have added text before, including header. This is
#' achieved by inserted new row(s) to the plot and will affect the row number.
#' A text vector can be inserted to multiple columns or rows.
#'
#' @param plot A forest plot object.
#' @param text A character or expression vector, see \code{\link[grid]{textGrob}}.
#' @param row Row to insert the text, this will be ignored if the \code{part} is
#' "header".
#' @param col A numeric value or vector indicating the columns the text will be
#' added. The text will span over the column if a vector is given.
#' @param part Part to insert text, body (default) or header.
#' @param just The justification of the text, \code{"center"} (default),
#' \code{"left"} or \code{"right"}.
#' @param before Indicating the text will be inserted before or after the row.
#' @param gp An object of class \code{"gpar"}, this is the graphical parameter
#'  settings of the text. See \code{\link[grid]{gpar}}.
#' @param padding Padding of the text, default is \code{unit(1, "mm")}
#' 
#' @return A \code{\link[gtable]{gtable}} object.
#'
#' @export
#'
insert_text <- function(plot,
                        text,
                        row = NULL,
                        col = NULL,
                        part = c("body", "header"),
                        just = c("center", "left", "right"),
                        before = TRUE,
                        gp = gpar(),
                        padding = unit(2, "mm")){
  
  if(!inherits(plot, "forestplot"))
    stop("plot must be a forestplot object.")
  
  if(!inherits(gp, "gpar"))
    stop("gp must be a gpar object.")
  
  if(!is.unit(padding))
  padding <- unit(padding, "mm")
  #browser()
  part <- match.arg(part)
  just <- match.arg(just)
  
  # Row must be provided for the body
  if(part == "body" & is.null(row))
    stop("Row must be defined if the text is interting to body.")
  
  # Check text length 
  if(length(text) > 1 && length(row) != length(text) && length(col) != length(text))
    stop("text must have same legnth with row or col.")
  
  # Align text
  tx_x <- switch(just,
                 right = unit(1, "npc") - padding,
                 left  = padding,
                 center = unit(.5, "npc"))
  
  l <- plot$layout
  
  # Header
  if(part == "header"){
    if(!is.null(row))
      row <- row + min(l$b[which(l$name == "colhead-fg")]) - 1
    else
      row <- max(l$b[which(l$name == "colhead-fg")])
  }else{
    row <- max(l$b[which(l$name == "colhead-fg")]) + row
  }
  
  # If the text will be put in columns
  if(!is.null(col) && length(text) == length(col) && length(row) == 1)
    by_col <- TRUE
  else
    by_col <- FALSE
  
  # Span to whole plot if col is missing
  if(is.null(col))
    col <- 2:(ncol(plot) - 1)
  else
    col <- 1 + col # Add 1 to account for padding of the plot
  
  # Order row
  if(!by_col){
    od_rw <- order(row)
    text <- text[od_rw]
    row <- row[od_rw]
  }
  
  nam_text <- paste(ifelse(part == "header", "colhead", "core"), "fg", sep = "-")
  
  for(i in seq_along(row)){
    if(before)
      row[i] <- row[i] - 1 # Account for padding of the plot and header
    
    if(i != 1)
      row[1] <- row[i] + i - 1 # The row number will change after adding one row
    
    if(by_col){
      # Get maximum height of text and add a row
      max_height <- max(convertHeight(stringHeight(text), "mm", valueOnly = TRUE))
      plot <- gtable_add_rows(plot, unit(max_height, "mm") + 2*padding, pos = row[i])
      
      for(j in seq_along(col)){
        txt_grob <- textGrob(label = text[j],
                             gp = gp,
                             x = tx_x,
                             just = just,
                             check.overlap = TRUE,
                             name = "custom-text.insert")
        
        plot <- gtable_add_grob(plot, txt_grob,
                                t = row[i] + 1,
                                b = row[i] + 1,
                                l = col[j],
                                r = col[j],
                                clip = "off",
                                name = nam_text)
      }
      
    }else{
      txt_grob <- textGrob(label = text[i],
                           gp = gp,
                           x = tx_x,
                           just = just,
                           check.overlap = TRUE,
                           name = "custom-text.insert")
      
      plot <- gtable_add_rows(plot,
                              grobHeight(txt_grob) + 2*padding,
                              pos = row[i])
      
      plot <- gtable_add_grob(plot, txt_grob,
                              t = row[i] + 1,
                              b = row[i] + 1,
                              l = min(col),
                              r = max(col),
                              clip = "off",
                              name = nam_text)
    }
    
  }
  
  return(plot)
  
}


#' Create legends
#'
#' This function used to create legends for the forest plot.
#'
#' @param name Character string, Legend name.
#' @param label legend labels (expressions).
#' @param color Colors for the group
#' @param pch Legend symbol.
#' @param lty Line type.
#' @param position Position of the legend, \code{"right"}, \code{"top"},
#' \code{"bottom"}.
#' @param hgap Horizontal gap between the legend entries,
#' see \code{\link[grid]{legendGrob}} for details.
#' @param vgap Vertical gap between the legend entries,
#' see \code{\link[grid]{legendGrob}} for details.
#' @param fontsize Font size of the legend.
#' @param fontfamily Font family of the legend.
#' @param ... Other parameters, not used currently.
#'
#' @return A frame grob
#'
#' @keywords internal
legend_grob <- function(name = "",
                        label,
                        color,
                        position = c("right", "top", "bottom"),
                        hgap = unit(0.1, "lines"), #horizontal gap
                        vgap = unit(0.5, "lines"), #vertical gap
                        pch = 15,
                        lty = 1,
                        fontsize = 12,
                        fontfamily = "",
                        ...
){
  
  position <- match.arg(position)
  
  # Legend title
  title_grob <- textGrob(label = name,
                         just = "left",
                         x = 0,
                         y = 0.5,
                         gp = gpar(fontsize = fontsize,
                                   fontfamily = fontfamily,
                                   fontface = 'bold',
                                   fill = 'black'))
  
  if(position %in% c("top", "bottom")){
    by_row <- FALSE
    ncol <- length(color)
    
  }else{
    by_row <- TRUE
    ncol <- 1
  }
  
  leg_grob <- legendGrob(label, pch = pch, ncol = ncol,
                         do.lines = TRUE, byrow = by_row,
                         hgap = hgap, vgap = vgap,
                         gp = gpar(col = color,
                                   fill = color,
                                   lty = lty,
                                   fontsize = fontsize,
                                   fontfamily = fontfamily))
  
  u0 <- unit(0, "npc")
  u1 <- unit(0.02, "npc")
  
  if(position %in% c("top", "bottom")){
    
    packGrob(frame = packGrob(frameGrob(name = "legend"), title_grob,
                              border = unit.c(u0, u1, u0, u0)),
             grob = leg_grob,
             side = "right")
    
  }else{
    packGrob(frame = packGrob(frameGrob(name = "legend"), title_grob,
                              border = unit.c(u0, u0, u1, u0)),
             grob = leg_grob,
             side = "bottom")
  }
}


#' Make arrow
#'
#' @param x0 Position of vertical line for 0 or 1.
#' @param xlim Limits for the x axis as a vector length 2, i.e. c(low, high)
#' @param arrow_lab Label for the arrow, left and right.
#' @param is_exp If values is exponential.
#' @param gp Graphical parameters for arrow.
#'
#' @keywords internal
make_arrow <- function(x0 = 1, arrow_lab, gp, xlim, is_exp = FALSE){
  
  if(is_exp)
    x0 <- log(x0)
  
  t_lft <- textGrob(arrow_lab[1],
                    x = unit(x0, "native") - unit(0.05, "inches"),
                    y = unit(0.5, "npc"), just = "right",
                    gp = gp,
                    name="arrow.text.left")
  
  t_rgt <- textGrob(arrow_lab[2],
                    x = unit(x0, "native") + unit(0.05, "inches"),
                    y = unit(0.5, "npc"), just = "left",
                    gp = gp,
                    name="arrow.text.right")
  
  t_cord_lft <- getCorners(t_lft)
  t_cord_rgt <- getCorners(t_rgt)
  
  s_lft <- segmentsGrob(t_cord_lft$xl,
                        t_cord_lft$yt + unit(.2, "lines"),
                        t_cord_lft$xr,
                        t_cord_lft$yt + unit(.2, "lines"),
                        gp = gp,
                        arrow = arrow(length=unit(0.05, "inches"),
                                      ends = "first"),
                        name="arrow.left")
  
  s_rgt <- segmentsGrob(t_cord_rgt$xl,
                        t_cord_rgt$yt + unit(.2, "lines"),
                        t_cord_rgt$xr,
                        t_cord_rgt$yt + unit(.2, "lines"),
                        gp = gp,
                        arrow = arrow(length=unit(0.05, "inches"),
                                      ends = "last"),
                        name="arrow.right")
  
  grobTree(gList(t_lft, s_lft, t_rgt, s_rgt),
           vp = viewport(xscale = xlim),
           name = "arrow")
  
}

#' Create xaxis
#'
#' This function used to xais for the forest plot.
#'
#' @param at Numerical vector, create ticks at given values.
#' @param xlab X-axis label.
#' @param is_exp If values is exponential.
#' @param x0 Position of vertical line for 0 or 1.
#' @param ticks_digits Number of digits for the x-axis, default is 1.
#' @param gp Graphical parameters for arrow.
#' @param xlim Limits for the x axis as a vector length 2, i.e. c(low, high)
#'
#' @return A grob
#'
#' @keywords internal
make_xaxis <- function(at, xlab = NULL, x0 = 1, is_exp = FALSE, ticks_digits = 1, gp = gpar(), xlim){
  
  if(is_exp){
    label_at <- log(round(exp(at), ticks_digits))
    x0 <- log(x0)
    labels <- format(round(exp(at), ticks_digits), nsmall = ticks_digits)
  }else {
    label_at <- round(at, ticks_digits)
    labels <- format(round(at, ticks_digits), nsmall = ticks_digits)
  }
  
  maj <- linesGrob(x = unit(c(min(xlim), max(xlim)), "native"),
                   y = unit(c(0.99, 0.99), "npc"),
                   gp = gp,
                   name="major")
  
  maj_cord <- getCorners(maj)
  
  tick <- segmentsGrob(x0 = unit(at, "native"), y0 = maj_cord$yb,
                       x1 = unit(at, "native"), y1 = maj_cord$yb - unit(.5, "lines"),
                       gp = gp,
                       name = "tick")
  
  lab <- textGrob(labels,
                  x = unit(label_at, "native"),
                  y = maj_cord$yb - unit(1, "lines"),
                  gp = gp,
                  # check.overlap=TRUE,
                  name = "label")
  
  if(!is.null(xlab)){
    xlab_gb <- textGrob(xlab,
                        x = unit((xlim[1]+xlim[2])/2, "native"),
                        y = maj_cord$yb - unit(2, "lines"),
                        gp = gp,
                        check.overlap=TRUE,
                        name = "xlab")
    
    grobTree(gList(maj, tick, lab, xlab_gb),
             vp = viewport(xscale = xlim),
             name = "xaxis")
    
  }else{
    grobTree(gList(maj, tick, lab),
             vp = viewport(xscale = xlim),
             name = "xaxis")
  }
  
  
}


#' Create confidence interval
#'
#' @param est Point estimates, numeric
#' @param lower Lower bound
#' @param upper Upper bound
#' @param size Size of the point
#' @param pch Numeric or character vector indicating what sort of plotting
#' symbol to use. See \code{\link[grid]{pointsGrob}}.
#' @param gp Grphical parameters.
#' @param t_height The height confidence interval line end vertices. If 
#' value is `NULL` (default), no vertices will be drawn.
#' @param xlim Limits for the x axis as a vector length 2, i.e. c(low, high)
#' @param nudge_y Offset Y coordinates.
#' @param color Color of the point and the line
#'
#' @keywords internal
makeci <- function(est, lower, upper, pch, size = 1, col = "black", gp = gpar(), 
                   t_height = NULL, xlim = c(0, 1), nudge_y = 0){
  
  rec <- pointsGrob(x = unit(est, "native"),
                    y = 0.5 + nudge_y,
                    pch = pch,
                    size = unit(size, "mm"),
                    gp = gp)
  
  if(upper > max(xlim) | lower < min(xlim)){
    # Both side arrow
    if(upper > max(xlim) & lower < min(xlim)){
      x_pos <- unit(c(0, 1), c("npc", "npc"))
      arrow_side <- "both"
      x_vert <- NULL
    }
    
    # Left side arrow
    else if(lower < min(xlim) & upper <= max(xlim)){
      x_pos <- unit(c(0, upper), c("npc", "native"))
      arrow_side <- "first"
      x_vert <- unit(upper, "native")
    }
    
    # Right side arrow
    else{
      x_pos <- unit(c(lower, 1), c("native", "npc"))
      arrow_side <- "last"
      x_vert <- unit(lower, "native")
    }
    
    lng <- linesGrob(x=x_pos, y = 0.5 + nudge_y,
                     arrow=arrow(length=unit(0.05, "inches"),
                                 ends = arrow_side),
                     gp=gp)
  } else {
    lng <- linesGrob(x=unit(c(lower, upper), "native"), y=0.5 + nudge_y,
                     gp=gp)
    
    x_vert <- unit(c(lower, upper),  "native")
  }
  
  # Draw T end to the CI
  if(!is.null(t_height) & !is.null(x_vert)){
    if(!is.unit(t_height))
      t_height <- unit(t_height, "npc")
    
    vert <- segmentsGrob(x0 = x_vert, y0 = unit(0.5 + nudge_y, "npc") - t_height/2,
                         x1 = x_vert, y1 = unit(0.5 + nudge_y, "npc") + t_height/2,
                         gp = gp,
                         name = "T.end")
  }else {
    vert <- NULL
  }
  
  
  # No dots if outside
  if(est > max(xlim) | est < min(xlim))
    grobTree(gList(lng, vert),
             vp = viewport(xscale = xlim),
             name = "ci")
  else
    grobTree(gList(rec, lng, vert),
             vp = viewport(xscale = xlim),
             name = "ci")
  
}

# Create pooled summary diamond shape
make_summary <- function(est, lower, upper, size = 1, gp, xlim){
  polygonGrob(x = unit(c(lower, est, upper, est), "native"),
              y = unit(0.5 + c(0, 0.5 * size, 0, -0.5*size), "npc"),
              gp = gp,
              vp = viewport(xscale = xlim),
              name = "pooled.diamond")
}



#' Forest plot default theme
#'
#' Default theme for the forest plot, but can pass other parameters. The
#' parameters will be passed to corresponding elements of the forest plot.
#' See \code{\link[grid]{gpar}} for details.
#'
#' @param base_size The size of text
#' @param base_family The font family
#' @param ci_pch Shape of the point estimation. It will be reused if the
#' forest plot is grouped.
#' @param ci_lty Line type of the CI. A vector of line type should be provided
#' for the grouped forest plot.
#' @param ci_lwd Line width of the CI. A vector of line type should be provided
#' for the grouped forest plot.
#' @param ci_Theight A unit specifying the height of the T end of CI. If set to
#' `NULL` (default), no T end will be drawn.
#' @param legend_name Title of the legend.
#' @param legend_position Position of the legend, \code{"right"}, \code{"top"},
#' \code{"bottom"}.
#' @param legend_value Legend labels (expressions). A vector should be provided
#' for the grouped forest plot. A "Group 1" etc will be created if not a vector
#' for a grouped forest plot.
#' @param xaxis_lwd Line width for x-axis.
#' @param xaxis_cex Multiplier applied to font size for x-axis.
#' @param refline_lwd Line width for reference line.
#' @param refline_lty Line type for reference line.
#' @param refline_col Line color for the reference line.
#' @param vertline_lwd Line width for extra vertical line. A vector can be provided
#' for each vertical line, and the values will be recycled if no enough values are
#' given.
#' @param vertline_lty Line type for extra vertical line. Works same as \code{vertline_lwd}.
#' @param vertline_col Line color for the extra vertical line. Works same as \code{vertline_lwd}.
#' @param summary_fill Color for filling the summary diamond shape.
#' @param summary_col Color for borders of the summary diamond shape.
#' @param footnote_cex Multiplier applied to font size for footnote.
#' @param footnote_fontface The font face for footnote.
#' @param footnote_col Color of the footnote.
#' @param title_just The justification of the title, default is \code{'left'}.
#' @param title_cex Multiplier applied to font size for title.
#' @param title_fontface The font face for title, default is \code{'bold'}.
#' @param title_col Color of title.
#' @param title_fontfamily Font family of title. 
#' @param ... Other parameters passed to table. See \code{\link[gridExtra]{tableGrob}}
#'  for details.
#'
#' @importFrom utils modifyList
#'
#' @return A list.
#'
#' @export
#'
forest_theme <- function(base_size = 12,
                         base_family = "",
                         # Confidence interval
                         ci_pch = 15,
                         ci_lty = 1,
                         ci_lwd = 1,
                         ci_Theight = NULL,
                         # Legend
                         legend_name = "Group",
                         legend_position = "right",
                         legend_value = "",
                         # X-axis
                         xaxis_lwd = 0.6,
                         xaxis_cex = 1,
                         # Reference line
                         refline_lwd = 1,
                         refline_lty = "dashed",
                         refline_col = "grey20",
                         # Vertical line
                         vertline_lwd = 1,
                         vertline_lty = "dashed",
                         vertline_col = "grey20",
                         # summary
                         summary_fill = "#4575b4",
                         summary_col = "#4575b4",
                         # Footnote
                         footnote_cex = 0.6,
                         footnote_fontface = "plain",
                         footnote_col = "black",
                         # Title
                         title_just = c("left", "right", "center"),
                         title_cex = 1.2,
                         title_fontface = "bold",
                         title_col = "black",
                         title_fontfamily = base_family,
                         # Legend
                         # legend_lwd = 0.6,
                         ...){
  
  legend_position <- match.arg(legend_position, c("right", "top", "bottom"))
  
  # Default color set
  col_set <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
               "#ffff33","#a65628","#f781bf","#999999")
  
  # Check length
  if(!is.null(ci_Theight) && length(ci_Theight) > 1)
    stop("`ci_Theight` must be of length 1.")
  
  # Recycle if one of the values
  max_len <- list(legend_value, ci_pch, ci_lty, ci_lwd)
  max_len <- max(vapply(max_len, length, FUN.VALUE = 1L), na.rm = TRUE)
  
  if(max_len > 1){
    if(length(legend_value) < max_len)
      stop("legend_value should be provided each groups.")
    
    ci_pch <- rep_len(ci_pch, max_len)
    ci_lty <- rep_len(ci_lty, max_len)
    ci_lwd <- rep_len(ci_lwd, max_len)
    
    
  }
  
  
  # Reference line
  refline_gp <- gpar(lwd = refline_lwd,
                     lty = refline_lty,
                     col = refline_col,
                     fontsize = base_size,
                     fontfamily = base_family)
  
  # Reference line
  vertline_gp <- gpar(lwd = vertline_lwd,
                      lty = vertline_lty,
                      col = vertline_col,
                      fontsize = base_size,
                      fontfamily = base_family)
  
  # Confidence interval
  ci_gp <- list(pch = ci_pch, lty = ci_lty, lwd = ci_lwd, t_height = ci_Theight)
  
  # X-axis
  xaxis_gp <- gpar(lwd = xaxis_lwd,
                   cex = xaxis_cex,
                   fontsize = base_size,
                   fontfamily = base_family)
  
  # Summary
  sum_gp <- gpar(col = summary_col,
                 fill = summary_fill)
  
  # Footnote
  footnote_gp <- gpar(fontsize = base_size,
                      fontfamily = base_family,
                      cex = footnote_cex,
                      fontface = footnote_fontface,
                      col = footnote_col)
  
  # Legend
  legend_gp <- list(fontsize = base_size,
                    fontfamily = base_family,
                    name = legend_name,
                    position = legend_position,
                    label = legend_value)
  
  # Title
  title_gp <- list(just = match.arg(title_just),
                   gp = gpar(cex = title_cex,
                             fontface = title_fontface,
                             col = title_col,
                             fontfamily = title_fontfamily))
  
  # Table body
  core <- list(fg_params = list(hjust = 0,
                                x = 0.05,
                                fontsize = base_size,
                                fontfamily = base_family),
               bg_params = list(fill=c(rep(c("#eff3f2", "white"),
                                           length.out=4))),
               padding = unit(c(4, 3), "mm"))
  
  # Table header
  #https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html
  #https://cran.r-project.org/web/packages/gridExtra/gridExtra.pdf
  colhead <- list(fg_params = list(hjust = 0, x = 0.05,
                                   fontface=2L,
                                   fontsize = base_size,
                                   fontfamily = base_family),
                  bg_params = list(fill = "white"),
                  padding = unit(c(4, 3), "mm"))
  
  default <- list(core = core,
                  colhead = colhead)
  
  tab_theme <- modifyList(default, list(...))
  tab_theme <- modifyList(ttheme_minimal(), tab_theme)
  
  return(list(legend = legend_gp,
              ci = ci_gp,
              xaxis = xaxis_gp,
              footnote = footnote_gp,
              title  = title_gp,
              refline = refline_gp,
              vertline = vertline_gp,
              summary = sum_gp,
              tab_theme  = tab_theme))
  
}


# 
make_group_theme <- function(theme, group_num){
  
  # Default color set
  col_set <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
               "#ffff33","#a65628","#f781bf","#999999")
  
  # If color is given and not have the same length as group number
  if(group_num > 1 & length(theme$ci$col) == 1)
    theme$ci$col <- col_set[1:group_num]
  
  # If line type is given and not have the same length as group number
  if(group_num > 1 & length(theme$ci$lty) == 1)
    theme$ci$lty <- rep_len(theme$ci$lty, group_num)
  
  # If line width is given and not have the same length as group number
  if(group_num > 1 & length(theme$ci$lwd) == 1)
    theme$ci$lwd <- rep_len(theme$ci$lwd, group_num)
  
  # Make legend multiple
  if(group_num > 1 & length(theme$ci$pch) == 1)
    theme$ci$pch <- rep_len(theme$ci$pch, group_num)
  
  if(group_num > 1 && length(theme$legend$label) == 1 && theme$legend$label == ""){
    theme$legend$label <- paste("Group", 1:group_num)
  }
  
  # Check for group and color
  if(group_num > 1 & length(theme$ci$col) < group_num & length(theme$ci$col) > 1)
    stop("More groups than colors.")
  
  # Check for group and legend label
  if(group_num > 1 & length(theme$legend$label) < group_num & length(theme$legend$label) > 1)
    stop("More groups than legend labels.")
  
  return(theme)
  
}

#' Set x-axis ticks
#'
#' Create ticks points.
#'
#' @param at Numerical vector, create ticks at given values.
#' @param is_exp If values is exponential.
#' @inheritParams forest
#'
#' @return A vector
#'
#' @keywords internal

make_ticks <- function(at = NULL,
                       xlim,
                       refline = 1,
                       is_exp = FALSE){
  
  if(is.null(at)){
    
    if(!is_exp){
      ticks_at <- pretty(xlim)
    }else {
      pt_cut <- pretty(range(c(xlim, log(refline))))
      pt_cut <- round(exp(pt_cut), 1) # Keep 1 digits
      ticks_at <- log(unique(pt_cut)) # avoid any duplicate
      # Limit values inside xlim
      ticks_at <- ticks_at[exp(ticks_at) <= max(exp(xlim))]
    }
    
  }else {
    if(is_exp)
      ticks_at <- log(at)
    else
      ticks_at <- at
  }
  
  ticks_at <- ticks_at[is.finite(ticks_at)]
  ticks_at[ticks_at >= min(xlim) & ticks_at <= max(xlim)]
  
}


#' Create xlim
#'
#' Create xlim based on value ranges.
#'
#' @param is_exp If values is exponential.
#' @inheritParams forest
#'
#' @return A list
#'
#' @keywords internal
#'
make_xlim <- function(xlim = NULL,
                      lower,
                      upper,
                      ref_line = ifelse(is_exp, 1, 0),
                      ticks_at = NULL,
                      is_exp = FALSE){
  
  # Use range if missing
  if(is.null(xlim)){
    xlim <- range(c(min(unlist(lower), na.rm = TRUE),
                    ref_line,
                    max(unlist(upper), na.rm = TRUE)),
                  na.rm = TRUE)
  }
  
  if(is_exp){
    if(min(xlim) == 0)
      xlim[which.min(xlim)] <- min(c(unlist(lower),
                                     ref_line,
                                     ticks_at),
                                   na.rm = TRUE)
    xlim <- log(xlim)
    
  }
  
  return(xlim)
  
}




