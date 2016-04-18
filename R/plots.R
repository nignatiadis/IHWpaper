#' helper function to create nice legends
#'
#' @param gg ggplot2 object
#' @param last_vals, data frame with columns 
#'    label, colour, last_vals (i.e. place label with colour at y-coordinate last_vals)
#' @param xmin Numeric, x axis position at which labels should be placed
#' @param fontsize Integer, fontsize
#'
#' @return Another ggplot2 object
#'
#  @details As described in the following blog post:
#'    http://www.r-bloggers.com/coloring-and-drawing-outside-the-lines-in-ggplot/
#'
#' @examples
#'   labels  <- c("A","B","C")
#'   mypoints <- rbind(data.frame(y=1:3, x=1, label=as.factor(labels)),
#'                      data.frame(y=2:4, x=2, label=as.factor(labels)))
#'   mycolours <- c("#F8766D","#00BA38","#619CFF")
#'   gg <- ggplot(mypoints,aes(x=x,y=y,color=label)) + 
#'                geom_line(size=2) + 
#'                scale_color_manual(values=mycolours) +
#'                xlim(c(0,2.2))
#'   gg
#'   annotation_df <- data.frame(colour=mycolours, last_vals=2:4, label=labels)
#'   pretty_legend(gg, annotation_df, 2.1)
#'
#' @export
#' @importFrom grid gpar textGrob
#' @importFrom ggplot2 annotation_custom ggplot_build theme ggplot_gtable
#' @importFrom cowplot ggdraw

pretty_legend <- function(gg,  last_vals, xmin, fontsize=13){
  #tmp <- ggplot_build(gg)
  #label_name <- tmp$plot$labels$colour
  #label <- levels(droplevels(tmp$plot$data[[label_name]]))

  #df <- tmp$data[[1]]
  #last_vals2 <- group_by(df, colour) %>% summarize(last_vals = y[which.max(x)])
  #last_vals$last_vals <- last_vals$last_vals + offset
  #last_vals$label <- label

  last_vals$colour <- as.character(last_vals$colour)

  gg <- gg + theme(legend.position="none")

  for (i in 1:nrow(last_vals)) {
    gg <- gg + annotation_custom(grob=textGrob(last_vals$label[i], hjust=0,
                                               gp=gpar(fontsize=fontsize, 
                                                       col=last_vals$colour[i])),
                                 xmin=xmin,
                                 xmax=xmin,
                                 ymin=last_vals$last_vals[i],
                                 ymax=last_vals$last_vals[i])
  }
  gg
  gb <- ggplot_build(gg)
  gt <- ggplot_gtable(gb)

  gt$layout$clip[gt$layout$name=="panel"] <- "off"
  gt <- ggdraw(gt)
  gt
}