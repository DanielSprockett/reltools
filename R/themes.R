#' @title A custom theme for the Relman Lab
#'
#' @description This theme is inspired by the unique aesthetic sensibilities of David Relman, M.D.
#'
#' @param NULL
#'
#' @examples
#'
#' Make a Plot
#' p <- ggplot(data.frame(A =rnorm(100), B = rnorm(100)))
#' p <- p + geom_point(aes(x = A, y = B), color = "#FF56D6", size =3)
#' p <- p + ggtitle("A plot")
#' p <- p + theme_relman()
#' p
#'
theme_relman <- function () {
  if (!requireNamespace("extrafont", quietly = TRUE)) {
    stop("The R package extrafont is needed for this function to work.
         Please install it and run the command font_import().
         This may take a few minutes",
         call. = FALSE)
  }

  theme_dark(base_size=12, base_family="")
  theme(
    text=element_text(family = "Comic Sans MS", face = "plain",
                      color = "black", size = 16,
                      hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                      margin = margin(), debug = FALSE),

    plot.background = element_rect(fill = "black", colour = NA),
    plot.title = element_text(size = rel(1.2), colour = "#01FFFF"),
    panel.background = element_rect(fill = "black", colour = NA),
    panel.grid.major = element_line(color = "white", size = 0.25),
    panel.grid.minor = element_line(color = "white", size = 0.25),
    panel.border = element_rect(colour = "white", fill=NA, size = 0.25),
    axis.ticks = element_line(colour = "white", size = 0.25),
    axis.text = element_text(size = rel(1.2), colour = "#CCFFCC"),
    axis.title.x = element_text(size = rel(1.2), colour = "#00FF02"),
    axis.title.y = element_text(size = rel(1.2), colour = "#00FF02", angle=90),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(size = rel(1.2), colour = "#00FF02"),
    legend.title = element_text(hjust = 0, colour = "#01FFFF"),
    strip.background = element_rect(fill = "black", colour = NA),
    strip.text = element_text(colour = "white", size = rel(0.8)),
    complete = TRUE)
}
