
################################
## ORDERED PREDICTOR FUNCTIONS #
################################

library(ggthemes)
library(GGally)

my_lower <- function(data, mapping, ...) {
  
  # get the x and y data to use the other code
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # compute the correlations
  corr <- cor(x, y, method = "p", use = "pairwise")
  abs_corr <- abs(corr)
  
  # plot the cor value
  ggally_text(
    label = formatC(corr, digits = 2, format = "f") %>% str_replace(., "0.", "."),
    mapping = aes(),
    size = 4,
    color = canva_pal("Green fields")(4)[2]) +
    scale_x_continuous(NULL, breaks = NULL) +
    scale_y_continuous(NULL, breaks = NULL)
}

my_diag <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) + 
    geom_density(fill = canva_pal("Green fields")(4)[1], linewidth = 0) +
    scale_x_continuous(NULL, breaks = NULL) +
    scale_y_continuous(NULL, breaks = NULL)
}

my_upper <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) + 
    geom_hex(bins = 18) +
    scale_fill_gradient(low = canva_pal("Green fields")(4)[4],
                        high = canva_pal("Green fields")(4)[3]) +
    scale_x_continuous(NULL, breaks = NULL) +
    scale_y_continuous(NULL, breaks = NULL) +
    theme(panel.background = element_rect(fill = canva_pal("Green fields")(4)[2]))
}