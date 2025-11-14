custom_theme <- function() {
  theme_bw() +
    theme(
      # Customize other theme elements here
      axis.text=element_text(size = 12),
      axis.text.x=element_text(angle=90),
      axis.title=element_text(size = 14, face = "bold"),
      legend.title=element_text(size = 14, face = "bold"),
      legend.text=element_text(size=8),
      strip.background = element_rect(fill = "white", colour = "black")
    )
}

large_font <- function() {
  theme_bw() +
    theme(
      # Customize other theme elements here
      axis.text=element_text(size = 14),
      axis.text.x=element_text(angle=90),
      axis.title=element_text(size = 18, face = "bold"),
      legend.title=element_text(size = 14, face = "bold"),
      legend.text=element_text(size=12),
      strip.background = element_rect(fill = "white", colour = "black"),
      axis.title.x = element_blank()
    )
}
