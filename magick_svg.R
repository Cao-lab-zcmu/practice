library(magick)

setwd("results/structure_2d")

img <- image_read_svg("3918.svg", width = 3000, height = 3000)

image_display(img)

image_write(img,path="test.svg", format="svg")


