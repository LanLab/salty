library(UpSetR)
df1 = read.csv('')
upset(df1, nsets = 16, point.size = 2,order.by = "freq", text.scale = 1, mb.ratio = c(0.5, 0.5),line.size = 0.5)