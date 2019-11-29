Rain <- list()
Rain$rain <- as.matrix(read.csv("Rain_rain.csv", header = FALSE)) + 1
Rain$months <- as.matrix(read.csv("Rain_months.csv", header = FALSE))
save(Rain, file="Rain.RData")
