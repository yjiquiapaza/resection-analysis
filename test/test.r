data <- read.csv("resection_log_20260924_153600.csv")
print(data)

data <- read.csv("resection_log_20260924_150107.csv")

barplot(height = data$Azimuth, names.arg = data$Azimuth)
