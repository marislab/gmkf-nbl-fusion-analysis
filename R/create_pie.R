# Author: Komal S. Rathi
# Date: 03/9/2019
# Function to create summary pie charts

load('data/clin_mapped.RData')

# Stage
stage <- plyr::count(clin$inss_stage)
stage$x <- paste(stage$x, "\n", stage$freq, sep="")
pie(stage$freq, labels = stage$x, main="Stage")

# Risk
risk <- plyr::count(clin$risk)
risk$x <- paste(risk$x, "\n", risk$freq, sep="")
pie(risk$freq, labels = risk$x, main="Risk")

# MYCN Status
clin$mycn_status <- as.character(clin$mycn_status)
clin$mycn_status[clin$mycn_status == ""] <- "Unknown"
mycn <- plyr::count(clin$mycn_status)
mycn$x <- paste(mycn$x, "\n", mycn$freq, sep="")
pie(mycn$freq, labels = mycn$x, main="MYCN Status")

# Callers
pie(x = c(4719, 542), labels = c("Arriba\n4719", "STAR-Fusion\n542"), main = "Fusion Callers")

