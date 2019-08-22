# Author: Komal S. Rathi
# Date: 08/19/2019
# Function: Combine TRANSFAC curated and predicted transcription factors
# URL: https://amp.pharm.mssm.edu/Harmonizome/dataset/TRANSFAC+Curated+Transcription+Factor+Targets

# curated
tf.curated <- read.delim('data/TRANSFAC_curated.txt', stringsAsFactors = F, header = F)
tf.curated$label <- 'Curated'
# predicted
tf.predicted <- read.delim('data/TRANSFAC_predicted.txt', stringsAsFactors = F, header = F)
tf.predicted$label <- 'Predicted'
tf <- rbind(tf.curated, tf.predicted)
tf <- tf %>% 
  dplyr::group_by(V1) %>%
  dplyr::summarise(label = toString(label))
write.table(tf, file = 'data/TRANSFAC_TFs.txt', quote = F, sep = "\t", row.names = F, col.names = F)
