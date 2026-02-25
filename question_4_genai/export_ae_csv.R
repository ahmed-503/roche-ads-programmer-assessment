# In R / RStudio
library(pharmaversesdtm)

# Create the Q4 folder if needed
dir.create("question_4_genai", showWarnings = FALSE)

# Export AE dataset to CSV
write.csv(pharmaversesdtm::ae, "question_4_genai/adae.csv", row.names = FALSE)