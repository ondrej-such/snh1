source("lda.R")

# df1 <- bc_dataset("mnist")
# df2 <- bc_dataset("satimage")

df <- map(files, function(f) bc_dataset(f)) |> list_rbind()

#write.csv(rbind(df1, df2), file = "data/exp4.csv")
write.csv(df, file = "data/exp4.csv")
