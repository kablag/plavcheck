descr <- trdml$AsTable() %>% filter(target == "CES1", exp.id =="deriv")
fdata <- trdml$GetFData(descr, dp.type = "mdp")
deucl <- diss(fdata[-1], "EUCL")
dist(fdata[-1])
hceucl <- hclust(deucl, "complete")
plot(hceucl)
deucl