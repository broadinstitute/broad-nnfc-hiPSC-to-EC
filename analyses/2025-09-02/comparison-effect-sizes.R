apa = read.table("results/2025-09-02/100-D0-r2_power_analysis.tsv", sep="\t", header = T)
bpa = read.table("results/2025-09-02/150-D0-r2_power_analysis.tsv", sep="\t", header = T)

apa$pair = paste(apa$response_id, apa$grna_target, sep = ".")
bpa$pair = paste(bpa$response_id, bpa$grna_target, sep = ".")
merged = merge(apa, bpa, by="pair", suffixes=c(".A", ".B"))

plot(apa$fold_change, bpa$fold_change, main="Positive Controls Fold-Changes",xlab="100-D0-r2", ylab="150-D0-r2")
abline(0,1,col="red")

plot(-10*log(apa$p_value), -10*log(bpa$p_value), main="Positive Controls p-values",xlab="100-D0-r2", ylab="150-D0-r2")
abline(0,1,col="red")
