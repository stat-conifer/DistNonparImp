# Figure
# d=5
library(latex2exp)
Res.comb.sg = data.frame(SGM10 = as.vector(L_10_d_5_Res_sg),
                      SGM20 = as.vector(L_20_d_5_Res_sg),
                      SGM50 = as.vector(L_50_d_5_Res_sg),
                      SGM100 = as.vector(L_100_d_5_Res_sg),
                      SGM200 = as.vector(L_200_d_5_Res_sg),
                      SGM500 = as.vector(L_500_d_5_Res_sg))


boxplot(Res.comb.sg, outline = FALSE, cex.axis = 1.8, ylim = c(1.5, 2.5), xaxt = "n", xlab = "", cex.lab = 2)
title(xlab = TeX("$\\mathit{L}$"), cex.lab = 2,cex.axis = 1.8,  font.lab = 3)
axis(1, 1:6, c(10, 20, 50, 100, 200, 500), cex.axis = 1.8)
abline(h=mu, lty = 5, col = "grey5")

Res.comb.ds = data.frame(KDI10 = as.vector(L_10_d_5_Res_ds),
                         KDI20 = as.vector(L_20_d_5_Res_ds),
                         KDI50 = as.vector(L_50_d_5_Res_ds),
                         KDI100 = as.vector(L_100_d_5_Res_ds),
                         KDI200 = as.vector(L_200_d_5_Res_ds),
                         KDI500 = as.vector(L_500_d_5_Res_ds))


boxplot(Res.comb.ds, outline = FALSE, cex.axis = 1.8, ylim = c(1.5, 2.5), xaxt = "n", xlab = "", cex.lab = 2)
axis(1, 1:6, c(10, 20, 50, 100, 200, 500), cex.axis = 1.8)
title(xlab = TeX("$\\mathit{L}$"), cex.lab = 2,cex.axis = 1.8,  font.lab = 3)
abline(h=mu, lty = 5, col = "grey5")

# d=15
Res.comb.sg = data.frame(SGM10 = as.vector(L_10_d_15_Res_sg),
                         SGM20 = as.vector(L_20_d_15_Res_sg),
                         SGM50 = as.vector(L_50_d_15_Res_sg),
                         SGM100 = as.vector(L_100_d_15_Res_sg),
                         SGM200 = as.vector(L_200_d_15_Res_sg),
                         SGM500 = as.vector(L_500_d_15_Res_sg))


boxplot(Res.comb.sg, outline = FALSE, cex.axis = 1.8, ylim = c(1.5, 2.5), xaxt = "n", xlab = "", cex.lab = 2)
title(xlab = TeX("$\\mathit{L}$"), cex.lab = 2,cex.axis = 1.8,  font.lab = 3)
axis(1, 1:6, c(10, 20, 50, 100, 200, 500), cex.axis = 1.8)
abline(h=mu, lty = 5, col = "grey5")

Res.comb.ds = data.frame(KDI10 = as.vector(L_10_d_15_Res_ds),
                         KDI20 = as.vector(L_20_d_15_Res_ds),
                         KDI50 = as.vector(L_50_d_15_Res_ds),
                         KDI100 = as.vector(L_100_d_15_Res_ds),
                         KDI200 = as.vector(L_200_d_15_Res_ds),
                         KDI500 = as.vector(L_500_d_15_Res_ds))


boxplot(Res.comb.ds, outline = FALSE, cex.axis = 1.8, ylim = c(1.5, 2.5), xaxt = "n", xlab = "", cex.lab = 2)
title(xlab = TeX("$\\mathit{L}$"), cex.lab = 2,cex.axis = 1.8,  font.lab = 3)
axis(1, 1:6, c(10, 20, 50, 100, 200, 500), cex.axis = 1.8)
abline(h=mu, lty = 5, col = "grey5")




# Table
# d=5 kernel
order = seq(6,1,-1)
Res.bias.sg = Res.bias.sg[order]
Res.sd.sg = Res.sd.sg[order]
Res.bias.ds = Res.bias.ds[order]
Res.sd.ds = Res.sd.ds[order]

# kernel & sieve
Bias.SD.d5 = cbind(Res.bias.sg,Res.sd.sg,Res.bias.ds,Res.sd.ds)
colnames(Bias.SD.d5)=NULL
print(round(Bias.SD.d5,3))

# d=15 kernel
order = seq(6,1,-1)
Res.bias.sg = Res.bias.sg[order]
Res.sd.sg = Res.sd.sg[order]
Res.bias.ds = Res.bias.ds[order]
Res.sd.ds = Res.sd.ds[order]

# kernel & sieve
Bias.SD.d15 = cbind(Res.bias.sg,Res.sd.sg,Res.bias.ds,Res.sd.ds)
colnames(Bias.SD.d15)=NULL
print(round(Bias.SD.d15,3))

# Table kernel time 
L_d5_ds_time = ds.time
L_d5_all_time = all.time
L_d15_ds_time = ds.time
L_d15_all_time = all.time
NL_d5_ds_time = ds.time
NL_d5_all_time = all.time
NL_d15_ds_time = ds.time
NL_d15_all_time = all.time

ds_kernel_time <- cbind(L_d5_ds_time,L_d15_ds_time,NL_d5_ds_time,NL_d15_ds_time)
print(round(ds_kernel_time,2))
all_ker_time <- c(L_d5_all_time,L_d15_all_time,NL_d5_all_time,NL_d15_all_time)
print(all_ker_time)

# Table sieve time 
load("D:/File/坚果云/喵的坚果云/工作/#3 NonDistributed/Revision/CodeRevise/Setting2-Linear/Sieve/linear_time_d=5.RData")
time <- rowMeans(Res)
all.time <- time[1]
ds.time <- time[2:7]
L_d5_ds_time = ds.time
L_d5_all_time = all.time
load("D:/File/坚果云/喵的坚果云/工作/#3 NonDistributed/Revision/CodeRevise/Setting2-Linear/Sieve/linear_time_d=15.RData")
time <- rowMeans(Res)
all.time <- time[1]
ds.time <- time[2:7]
L_d15_ds_time = ds.time
L_d15_all_time = all.time
load("D:/File/坚果云/喵的坚果云/工作/#3 NonDistributed/Revision/CodeRevise/Setting1-Nonlinear/Sieve/nonlinear_time_d=5.RData")
time <- rowMeans(Res)
all.time <- time[1]
ds.time <- time[2:7]
NL_d5_ds_time = ds.time
NL_d5_all_time = all.time
load("D:/File/坚果云/喵的坚果云/工作/#3 NonDistributed/Revision/CodeRevise/Setting1-Nonlinear/Sieve/nonlinear_time_d=15.RData")
time <- rowMeans(Res)
all.time <- time[1]
ds.time <- time[2:7]
NL_d15_ds_time = ds.time
NL_d15_all_time = all.time

ds_kernel_time <- cbind(L_d5_ds_time,L_d15_ds_time,NL_d5_ds_time,NL_d15_ds_time)
print(round(ds_kernel_time,2))
all_ker_time <- c(L_d5_all_time,L_d15_all_time,NL_d5_all_time,NL_d15_all_time)
print(round(all_ker_time,2))

# dbw plot d=5 & 15
library(latex2exp)
# library(showtext)
ds_rmse <- c()
for(j in 1:length(hc_seq)){
  ds_rmse[j] <- sqrt(mean((abn_rmv(Res.ds[j,])-mu)^2))
}
par(mar=c(4,5,2.5,1),mgp=c(3,1,0))
plot(ds_rmse, ylim = c(0,0.35),lty = 3, pch=15, type = 'o',col="blue",lwd = 5,cex = 2,xaxt = "n",yaxt = "n", xlab = "", ylab = "")
axis(1, at = 1:6, labels = as.character(hc_seq),cex.axis = 1.8)
axis(2, at = c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35), labels = c("0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35"),cex.axis = 1.8)
# showtext_auto(enable = TRUE)
# font_add("CENTURY", "CENTURY.ttf")
# plot(ds_rmse, ylim = c(0,max(ds_rmse)),lty = 1, pch=1, type = 'o',cex = 1,xaxt = "n",yaxt = "n",xlab = "", ylab = "")
title(xlab = TeX("$\\mathit{h_c}$"), cex.lab = 2,cex.axis = 1.8,  font.lab = 3)
rmse_bws = sqrt(mean((abn_rmv(Res.ds[1,])-mu)^2))
abline(h = rmse_bws, lwd = 5, lty = 3, col = "red")
# points(5,rmse_bws,pch = 4, cex = 1, lwd = 1,col = "red")


# dkc plot d=5
library(latex2exp)
kc_seq = c(0.1,0.5,0.9,1.3,1.7,2.1)
ds_rmse <- c()
ds_rmse[1] <- sqrt(mean((abn_rmv(kc_0.1_d_5_Res_ds)-mu)^2))
ds_rmse[2] <- sqrt(mean((abn_rmv(kc_0.5_d_5_Res_ds)-mu)^2))
ds_rmse[3] <- sqrt(mean((abn_rmv(kc_0.9_d_5_Res_ds)-mu)^2))
ds_rmse[4] <- sqrt(mean((abn_rmv(kc_1.3_d_5_Res_ds)-mu)^2))
ds_rmse[5] <- sqrt(mean((abn_rmv(kc_1.7_d_5_Res_ds)-mu)^2))
ds_rmse[6] <- sqrt(mean((abn_rmv(kc_2.1_d_5_Res_ds)-mu)^2))

par(mar=c(4,5,2.5,1),mgp=c(3,1,0))
plot(ds_rmse, ylim = c(0,0.02),lty = 3, pch=15, type = 'o',col="blue",lwd = 5,cex = 2,xaxt = "n",yaxt = "n", xlab = "", ylab = "")
axis(1, at = 1:6, labels = as.character(kc_seq),cex.axis = 1.8)
title(xlab = TeX("$\\mathit{K_c}$"), cex.lab = 2,cex.axis = 1.8,  font.lab = 3)
axis(2, at = c(0.00,0.01,0.02), labels = c("0.00","0.01","0.02"),cex.axis = 1.8)
# axis(1, at = c(1,2), labels = c("1","2"),cex.axis = 1)
rmse_bws = sqrt(mean((abn_rmv(Res.ds[1,])-mu)^2))
abline(h = rmse_bws, lwd = 5, lty = 3, col = "red")
# points(2,rmse_bws,pch = 4, cex = 1, lwd = 1,col = "red")

# d=15
library(latex2exp)
kc_seq = c(0.1,0.5,0.9,1.3,1.7,2.1)
ds_rmse <- c()
ds_rmse[1] <- sqrt(mean((abn_rmv(kc_0.1_d_15_Res_ds)-mu)^2))
ds_rmse[2] <- sqrt(mean((abn_rmv(kc_0.5_d_15_Res_ds)-mu)^2))
ds_rmse[3] <- sqrt(mean((abn_rmv(kc_0.9_d_15_Res_ds)-mu)^2))
ds_rmse[4] <- sqrt(mean((abn_rmv(kc_1.3_d_15_Res_ds)-mu)^2))
ds_rmse[5] <- sqrt(mean((abn_rmv(kc_1.7_d_15_Res_ds)-mu)^2))
ds_rmse[6] <- sqrt(mean((abn_rmv(kc_2.1_d_15_Res_ds)-mu)^2))

par(mar=c(4,5,2.5,1),mgp=c(3,1,0))
plot(ds_rmse, ylim = c(0,0.02),lty = 3, pch=15, type = 'o',col="blue",lwd = 5,cex = 2,xaxt = "n",yaxt = "n", xlab = "", ylab = "")
axis(1, at = 1:6, labels = as.character(kc_seq),cex.axis = 1.8)
title(xlab = TeX("$\\mathit{K_c}$"), cex.lab = 2,cex.axis = 1.8,  font.lab = 3)
axis(2, at = c(0.00,0.01,0.02), labels = c("0.00","0.01","0.02"),cex.axis = 1.8)
# axis(1, at = c(1,2), labels = c("1","2"),cex.axis = 1)
rmse_bws = sqrt(mean((abn_rmv(Res.ds[1,])-mu)^2))
abline(h = rmse_bws, lwd = 5, lty = 3, col = "red")


# Figure real data kernel
library(latex2exp)
library(showtext)
# ds_rmse <- c()
# for(j in 1:length(hc_seq)){
#   ds_rmse[j] <- sqrt(mean((abn_rmv(Res.ds[j,])-mu)^2))
# }
par(mar=c(4,5,2.5,1),mgp=c(3,1,0))
plot(abs(ds_L-whole), ylim = c(0,0.35),lty = 3, pch=15, type = 'o',col="blue",xaxt = "n",yaxt = "n", xlab = "L", ylab = "absolute difference", cex = 2,lwd = 5,font.lab = 3,cex.lab = 1.5)
axis(1, at = 1:6, labels = as.character(L_seq),cex.axis = 1.8)
axis(2, at = c(0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35), labels = c("0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35"),cex.axis = 1.8)
lines(abs(sg_L-whole), lty = 3, col = 3, lwd = 5)
points(abs(sg_L-whole), pch = 17, col = 3, cex = 2)
# showtext_auto(enable = TRUE)
# font_add("CENTURY", "CENTURY.ttf")
# plot(ds_rmse, ylim = c(0,max(ds_rmse)),lty = 1, pch=1, type = 'o',cex = 1,xaxt = "n",yaxt = "n",xlab = "", ylab = "")
abline(h = abs(cc-whole), lwd = 5, lty = 3, col = "red")

# Figure real data sieve
ResRD <- as.data.frame(t(c(all, cc, sg, ds)))
names(ResRD) <- c("all", "cc","sg10", "sg20", "sg50", "sg100", "sg200", "sg500", "ds10", "ds20", "ds50", "ds100", "ds200", "ds500")
cent <- ResRD$all
res <- as.vector(as.matrix(ResRD[, -1]))
z <- abs(unlist(res) - cent)
plot(rep(z[1], 6), type = "l", lty = 3, col = "red", ylim = c(0, 0.35), xaxt = "n", xlab = "L", ylab = "absolute difference", 
     cex = 2,lwd = 5,font.lab = 3,cex.lab = 1.5)
axis(1, at = 1:6, labels = L_seq)
lines(z[2:7], lty = 3, col = 3, lwd = 5)
points(z[2:7], pch = 17, col = 3, cex = 2)
lines(z[8:13], lty = 3, col = "blue", lwd = 5)
points(z[8:13], pch = 15, col = "blue",cex = 2)


