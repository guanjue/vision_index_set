

d_0 = read.table('1575_peak_rc.txt', header = F)
d_1 = read.table('1576_peak_rc.txt', header = F)
d_2 = read.table('1577_peak_rc.txt', header = F)
d_3 = read.table('1578_peak_rc.txt', header = F)
d_4 = read.table('1579_peak_rc.txt', header = F)
d_5 = read.table('1580_peak_rc.txt', header = F)
d_6 = read.table('1581_peak_rc.txt', header = F)
d_7 = read.table('1582_peak_rc.txt', header = F)
d_8 = read.table('1583_peak_rc.txt', header = F)
d_9 = read.table('1584_peak_rc.txt', header = F)
d_10 = read.table('1585_peak_rc.txt', header = F)
d_11 = read.table('1586_peak_rc.txt', header = F)



plot(density(log2(d_0[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_3[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_6[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_9[,4]+1)), xlim=c(0, 12))

plot(density(log2(d_1[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_2[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_4[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_5[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_7[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_8[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_10[,4]+1)), xlim=c(0, 12))
plot(density(log2(d_11[,4]+1)), xlim=c(0, 12))



d_0r = read.table('1575_random_rc.txt', header = F)
d_1r = read.table('1576_random_rc.txt', header = F)
d_2r = read.table('1577_random_rc.txt', header = F)
d_3r = read.table('1578_random_rc.txt', header = F)
d_4r = read.table('1579_random_rc.txt', header = F)
d_5r = read.table('1580_random_rc.txt', header = F)
d_6r = read.table('1581_random_rc.txt', header = F)
d_7r = read.table('1582_random_rc.txt', header = F)
d_8r = read.table('1583_random_rc.txt', header = F)
d_9r = read.table('1584_random_rc.txt', header = F)
d_10r = read.table('1585_random_rc.txt', header = F)
d_11r = read.table('1586_random_rc.txt', header = F)


plot(density(log2(rbind(d_1, d_1r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='WT')
plot(density(log2(rbind(d_2, d_2r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='WT')
plot(density(log2(rbind(d_4, d_4r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='0hr')
plot(density(log2(rbind(d_5, d_5r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='0hr')
plot(density(log2(rbind(d_7, d_7r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='2hr')
plot(density(log2(rbind(d_8, d_8r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='2hr')
plot(density(log2(rbind(d_10, d_10r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='4hr')
plot(density(log2(rbind(d_11, d_11r[c(1:50000),])[,4]+1)), xlim=c(0, 12), main='4hr')





