dat = read.csv('firsthh.csv', sep='\t', head=FALSE)
pdf('err.pdf', width=8, height=6)
#par(mar=c(4,4,0.5,0.5))
plot(dat[,1], log2(dat[,2]), pch='.',
	xlab='Bit-width', ylab='log2(err)', main='Error in the result')
dev.off()

# Compute histogram for 100 bits
dat = read.csv('firsthh-10000.csv', sep='\t', head=FALSE)
err100 = dat[dat[1] == 100,2]
logerr = log2(err100)


pdf('err100.pdf', width=8, height=6)
#par(mar=c(4,4,0.5,0.5))
hist(logerr, 100, prob=TRUE,
	xlab='log2(err)', main='Dispersion of log2(err) for 100 bits')
dev.off()

err50 = dat[dat[1] == 50,2]
logerr = log2(err50)

pdf('err50.pdf', width=8, height=6)
#par(mar=c(4,4,0.5,0.5))
hist(logerr, 100, prob=TRUE,
	xlab='log2(err)', main='Dispersion of log2(err) for 50 bits')
dev.off()

err10 = dat[dat[1] == 10,2]
logerr = log2(err10)

pdf('err10.pdf', width=8, height=6)
#par(mar=c(4,4,0.5,0.5))
hist(logerr, 100, prob=TRUE,
	xlab='log2(err)', main='Dispersion of log2(err) for 10 bits')
dev.off()
