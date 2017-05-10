
#dat = read.csv('data/exp8a.csv', head=FALSE)
#hist(dat, 100, prob=TRUE)

# Read all runned configuration for steps
steps = read.csv('data/8c-steps.csv', head=FALSE)[[1]]
runs = 100000
bits = 10

for(s in steps)
{
	#csv_file = sprintf("data/8c-%d-%d-%d.csv", runs, s, bits)
	#print(csv_file)
	#csv = read.csv(csv_file, sep=' ', head=TRUE)

	rda_file = sprintf("data/8c-%d-%d-%d.rda", runs, s, bits)
	print(rda_file)
	load(rda_file) # loads 'csv'

	#print(summary(csv))
	sha = shapiro.test(sample(csv$err, 5000))

	img_file = sprintf("img/8c-gsum-err-%d-%d-%d.png", runs, s, bits)
	png(img_file, width=640*1.5,height=480*1.5,units="px")
	title = sprintf("Error after %d sums", s)
	plot(csv$gsum, csv$err, pch='.', main=title, xlab='Sum', ylab='Error')
	dev.off()

	img_file = sprintf("img/8c-hist-err-%d-%d-%d.png", runs, s, bits)
	png(img_file, width=640*1.5,height=480*1.5,units="px")
	title = sprintf("Histogram of the error after %d sum. Shapiro p-value = %g",
		s, sha$p.value)
	hist(csv$err, 100, main=title, xlab='Error')
	dev.off()

	img_file = sprintf("img/8c-qqnorm-err-%d-%d-%d.png", runs, s, bits)
	png(img_file, width=640*1.5,height=480*1.5,units="px")
	title = sprintf("Normal Q-Q plot with %d steps. Shapiro p-value = %g",
		s, sha$p.value)
	qqnorm(csv$err, main=title)
	qqline(csv$err, col='red')
	dev.off()

	rm(csv)
}
