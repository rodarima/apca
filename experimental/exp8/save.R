
#dat = read.csv('data/exp8a.csv', head=FALSE)
#hist(dat, 100, prob=TRUE)

# Read all runned configuration for steps
steps = read.csv('data/8c-steps.csv', head=FALSE)[[1]]
runs = 100000
bits = 10

for(s in steps)
{
	csv_file = sprintf("data/8c-%d-%d-%d.csv", runs, s, bits)
	rda_file = sprintf("data/8c-%d-%d-%d.rda", runs, s, bits)
	print(csv_file)
	csv = read.csv(csv_file, sep=' ', head=TRUE)
	save(csv, file=rda_file)
}
