import docopt
import process
import analysis
import PhiXpatch
import comparison

if __name__ == "__main__":
	arguments = docopt.docopt("\n".join(open("diffseq.txt").readlines()), version="0.1")
	if arguments["merge"]:
		process.merge(arguments)
	if arguments["dedup"]:
		process.dedup(arguments)
	if arguments["trimadaptors"]:
		process.trimadaptors(arguments)
	if arguments["filter"]:
		process.filter(arguments)
	if arguments["align"]:
		process.align(arguments)
	if arguments["coverage"]:
		process.coverage(arguments)
	if arguments["coverageplots"]:
		process.coverageplots(arguments)
	if arguments["stacked"]:
		analysis.stacked(arguments)
	if arguments["heatmap"]:
		analysis.heatmap(arguments)
	if arguments["diffsignal"]:
		analysis.diffsignal(arguments)
	if arguments["model"]:
		analysis.model(arguments)
	if arguments["removePhiX"]:
		PhiXpatch.preprocess(arguments)
	if arguments["prepareDS"]:
		comparison.prepareDS(arguments)
	if arguments["prepareNX"]:
		comparison.prepareNX(arguments)
	if arguments["dilutioncomparison"]:
		comparison.dilutioncomparison(arguments)
	if arguments["posfreq"]:
		analysis.posfreq(arguments)
