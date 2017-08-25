import docopt
import processNX

if __name__ == "__main__":
	arguments = docopt.docopt("\n".join(open("nextera.txt").readlines()), version="0.1")
	if arguments["trim"]:
		processNX.trim(arguments)
	if arguments["align"]:
		processNX.align(arguments)
	if arguments["coverage"]:
		processNX.coverage(arguments)
	
