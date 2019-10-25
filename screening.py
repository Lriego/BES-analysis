#!/usr/bin/env python

from sys import argv, exit

def _get_size(value):
	"""

	"""
	max_value = max(int(value[0][2]), int(value[0][3]), int(value[1][2]), int(value[1][3]))
	min_value = min(int(value[0][2]), int(value[0][3]), int(value[1][2]), int(value[1][3]))
	size = max_value - min_value + 1
	return max_value, min_value, size

def screening(inFile):
	"""

	"""
	header2value = {}
	# Creates the dictionary data structure
	IN = open(inFile, "r")
	for line in IN:
		line = line.strip()
		fields = line.split("\t")
		scaffold, BESname, start, end, ori = \
		fields[0], fields[2], fields[3], fields[4], fields[6]
		identifier = BESname.split(".")
		# This way to get the scaffold size is specific for the
		# sequence header format "scaffoldxx|sizexxxxxx".
		# If this is the case, and we are using another reference genome
		# with other headers format, it will crash!
		size = int("".join([s for s in scaffold.split("|")[1] if s.isdigit()]))
		try:
			header2value[identifier[0]].append([scaffold, BESname, start, end, ori, size])
		except Exception:
			header2value[identifier[0]] = [[scaffold, BESname, start, end, ori, size]]
	IN.close()
	# Selects BES accordingly and writes to OUTput
	pairedBES = 0 #Paired-end
	oppositeBES = 0 #Opposite
	singleBES = 0 #Singletons
	positiveBES = 0 #Positive
	negativeBES = 0	#Negative
	unpairedBES = 0 #Unpaired-end
	big_unpair = 0
	small_unpair = 0
	totalBES = 0
	typeOfBES = {
		"pairedBES": [],
		"oppositeBES": [],
		"singleBES": [],
		"positiveBES": [],
		"negativeBES": [],
		"unpairedBES": [],
		"big_unpair": [],
		"small_unpair": [],
		"totalBES": []
	}

	OUT = open("%s_final.gff3" % inFile[:-5], "w+")
	for key, value in header2value.items():
	#	Checking data structure... [OK]
	#	print key, len(value)
	#	print key, value
		if len(value) == 2 and value[0][0] == value[1][0]:
			if (value[0][4] == "+" and value[1][4] == "-" and int(value[0][2]) < int(value[1][2])) or \
			(value[0][4] == "-" and value[1][4] == "+" and int(value[0][2]) > int(value[1][2])):
				pairedBES += 2
				typeOfBES["pairedBES"].append(key)
				end, start, BESsize = _get_size(value)
				scaffold = value[0][0]
				orientation = "."
				color = "#73d216"
				OUT.write("%s\t.\t%s\t%s\t%s\t.\t%s\t.\tID=%s;color=%s;BESsize=%s;BEStype=Paired_BES\n" \
				% (scaffold, key, start, end, orientation, key, color, BESsize))
			if (value[0][4] == "+" and value[1][4] == "-" and int(value[0][2]) > int(value[1][2])) or \
			(value[0][4] == "-" and value[1][4] == "+" and int(value[0][2]) < int(value[1][2])):
				oppositeBES += 2
				typeOfBES["oppositeBES"].append(key)
				end, start, BESsize = _get_size(value)
				scaffold = value[0][0]
				orientation = "."
				color = "#cc0000"
				OUT.write("%s\t.\t%s\t%s\t%s\t.\t%s\t.\tID=%s;color=%s;BESsize=%s;BEStype=Opposite_BES\n" \
				% (scaffold, key, start, end, orientation, key, color, BESsize))
	#			raw_input("Press any key to continue")
			if (value[0][4] == "+" and value[1][4] == "+"):
				positiveBES += 2
				typeOfBES["positiveBES"].append(key)
				end, start, BESsize = _get_size(value)
				scaffold = value[0][0]
				orientation = "+"
				color = "#ff009d"
				OUT.write("%s\t.\t%s\t%s\t%s\t.\t%s\t.\tID=%s;color=%s;BESsize=%s;BEStype=Positive_BES\n" \
				% (scaffold, key, start, end, orientation, key, color, BESsize))

	#			raw_input("Press any key to continue")
			if (value[0][4] == "-" and value[1][4] == "-"):
				negativeBES += 2
				typeOfBES["negativeBES"].append(key)
				end, start, BESsize = _get_size(value)
				scaffold = value[0][0]
				orientation = "-"
				color = "#00ffed"
				OUT.write("%s\t.\t%s\t%s\t%s\t.\t%s\t.\tID=%s;color=%s;BESsize=%s;BEStype=Negative_BES\n" \
				% (scaffold, key, start, end, orientation, key, color, BESsize))

	#			raw_input("Press any key to continue")
		elif len(value) < 2:
			singleBES += 1
			typeOfBES["singleBES"].append(key)
			scaffold = value[0][0]
			BESname = value[0][1]
			start = value[0][2]
			end = value[0][3]
			orientation = value[0][4]
			color = "#878787"
			OUT.write("%s\t.\t%s\t%s\t%s\t.\t%s\t.\tID=%s;color=%s;BEStype=Single_BES\n" \
			% (scaffold, BESname, start, end, orientation, BESname, color))

	#		raw_input("Press any key to continue")
		elif len(value) == 2 and value[0][0] != value[1][0]:
			unpairedBES += 2
			typeOfBES["unpairedBES"].append(key)
			scaffold_0 = value[0][0]
			scaffold_1 = value[1][0]
			BESname_0 = value[0][1]
			BESname_1 = value[1][1]
			start_0 = value[0][2]
			start_1 = value[1][2]
			end_0 = value[0][3]
			end_1 = value[1][3]
			orientation_0 = value[0][4]
			orientation_1 = value[1][4]
			scaffold_size_0 = value[0][5]
			scaffold_size_1 = value[1][5]
			color_0 = None
			color_1 = None
			# Gets the hypothetical BESsize
			if orientation_0 == "+":
				BESremind_size_0 = int(scaffold_size_0) - int(start_0) + 1
				color_0 = "#a89807"
			else:
				BESremind_size_0 = int(end_0) + 1
				color_0 = "#830376"
			if orientation_1 == "+":
				BESremind_size_1 = int(scaffold_size_1) - int(start_1) + 1
				color_1 = "#a89807"
			else:
				BESremind_size_1 = int(end_1) + 1
				color_1 = "#830376"

			# Hypotetical BES size
			hyBESsize = BESremind_size_0 + BESremind_size_1

			if BESremind_size_0 > 100000:
				big_unpair += 1
				flag0 = "flagged"
			else:
				small_unpair += 1
				flag0 = "none"

			if BESremind_size_1 > 100000:
				big_unpair += 1
				flag1 = "flagged"
			else:
				small_unpair += 1
				flag1 = "none"

			OUT.write("%s\t.\t%s\t%s\t%s\t.\t%s\t.\tID=%s;color=%s;Partner=%s;Points=%s;BESsize=%s;BEStype=Unpaired_BES;Flag=%s\n" % \
			(scaffold_0, BESname_0, start_0, end_0, orientation_0, BESname_0, color_0, \
			BESname_1, scaffold_1, abs(hyBESsize), flag0))

			OUT.write("%s\t.\t%s\t%s\t%s\t.\t%s\t.\tID=%s;color=%s;Partner=%s;Points=%s;BESsize=%s;BEStype=Unpaired_BES;Flag=%s\n" % \
			(scaffold_1, BESname_1, start_1, end_1, orientation_1, BESname_1, color_1, \
			BESname_0, scaffold_0, abs(hyBESsize), flag1))
		else:
			pass

	OUT.close()
	# Writes summary to file
	totalBES = pairedBES + oppositeBES + positiveBES + negativeBES + unpairedBES + singleBES
	SUMMARY = open("%s.summary" % inFile[:-5], "w+")
	#		raw_input("Press any key to continue!")
	SUMMARY.write("BES type\t%s\n" % inFile[:-5])
	SUMMARY.write("Paired-end\t%d\n" % pairedBES)
	SUMMARY.write("Opposite\t%d\n" % oppositeBES)
	SUMMARY.write("Positive\t%d\n" % positiveBES)
	SUMMARY.write("Negative\t%d\n" % negativeBES)
	SUMMARY.write("Unpaired-end\t%d\n" % unpairedBES)
	SUMMARY.write("  >100 Kb\t%d\n" % big_unpair)
	SUMMARY.write("  <100 Kb\t%d\n" % small_unpair)
	SUMMARY.write("Singletons\t%d\n" % singleBES)
	SUMMARY.write("BES Total\t%d\n" % totalBES)

	SUMMARY.close()

def main():
	"""
		Main function
	"""
	inFile = argv[1]
	screening(inFile)

if __name__ == "__main__": main()
