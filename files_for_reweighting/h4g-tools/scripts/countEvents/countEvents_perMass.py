from ROOT import RDataFrame

for mass in range(15,65,5):
	print(f"\n{mass} GeV sample:")
	true_count, false_count = (0,0)
	event_counter = 0
	files = set()
	with open(f"files_{mass}GeV.txt", "r") as f:
		for line in f:
			files.add(line)

	for idx, path in enumerate(files):
		true_count += 1
		try:
			f = RDataFrame("Events", path.replace("\n",""))
			event_counter += f.Count().GetValue()
			false_count += 1
			#if idx % 20 == 0:
				#print(event_counter)
		except:
			pass

	print("Final event count:  ", end="")
	print(event_counter)
	print(f"Bad Files: {(true_count-false_count)/true_count:.2%}")
