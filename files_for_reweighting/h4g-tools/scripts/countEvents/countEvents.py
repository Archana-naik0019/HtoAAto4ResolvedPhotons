from ROOT import TChain, TFile

print(f"\nBackground samples:")
true_count = 0
files = set()
with open(f"bkg_files.txt", "r") as f:
    for line in f:
        files.add(line)

true_count = len(files)
tree = TChain("Events")
for idx, path in enumerate(files):
    true_count += 1
    try:
        tree.Add(path.replace("\n",""))
        true_count += 1
    except:
        pass

event_counter = tree.GetEntries()
print("Final event count:  ", end="")
print(event_counter)
print(f"Bad Files: {(len(files) - true_count)/len(files) : .2%}")
