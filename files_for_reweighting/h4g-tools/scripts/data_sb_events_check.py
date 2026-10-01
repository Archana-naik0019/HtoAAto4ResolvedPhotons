import numpy as np
import awkward as ak
import glob 

data_files = glob.glob("/eos/home-t/takumar/gitcode/higgsdna_finalfits_tutorial_24/sergi-bdt-files/2022/data/_GeV/nominal/.parquet")
for i in range(10):
    mass = i*5+15
    data = ak.from_parquet(data_files[i])
    data = data[(data.mass_gggg < 115.0) | (data.mass_gggg > 135.0)]
    scores = data.BDT_score
    scores_sorted = np.sort(scores)[::-1]
    print("mass : %d"%mass)
    print("Scores : ",scores_sorted[:10])
    print("Score threshold : %.4f"%scores_sorted[10])
    print("\n")


