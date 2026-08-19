import awkward as ak

def event_mask(collections):
    """
    collections is a dictionary containing the selected objects.
    """

    #jets = collections["Jet"]
    electrons = collections["Electron"]
    #muons = collections["Muon"]
    photons = collections["Photon"]

    #nJet = ak.num(jets)
    nEle = ak.num(electrons)
    #nMuon = ak.num(muons)
    nPho = ak.num(photons)

    mask = (
        #(nJet >= 1)
        #&
        (nPho >= 4)
        #&
        #((nEle + nMuon) >= 1)
    )

    return mask