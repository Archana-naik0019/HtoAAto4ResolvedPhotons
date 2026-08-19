def photon_mask(photons):

    return (
        (photons.pt > 12)
        &
        (abs(photons.eta) < 2.5)
    )
