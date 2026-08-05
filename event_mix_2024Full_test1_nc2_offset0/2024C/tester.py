
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

events = NanoEventsFactory.from_root(
    {"cfafd278-00ad-4f0c-b6f0-1fa4ccb94c8a_nc2_offset0.root": "Events"},
    schemaclass=NanoAODSchema,
).events()

print(events.Photon.fields)
