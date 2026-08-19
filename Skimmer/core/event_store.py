import awkward as ak

class EventStore:

    def __init__(self):
        self.collections = {}
        self.scalars = {}
        self.weights = {}
        self.metadata = {}
        self.temp = {}

    def add_collection(self, name, collection):

        self.collections[name] = collection

    def add_scalar(self, name, value):

        self.scalars[name] = value

    def add_weight(self, name, value):
        self.weights[name] = value

    def add_metadata(self, name, value):
        self.metadata[name] = value
        
    def add_temp(self, name, value):
        self.temp[name] = value