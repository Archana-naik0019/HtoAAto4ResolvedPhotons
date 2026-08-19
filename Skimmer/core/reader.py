import uproot
import awkward as ak

class NanoReader:

    def __init__(self, filename):

        self.tree = uproot.open(filename)["Events"]

    def read_scalar(self, branch):

        return self.tree[branch].array(library="ak")

    def read_weight(self, branch):
        return self.tree[branch].array(library="ak")

    def read(self, collection):

        branches = [
            b for b in self.tree.keys()
            if b.startswith(collection + "_")
        ]

        arrays = self.tree.arrays(
            branches,
            library="ak",
        )

        return ak.zip(
            {
                b[len(collection)+1:]: arrays[b]
                for b in branches
            }
        )
