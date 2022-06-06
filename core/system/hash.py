class hash:
    def __init__(self,object, id):
        self.object = object
        self.id = id

    def __eq__(self, other):
        return self.object == other.object and self.id == other.id

    def __hash__(self):
        return hash((self.object, self.id))
