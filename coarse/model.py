class CarPhysicsModel:
    def __init__(self, array):
        self.array = array

    def __call__(self):
        self.set_simple_parameters()

    def __getitem__(self, key):
        return self.array[key]

    def set_simple_parameters(self):
        self['fuel cell system efficiency'] = (
            self['fuel cell stack efficiency'] / self['fuel cell own consumption']
        )
