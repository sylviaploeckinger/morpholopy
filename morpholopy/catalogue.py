class Galaxy_data:
    def __init__(self,stellar_mass,num_galaxies):
        
        self.stellar_mass = stellar_mass
        self.len_sample = num_galaxies
        self.kappa_co = [ None for i in range(num_galaxies) ]
        self.momentum = [ None for i in range(num_galaxies) ]
        self.axis_ca = [ None for i in range(num_galaxies) ]
        self.axis_cb = [ None for i in range(num_galaxies) ]
        self.axis_ba = [ None for i in range(num_galaxies) ]
        
        self.gas_kappa_co = [ None for i in range(num_galaxies) ]
        self.gas_momentum = [ None for i in range(num_galaxies) ]
        self.gas_axis_ca = [ None for i in range(num_galaxies) ]
        self.gas_axis_cb = [ None for i in range(num_galaxies) ]
        self.gas_axis_ba = [ None for i in range(num_galaxies) ]

    def add_morphology(self,data,i):
        
        self.kappa_co[i] = data[0]
        self.momentum[i] = data[1]
        self.axis_ca[i] = data[2]
        self.axis_cb[i] = data[3]
        self.axis_ba[i] = data[4]

        self.gas_kappa_co[i] = data[5]
        self.gas_momentum[i] = data[6]
        self.gas_axis_ca[i] = data[7]
        self.gas_axis_cb[i] = data[8]
        self.gas_axis_ba[i] = data[9]
