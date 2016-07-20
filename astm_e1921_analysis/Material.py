class Material(object):
    '''
    Class describes steel material properties.
    '''
    
    def __init__(self, description, fu, fy, E, poisson_ratio = 0.3):
        '''
        Constructor Arguments
        ---------------------
        description: str
            Description of steel material.            
        fu: float
            Ultimate strength of steel (MPa).            
        fy: float
            Static yield strength of steel (MPa).          
        E: float
            Elastic modulus of steel (Mpa).           
        poisson_ratio: float
            Poisson's ratio of steel material.
        '''
        self.description = description
        self.fu = fu
        self.fy = fy
        self.E = E
        self.poisson_ratio = poisson_ratio
        
    def __str__(self):
        string = ""
        string += "Material description: %s\n" % (self.description,)
        string += "fu = %.1f MPa\n" % (self.fu,)
        string += "fy = %.1f MPa\n" % (self.fy,)
        string += "E  = %d MPa\n" % (self.E,)
        string += "PR = %s (Poisson's ratio)\n" % (str(self.poisson_ratio),)
        return string
        
    def __repr__(self):
        return "<Material at %s>" % (hex(id(self)),)