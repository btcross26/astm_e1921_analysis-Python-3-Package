class Specimen(object):
    '''
    Class describes test specimen geometry and material for SEB specimen.
    '''
    
    def __init__(self, id, W, B, BN, a0, S):
        '''
        Constructor Arguments
        ---------------------
        id: str
            ID of specimen.
        W: float
            Specimen depth (length).
        B: float
            Specimen width at gross section (length).
        BN: float
            Net specimen width at side-grooves, if applicable (length).
        a0: float
            Initial crack length (float).
        S: float
            Span length of SEB specimen (float).
            
        '''
        self.id = id
        self.W = W
        self.B = B
        self.BN = BN
        self.a0 = a0
        self.S = S
        
    def __str__(self):
        string = ""
        string += "Specimen ID: %s\n" % (self.id,)
        string += "W  = %.4f mm\n" % (self.W,)
        string += "B  = %.4f mm\n" % (self.B,)
        string += "BN = %.4f mm\n" % (self.BN,)
        string += "a0 = %.3f mm\n" % (self.a0,)
        string += "S  = %.3f mm\n" % (self.S,)
        return string
        
    def __repr__(self):
        return "<Specimen at %s>" % (hex(id(self)),)