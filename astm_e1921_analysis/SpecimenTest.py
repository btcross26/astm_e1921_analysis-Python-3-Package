from .helper import in_to_mm, lbf_to_newtons
import matplotlib.pyplot as plt
import numpy as np
   
class SpecimenTest(object):
    '''
    SpecimenTest object - combines a specimen with a material and test data
    '''
    
    class TestResults(object):
        def __init__(self, Jc = None, KJc = None, KJc_valid_1T = None):
            self.Jc = Jc
            self.KJc = KJc
            self.KJc_valid_1T = KJc_valid_1T
                    
    def __init__(self, specimen, material, temperature, data, elastic_ul):
        '''
        Constructor arguments:
            specimen: Specimen object
            material: Material object
            temperature: temperature of test run
            data: dictionary with keys "cod", "disp", and "force" containing equal length numpy float arrays
            elastic_ul: upper limit of elastic behavior in force units to determine C0 for analysis
        '''
        self.specimen = specimen
        self.material = material
        self.temperature = temperature
        assert len(data["cod"]) == len(data["force"]) and len(data["cod"]) == len(data["disp"])
        self.data = data
        self.elastic_ul = elastic_ul
        self.results = SpecimenTest.TestResults()
        self._analyze()
        
    def __str__(self):
        return "<SpecimenTest at %s>" % (hex(id(self)),)
        
    def __repr__(test):
        return self.__str__()
               
    def _analyze(self):
        '''
        KJc calculations on test data
        '''
        # create and initialize test variables for cleaner code (convert lengths to m)
        P = np.array(self.data["force"])
        cod = np.array(self.data["cod"]) / 1000.0
        a0_over_W = self.specimen.a0 / self.specimen.W
        a_over_W = (self.specimen.a0 + cod) / self.specimen.W
        S = self.specimen.S / 1000.0
        B = self.specimen.B / 1000.0
        BN = self.specimen.BN / 1000.0
        W = self.specimen.W / 1000.0
        b0 = (self.specimen.W - self.specimen.a0) / 1000.0
        
        # calculate cod area <A> using trapezoidal rule
        n = len(P)        
        force_int_avg = 0.5 * P[:n - 1] + 0.5 * P[1:]
        cod_diff = np.diff(cod)
        A = np.concatenate([np.array([0]), np.cumsum(cod_diff * force_int_avg)])
        
        # calculate Ke (convert Pa root m to MPa root m in final Ke calculation)
        k1 = P * S / ((B * BN)**(0.5) * W**(1.5))
        f1 = 1.5 * a0_over_W**(0.5) / (1 + 2 * a0_over_W)
        f2 = 1.99 - a0_over_W * (1.0 - a0_over_W) * (2.15 - 3.93 * a0_over_W + 2.7 * a0_over_W**2)
        f3 = (1.0 - a0_over_W)**(1.5)
        Ke = 1e-6 * k1 * f1 * f2 / f3
        
        # calculate Je
        Je = (1.0 - self.material.poisson_ratio**2) * Ke**2 / self.material.E
        
        # calculate Jp (convert final Jp from 
        trimmed_index = P <= self.elastic_ul
        C0 = np.var(cod[trimmed_index]) / np.cov(cod[trimmed_index], P[trimmed_index])[0, 1]
        nu = 3.667 - 2.199 * a_over_W + 0.4376 * a_over_W**2
        Ap = A - 0.5 * C0 * P**2
        Ap[:max(np.arange(len(Ap))[Ap < 0]) + 1] = 0.0
        Jp = 1e-6 * nu * Ap / (BN * b0)
        
        # calculate Jc (units MPa * m)
        self.results.Jc = (Je + Jp)[-1]
        
        # calculate KJc
        self.results.KJc = np.sqrt(self.results.Jc * self.material.E / (1.0 - self.material.poisson_ratio**2))
        
        # calculate KJc validity limit
        self.results.KJc_valid_1T = np.sqrt(self.material.E * b0 * self.material.fy / (30 * (1 - self.material.poisson_ratio**2)))

    def update(self, specimen = None, material = None, data = None, elastic_ul = None):
        '''
        Update test results based on a change in Specimen Test parameters.  Arguments
        set to equal None will remain unchanged in the updated results analysis.
        
        Arguments:
            specimen: new test specimen object
            material: new material object
            data: new data
            elastic_ul: process data using new elastic_ul for determining initial slope
        '''
        if isinstance(specimen, Specimen):
            self.specimen = specimen
        if isinstance(material, Material):
            self.material = material
        if data:
            assert len(data["cod"]) == len(data["force"]) and len(data["cod"]) == len("disp")
            self.data = data
        if elastic_ul:
            self.elastic_ul = elastic_ul
        self._analyze()
        
    def plot_data(self, grid = True, show = True, english_units = False):
        '''
        Plot specimen test data
        
        Arguments:
            grid: include xy grid on plot
            show: show plot in window
            english_units: plot in English units in. and lbf (default is SI -- mm and N)
            
        Returns:
            A tuple containing the plot figure object and the plot axes object
        '''
        if english_units == True:
            funits, lunits = "lbf", "in."
            c, f = self.data["cod"] / in_to_mm(1.0), self.data["force"] / lbf_to_newtons(1)
        else:
            funits, lunits = "N", "mm"
            c, f = self.data["cod"], self.data["force"]
        fig = plt.figure(self.specimen.id + " Test Data")
        axes = fig.add_subplot(111)
        axes.plot(c, f, "r-")
        plt.autoscale()
        plt.xlim([0.0, plt.xlim()[1]])
        plt.xlabel("Crack Opening Displacement (%s)" % (lunits,))
        plt.ylim([0.0, plt.ylim()[1]])
        plt.ylabel("Load (%s)" % (funits,))
        plt.grid(grid)
        if show:
            plt.show()
        return fig, axes