from .helper import stress_intensity_in_SI_to_english, fahrenheit_to_celsius
from .Material import Material
from .Specimen import Specimen
from .SpecimenTest import SpecimenTest
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from scipy.stats import norm
from scipy.integrate import simps
from scipy.optimize import minimize
    
class Sample(object):
    '''
    Sample object - combines a series of Specimen object/test data/
    temperature combinations and provides analysis tools to calculate Master
    curve reference temperature via the standard ASTM E1921 methods as well
    as methods as provided in the inhomogeneous steel paper by Kim Wallin
    titled "Master Curve analysis of inhomogeneous ferritic steels" (2003).
    '''
    
    class SampleResults(object):
        def __init__(self, num_valid = None, rmat = None, T0Q = None, sigma_T0Q = None, \
                     req = None, KJc_median = None, K0 = None, sigma = None, tol = None):
        '''
        Constructor Arguments
        ---------------------        
        num_valid: int
            Number of valid SpecimenTest objects in the sample.            
        rmat: numpy.matrix 
            Matrix with columns temperature, KJc, and valid (1 if valid, 0 if
            not valid). Non-tests should not be added to the sample to begin
            with and should be screened prior to buiding the Sample object.            
        T0Q: float
            Reference temperature as per ASTM E1921 calculations.            
        sigma_T0Q: float
            Reference temperature standard deviation as per ASTM E1921
            calculations.            
        req: float
            Calculation for valid sample as per ASTM E1921 calculations --
            a value greater than or equal to 1.0 equates to a valid sample.            
        KJc_median: callable
            Function of T that calculates median KJc value based on reference
            temperature.            
        sigma: callable
            Function of T that calculates standard deviation.            
        tol: callable
            Function of T, alpha that calculates predication interval line at
            T for given alpha, e.g. 97.5 percentile prediction interval
            percentile at -70 celsius.
        '''
            self.num_valid = num_valid
            self.rmat = rmat
            self.T0Q = T0Q
            self.sigma_T0Q = sigma_T0Q
            self.req = req
            self.KJc_median = KJc_median
            self.K0 = K0
            self.sigma = sigma
            self.tol = tol 
    
    def __init__(self, description, material):
        '''
        Constructor arguments
        ---------------------
        description: str
            Sample description, e.g. "Steel Heat #1 LT Orientation".           
        material: Material
            Object specifying test material of the specimens.
        '''
        self.description = description
        self.material = material
        self.tests = dict()
        self.results = None
        
    def __str__(self):
        string = ""
        string += "Sample Description: %s\n" % (self.description,)
        string += "Material: %s" % (self.material.description,)
        return string
        
    def __repr__(self):
        return "<Sample at %s>" % (hex(id(self)),)
        
    def __getitem__(self, id):
        return self.tests[id]
        
    def __iter__(self):
        for value in self.tests.values():
            yield value
            
    def __len__(self):
        return len(self.tests)
            
    def specimen_id_list(self):
        '''
        Returns
        -------
        A list of specimen ids currently contained in the sample.
        '''
        return [key for key in self.tests.keys()]
        
    def add_test(self, specimen, temperature, data, elastic_ul):
        '''
        Adds a SpecimenTest object to the sample.
        
        Arguments
        ---------
        specimen: Specimen
            A Specimen object.           
        temperature: float
            The temperature at which the test was conducted.       
        data: dict
            Test data dictionary in the format required by the SpecimenTest
            object.    
        elastic_ul: float
            Upper load limit in N to be used in calculating initial elastic
            slope.
        '''
        self.tests[specimen.id] = SpecimenTest(specimen, self.material, temperature, \
                                               data = data, elastic_ul = elastic_ul)
        
    def remove_test(self, specimen_id):
        '''
        Remove a SpecimenTest object from the sample by specimen id.  If
        specimen id is not in the sample, nothing happens.
        
        Arguments
        ---------
        specimen_id: str
            Specimen id, e.g. "2LS-2".
        '''
        if specimen_id in self.tests:
            del self.tests[specimen_id]
            
    def change_material(self, material = None):
        '''
        Changes material object of all specimens in sample.  Analysis must be
        re-run manually recalculate the master curve.
        
        Arguments
        ---------
        material: Material
            New material object.
        '''
        if isinstance(material, Material):
            self.material = material
            for st in self.tests:
                self.tests[st].update(material = material)

    def analyze(self, ignore_censored = False):
        '''
        Runs standard Master Curve analysis as specified in ASTM E1921.
        Calculated results are stored in instance 'results' object.
        '''
        self.results = Sample.SampleResults()
        results = np.zeros((len(self.tests), 3))
        for index, st in enumerate(self.tests.values()):
            results[index, 0] = st.temperature
            if st.results.KJc <= st.results.KJc_valid_1T:
                results[index, 1] = adjust_KJc(st.results.KJc, st.specimen.B, 25.4)
                results[index, 2] = 1
            else:
                results[index, 1] = adjust_KJc(st.results.KJc_valid_1T, st.specimen.B, 25.4)
                results[index, 2] = 0
        if ignore_censored == True:
            results = results[results[:,2] == 1,:]
        T0Q = fsolve(self._create_ofunc(results), np.max(results[:,0]), full_output = False)[0]
        self.results.num_valid = np.sum(results[:,2])
        self.results.rmat = results
        self.results.T0Q = T0Q
        self.results.sigma_T0Q = self._T0_sd(results, T0Q)
        self.results.req = np.sum(np.array([0, 8**-1, 7**-1, 6**-1, 0])[np.digitize(results[:,0] - T0Q * \
                                  results[:,2], np.array([-50, -35, -14, 50]))])
        self.results.KJc_median = lambda T: 30.0 + 70.0 * np.exp(0.019 * (T - self.results.T0Q))
        self.results.K0 = lambda T: (self.results.KJc_median(T) - 20.0) / np.log(2.0)**(0.25) + 20.0
        self.results.sigma = lambda T: 0.28 * self.results.KJc_median(T) * (1.0 - 20.0 / self.results.KJc_median(T))
        self.results.tol = lambda T, alpha: 20.0 + np.log(1.0 / (1.0 - alpha))**(0.25) * (11.0 + 77.0 * np.exp(0.019 * (T - self.results.T0Q)))
        
    def _create_ofunc(self, results):
        '''
        Helper function that returns objective function to be solved for as
        per ASTM E1921 calculations.  Not meant to be called externally.
        '''
        temp, KJc, delta = results[:,0], results[:,1], results[:,2]
        def ofunc(T0):
            C1 = np.exp(0.019 * (temp - T0))
            C2 = 11.0 + 76.7 * C1
            A = np.sum(delta * C1 / C2)
            B = np.sum((KJc - 20.0)**4 * C1 / C2**5)
            return A - B
        return ofunc
        
    def _T0_sd(self, results, T0):
        '''
        Helper function that calculates reference temperature standard
        deviation as per ASTM E1921 calculations. sigma_exp is assumed 4.0 as
        recommended in the specification. Not meant to be called externally.
        '''
        KJc_median_eq = np.mean((30.0 + 70.0 + np.exp(0.019 * (results[:,0] - T0)))[results[:,2] == 1])
        if KJc_median_eq >= 83:
            Beta = 18.0
        elif KJc_median_eq >= 66:
            Beta = 18.8
        elif KJc_median_eq >= 58:
            Beta = 20.1
        else:
            return "Invalid Beta"
        return np.sqrt(4.0**2 + Beta**2 / self.results.num_valid)
                
    def MML_estimation(self, optimize_kwargs = None, ignore_censored = False):
        '''
        MML estimation as per paper by Kim Wallin through discussions with
        Scibetta.
        
        Arguments
        ---------
        optimize_kwargs: dict
            Keyword arguments to pass on to scipy.optimize.minimize.  Default
            values are method = "L-BFGS-B" and selected values of x0 and
            bounds that are based on computed results.  If the default
            variables are included in optimize_kwargs, they will be
            overridden.  The dictionary is passed to the optimization routine
            in kwargs format, aka **kwargs.
        ignore_censored: bool
            If True, censored values are ignored.
            
        Returns
        -------
        A dictionary containing various values of the analysis.  The "x"
        entry contains an ndarray in which the first entry is the estimated
        mean of the reference temperature and the second entry is the
        estimated standard deviation of the reference temperature.  Note that
        this is different from the standard ASTM E1921 analysis in that the
        mean and standard deviation are assumed to be a part of the
        generative process for the fracture toughness of an individual
        specimen drawn from the material.  The mean temperature in this case
        is indicative of the overall average toughness of a material, and the
        standard deviation is indicative of the inhomogeneity of the
        material.
        '''
        # run analysis of sample to insure self.results.rmat has been computed
        self.analyze(ignore_censored)
        
        # define helper functions for creating the final objective function as per Wallin paper
        K0 = lambda T0, T: 31.0 + 77.0 * np.exp(0.019 * (T - T0))
        fT = lambda T0, muT0, sigT0: norm.pdf(T0, loc = muT0, scale = sigT0)
        ST0 = lambda T0, T, KJc: np.exp(-((KJc - 20.0) / (K0(T0, T) - 20.0))**4)
        fT0 = lambda T0, T, KJc: 4.0 * (KJc - 20.0)**3 / (K0(T0, T) - 20.0)**4 * ST0(T0, T, KJc)
        S_int = lambda T0, muT0, sigT0, KJc, T: fT(T0, muT0, sigT0) * ST0(T0, T, KJc)
        # Si = lambda muT0, sigT0, KJc, T: quad(S_int, -200.0, 200.0, args = (muT0, sigT0, KJc, T), limit = 400)[0]
        f_int = lambda T0, muT0, sigT0, KJc, T: fT(T0 , muT0, sigT0) * fT0(T0, T, KJc)
        # fi = lambda muT0, sigT0, KJc, T: quad(f_int, -200.0, 200.0, args = (muT0, sigT0, KJc, T), limit = 400)[0]

        # define objective function to be optimized as per Wallin paper
        def negLL(params):
            integration_pts, half_width_sds = 1001, 12.0
            muT0, sigT0 = params
            half_width = half_width_sds * sigT0
            h = 2.0 * half_width / integration_pts
            T0 = np.linspace(muT0 - half_width, muT0 + half_width, integration_pts)
            T, KJc, delta = [self.results.rmat[:,col] for col in range(3)]
            fi = lambda kjc, t: simps(f_int(T0, muT0, sigT0, kjc, t), dx = h)
            fci = np.vectorize(fi)(KJc, T)
            fci[fci == 0] = 1e-300
            Si = lambda kjc, t: simps(S_int(T0, muT0, sigT0, kjc, t), dx = h)
            Sci = np.vectorize(Si)(KJc, T)
            Sci[Sci == 0] = 1e-300
            ln_L = np.sum(delta * np.log(fci) + (1 - delta) * np.log(Sci))
            return -ln_L
            
        # optimize objective function
        kwargs = dict(optimize_kwargs) if isinstance(optimize_kwargs, dict) else dict()
        if "method" in kwargs:
            method = kwargs["method"]
            del kwargs["method"]
        else:
            method = "L-BFGS-B"
        if "x0" in kwargs:
            x0 = kwargs["x0"]
            del kwargs["x0"]
        else:
            x0 = np.array([self.results.T0Q + 3 * self.results.sigma_T0Q, self.results.sigma_T0Q])
        if "bounds" in kwargs:
            bounds = kwargs["bounds"]
            del kwargs["bounds"]
        else:
            bounds = [(None, None), (1e-100, None)]
        mle_result = minimize(negLL, x0, method = method, bounds = bounds, **kwargs)
        
        # return final answer for T0, sigT0
        return {key:mle_result[key] for key in ["x", "fun", "success", "status", "message", "nit"]}
        
    def plot_data(self, alpha = 0.05, grid = True, legend = True, show = True, \
                  english_units = False):
        '''
        Plot the sample data based on standard ASTM E1921 analysis
        
        Arguments
        ---------
        alpha: float
            Prediction interval two-tailed tolerance bounds to be shown on
            the plot. Value should be in the interval (0, 1).  If None, then
            tolerance bounds will not be plotted.            
        grid: bool
            If True, a grid will be shown on the plot.            
        legend: bool
            If True, a legend will be shown on the plot.      
        show: bool
            If True, the plot will be shown on the current plotting device.        
        english_units: bool
            If True, the plot will be shown in English units.  Default is SI
            units.
        '''
        assert alpha is None or (alpha < 1.0 and alpha > 0.0), "Invalid alpha specification"
        if self.results is None:
            raise RuntimeError("Sample.analyze method must be called prior to Sample.plot_data")
        results = np.array(self.results.rmat)
        if english_units == True:
            tunits, Kunits = "$^{\circ}F$", "$ksi \cdot \sqrt{in}$"
            tv = celsius_to_fahrenheit(results[results[:,2] == 1,0])
            tc = celsius_to_fahrenheit(results[results[:,2] == 0,0])
            KJc_v = stress_intensity_in_SI_to_english(results[results[:,2] == 1,1])
            KJc_c = stress_intensity_in_SI_to_english(results[results[:,2] == 0,1])
            KJc_median_func = lambda T: stress_intensity_in_SI_to_english(self.results.KJc_median(fahrenheit_to_celsius(T)))
        else:
            tunits, Kunits = "$^{\circ}C$", "$MPa \cdot \sqrt{m}$"
            tv = results[results[:,2] == 1,0]
            tc = results[results[:,2] == 0,0]
            KJc_v = results[results[:,2] == 1,1]
            KJc_c = results[results[:,2] == 0,1]
            KJc_median_func = self.results.KJc_median
        fig = plt.figure(self.description + " Sample Data")
        axes = fig.add_subplot(111)
        axes.plot(tv, KJc_v, "ko", markerfacecolor = "none", label = "Valid Data", zorder = 5)
        axes.plot(tc, KJc_c, "kx", label = "Censored Data", zorder = 4)
        axes.set_xlim(xmin = min(np.concatenate([tv, tc])) - 10, xmax = max(np.concatenate([tv, tc])) + 10)
        plt.xlabel("Temperature (%s)" % (tunits,))
        plt.ylabel("$K_{J_c}$ (%s)" % (Kunits,))
        t = np.linspace(plt.xlim()[0], plt.xlim()[1], 101)
        KJc_median = KJc_median_func(t)
        axes.plot(t, KJc_median_func(t), "k-", zorder = 3)
        axes.set_ylim(ymin = 0.0, ymax = max(np.concatenate([KJc_median, KJc_v, KJc_c])) + 10)
        if alpha:
            tol = self.results.tol
            if english_units == True:
                lower = stress_intensity_in_SI_to_english(tol(fahrenheit_to_celsius(t), alpha / 2.0))
                upper = stress_intensity_in_SI_to_english(tol(fahrenheit_to_celsius(t), 1.0 - alpha / 2.0))
            else:
                lower = tol(t, alpha / 2.0)
                upper = tol(t, 1.0 - alpha / 2.0)
            axes.plot(t, lower, "k--", label = "{0:.0f}% Interval".format(100.0 - 100.0 * alpha), zorder = 2)
            axes.plot(t, upper, "k--", zorder = 1)
        plt.grid(grid)
        if legend == True:
            axes.legend(loc = "best", numpoints = 1)
        if show:
            plt.show()
        return fig, axes        