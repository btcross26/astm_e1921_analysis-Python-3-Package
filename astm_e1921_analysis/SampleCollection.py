from .helper import stress_intensity_in_SI_to_english, fahrenheit_to_celsius
from .Sample import Sample
import matplotlib.pyplot as plt
import numpy as np
    
class SampleCollection(object):
    '''
    SampleCollection object
    Object that contains a collection of test sample and means to plot all samples
    on the same axes after normalizing the reference temperature to zero.
    '''
    def __init__(self, description):
        '''
        Constructor Arguments:
            description = string describing the sample collection
        '''
        self.description = description
        self.samples = dict()
        
    def __str__(self):
        return "SampleCollection(description=\"{0:s}\")".format(self.description)
        
    def __repr__(self):
        return "<SampleCollection at {0:s}>".format(hex(id(self)))
        
    def __iter__(self):
        for key in self.samples:
            yield self.samples[key]
            
    def __getitem__(self, sample_description):
        return self.samples[sample_description]
        
    def __len__(self):
        return len(self.samples)
        
    def add_samples(self, *samples):
        '''
        Method to add Sample objects to the collection
        Arguments:
            *samples = a comma separated list of 1 or more sample objects
        '''
        for sample in samples:
            assert isinstance(sample, Sample), "<arg>:*samples must be of type Sample"
        for sample in samples:    
            self.samples[sample.description] = sample
        
    def remove_samples(self, *sample_descriptions):
        '''
        Method to remove Sample objects from the collection based on sample description strings
        Arguments:
            *sample_description = a comma separated list of 1 or more sample description strings
        '''
        for sample_description in sample_descriptions:
            if sample_description in self.samples:
                del self.samples[sample_description]
            else:
                raise ValueError("Sample(description=\"{0:s}\") is not part of the sample collection".format(sample_description))
        
    def plot_data(self, alpha = 0.05, grid = True, legend = True, show = True, english_units = False, ignore_censored = False):
        '''
        Method to plot all samples in the collection against normalized temperature.  Note that
        all samples will be plotted with the same marker and no legend will be included on the plot,
        although each line object is labeled with the sample description such that a legend can
        be added accordingly.
        Returns:
            A tuple of figure, axes for the given plot
        '''
        assert alpha is None or (alpha < 1.0 and alpha > 0.0), "Invalid alpha specification"
        assert len(self) > 0, "SampleCollection is empty"
        
        # create figure object
        fig = plt.figure(self.description + " SampleCollection")
        
        # add axes to figure and label axes
        axes = fig.add_subplot(111)
        if english_units == True:
            tunits, Kunits = "$^{\circ}F$", "$ksi \cdot \sqrt{in}$"
        else:
            tunits, Kunits = "$^{\circ}C$", "$MPa \cdot \sqrt{m}$"
        plt.xlabel("$T - T_0$ (%s)" % (tunits,))
        plt.ylabel("$K_{J_c}$ (%s)" % (Kunits,))
        
        # plot valid data points for each sample
        zorder = len(self.samples) + 3
        temps, kjc = np.array([]), np.array([])
        for sample in self:
            sample.analyze(ignore_censored = ignore_censored)
            results = np.array(sample.results.rmat)
            tv = results[results[:,2] == 1,0] - sample.results.T0Q   # represents a change in temp rather than actual temp
            KJc_v = results[results[:,2] == 1,1]
            if english_units == True:
                tv = tv * 1.80   # convert a CHANGE in temp to Fahrenheit rather than actual temp
                KJc_v = stress_intensity_in_SI_to_english(KJc_v)
            temps = np.concatenate([temps, tv])
            kjc = np.concatenate([kjc, KJc_v])
            axes.set_xlim(xmin = np.min(temps) - 10.0, xmax = np.max(temps) + 10.0)
            axes.plot(tv, KJc_v, "ko", zorder = zorder, markerfacecolor = "none", label = "_" + sample.description)
            zorder -= 1
        
        # add master curve to plot
        t = np.linspace(axes.get_xlim()[0], axes.get_xlim()[1], 101)
        if english_units == True:
            T0Q = celsius_to_fahrenheit(sample.results.T0Q)
            KJc_median_func = lambda T: stress_intensity_in_SI_to_english(sample.results.KJc_median(fahrenheit_to_celsius(T + T0Q)))
        else:
            T0Q = sample.results.T0Q
            KJc_median_func = lambda T: sample.results.KJc_median(T + T0Q)
        KJc_median = KJc_median_func(t)
        kjc = np.concatenate([kjc, KJc_median])
        axes.plot(t, KJc_median, "k-", zorder = zorder)
        zorder -= 1
        axes.set_ylim(ymin = 0.0, ymax = np.max(kjc) + 10.0)
        if alpha:
            tol = sample.results.tol
            if english_units == True:
                lower = stress_intensity_in_SI_to_english(tol(fahrenheit_to_celsius(t + T0Q), alpha / 2.0))
                upper = stress_intensity_in_SI_to_english(tol(fahrenheit_to_celsius(t + T0Q), 1.0 - alpha / 2.0))
            else:
                lower = tol(t + T0Q, alpha / 2.0)
                upper = tol(t + T0Q, 1.0 - alpha / 2.0)
            axes.plot(t, lower, "k--", label = "{0:.0f}% Interval".format(100.0 - 100.0 * alpha))
            axes.plot(t, upper, "k--")
        plt.grid(grid)
        if legend == True and alpha:
            axes.legend(loc = "best", numpoints = 1)
        if show == True:
            plt.show()
        return fig, axes
 