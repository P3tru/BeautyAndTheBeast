import numpy as np

### fitter implementations
from scipy.optimize import curve_fit


class NLLPDF:
    '''
    A callable object that evaluates the -log(likelihood) values given some binned data
    Does a few nonstandard things:
        - Forces a 'min_val' for each bin to avoid -log(0) to make low stats behave and avoid NaNs
        - Rolls off the probability at PDF boundaries to 'pull' outliers into the defined range
            (probability e^(-Ax) with x being distance from pdf boundary and A being 'pull')
    '''

    def __init__(self, counts, edges, min_val=None, pull=None, bounds=None):
        '''
        Set pull to None to estimate exponential falloff from data (if desired)
        Set min_val high enough to avoid poisson fluctuations near edges of PDF (if desired)
        '''
        self.nll = counts[:]
        if bounds is not None:
            self.bounds(bounds[0], bounds[1])
        if min_val is None:
            min_val = np.max(self.nll) / 100
        pull_bounds = np.argwhere(self.nll > min_val)
        left, right = pull_bounds[0][0], pull_bounds[-1][0]
        widths = edges[1:] - edges[:-1]
        norm = np.sum(self.nll * widths)
        self.nll = -np.log(self.nll[left:right] / norm)
        self.min_val = np.min(self.nll)
        self.centers = (edges[:-1] + edges[1:])[left:right] / 2
        if pull is None:
            npts = min(len(self.nll) // 2, 10)
            left_pull = -np.mean(
                (self.nll[:npts - 1] - self.nll[1:npts]) / (self.centers[:npts - 1] - self.centers[1:npts]))
            right_pull = np.mean(
                (self.nll[-npts:-1] - self.nll[-npts + 1:]) / (self.centers[-npts:-1] - self.centers[-npts + 1:]))
            self.pull = (left_pull, right_pull)
        else:
            self.pull = (pull, pull)
        min_val = -np.log(min_val / norm)
        self.nll[self.nll > min_val] = min_val

    def __call__(self, x):
        vec = np.interp(x, self.centers, self.nll)
        greater = x > self.centers[-1]
        vec[greater] = self.nll[-1] + (x[greater] - self.centers[-1]) * self.pull[1]
        lesser = x < self.centers[0]
        vec[lesser] = self.nll[0] + (self.centers[0] - x[lesser]) * self.pull[0]
        return vec

    def mean(self):
        weights = np.exp(-self.nll)
        return np.sum(self.centers * weights) / np.sum(weights)


def polyform(x, *args):
    ''' args are cher_ang,cher_off,[poly_left],[poly_right]
        poly_ lacks the offset which is constrained to be cher_off for both
        poly_ is expanded around cher_ang'''
    total = len(args)
    poly = int((total - 1) / 2)

    cher_ang, cher_off = args[:2]
    poly_left = np.concatenate([args[1:1 + poly], [cher_off]])
    poly_right = np.concatenate([args[1 + poly:], [cher_off]])
    results = np.empty_like(x)
    mask = x <= cher_ang
    results[mask] = np.polyval(poly_left, x[mask] - cher_ang)
    mask = np.logical_not(mask)
    results[mask] = np.polyval(poly_right, x[mask] - cher_ang)

    return results


class CosThetaPDF:
    '''
        A callable object that evaluates the -log(likelihood) values by fitting
        Nth degree polynomials above and below the cherenkov peak of binned data
    '''

    def __init__(self, counts, edges, order=5):
        nll = counts[:]
        widths = edges[1:] - edges[:-1]
        norm = np.sum(nll * widths)
        nll = -np.log(nll / norm)
        min_val = np.min(nll)
        centers = (edges[:-1] + edges[1:]) / 2

        cher_peak = centers[np.argmin(nll)]
        poly_right = np.zeros(order)
        poly_left = np.zeros(order)

        p, p_cov = curve_fit(polyform, centers, nll, p0=np.concatenate([[cher_peak], poly_left, poly_right]), maxfev=20000)

        self.params = p

    def __call__(self, x):
        return polyform(x, *self.params)
