#!/usr/bin/env pytom

"""
Created on Jun 27, 2011
Updated on Jun 21, 2021

@author: yuxiangchen, Marten Chaillet
"""

# plotting
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# main functionality
import sys
import numpy as np
import scipy.interpolate as interpolate
from scipy.optimize import curve_fit
from scipy.special import erf
from pytom.localization.structures import readParticleFile


if __name__ == '__main__':

#     # parse command line arguments
#     import sys
#     from pytom.tools.script_helper import ScriptHelper, ScriptOption
#     from pytom.tools.parse_script_options import parse_script_options
#     helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
#                           description='Do the Gaussian fitting on the found particle list.',
#                           authors='Yuxiang Chen',
#                           options=[ScriptOption(['-f','--file'], 'Particle list after extracting candidates.', True, False),
#                                    ScriptOption(['-n','--numberBins'], 'Number of bins of histogram. Default is 10.', True, True),
#                                    ScriptOption(['-p','--gaussianPeak'], 'The correspondent index of the gaussian peak.', True, False),
#                                    ScriptOption(['-c','--min'], 'Only use particles with scores higher than cutoff.', True, True),
#                                    ScriptOption(['-m','--max'], 'Only use particles with scores lower than cutoff.', True, True),
#                                    ScriptOption(['-i','--imageFile'], 'Save plot to a image file.', True, True),
#                                    ScriptOption(['-h', '--help'], 'Help.', False, True)])
#     if len(sys.argv) == 1:
#         print(helper)
#         sys.exit()
#     try:
#         pl_filename, num_steps, peak_index, min_cut, max_cut, imageFile, help = parse_script_options(sys.argv[1:], helper)
#     except Exception as e:
#         print(e)
#         sys.exit()
#     if help is True:
#         print(helper)
#         sys.exit()
#
#     # process the arguments
#     num_steps = int(num_steps)
#     if not num_steps:
#         num_steps = 10
#
#     peak_index = int(peak_index)
#
#     scores = []
#
    # parse command line arguments with ScriptHelper2

    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Fit a Gaussian to particle list peak.',
        authors='Marten Chaillet (based on Yuxiang Chen)',
        options=[
            ScriptOption2(['-f', '--file'], 'Particle list xml file', 'file', 'required'),
            ScriptOption2(['-n', '--numberBins'], 'Numbers of bins to plot the score distribution', 'int',
                          'optional', 10),
            ScriptOption2(['-p', '--gaussianPeak'], 'Index of the peak (p < nbins) to be taken as a starting point '
                                                    'for searching the mean of the particle population', 'int',
                          'optional'),
            ScriptOption2(['--forcePeak'], 'Force the given peak index to be the mean', 'no arguments', 'optional'),
            ScriptOption2(['-c', '--numberParticles'], 'Return the number of particles up to this correlation score',
                          'float', 'optional'),
            ScriptOption2(['-i', '--imageFile'], 'Save plot to an image file; if not provided plot will be shown in '
                                                 'window', 'string', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    pl_filename, num_bins, peak_index, force_peak, ccc_value, imageFile = options

    # parse peak index
    if (peak_index is None) or (peak_index > num_bins):
        print('fall back to default peak index estimate, nbins / 2')
        peak_index = num_bins // 2


    # read out the scores
    foundParticles = readParticleFile(pl_filename)
    scores = []
    for f in foundParticles:
        scores.append(float(f.score.getValue()))

    #
    # # only use scores larger than ccc_cutoff
    # if min_cut: min_cut = float(min_cut)
    # else: min_cut = 0
    # if max_cut : max_cut = float(max_cut)
    # else: max_cut = 1
    # scores = [s for s in scores if min_cut < s < max_cut]
    #
    # # construct x and y array according to the given peak index
    # scores.sort()
    # min = scores[0]
    # max = scores[-1]
    #
    # step = (max-min)/num_steps
    # x = []
    # for i in range(num_steps):
    #     x.append(min+i*step)
    # x.append(max)
    #
    # y = []
    # for i in range(num_steps):
    #     lower = x[i]; upper = x[i+1]
    #     n = len([v for v in scores if lower<=v<=upper])
    #     y.append(n)
    #
    #
    # # plot
    # import matplotlib
    # matplotlib.use('Qt5Agg')
    # from matplotlib import pyplot
    # import numpy
    # matplotlib.rc('lines', linewidth=2)
    # matplotlib.rc('font', size=24)
    # fig = pyplot.figure()
    # plt = fig.add_subplot(111)
    # plt.plot(x[1:],y,'ro-')
    # plt.set_xlabel('Score')
    # plt.set_ylabel('Frequency')
    #
    # # do the fitting
    # from pytom.tools.maths import gaussian_fit
    # sigma, mu, a = gaussian_fit(x[peak_index:], y[peak_index-1:])
    #
    # if sigma.__class__ in [complex,numpy.complex128]:
    #     sigma = sigma.real
    # if mu.__class__ in [complex,numpy.complex128]:
    #     mu = mu.real
    # print('sigma: %f, mu: %f, a: %f' % (sigma, mu, a))
    #
    # # plot the Gaussian fitting
    # from math import exp
    # gaussian_fnc = lambda x: a*exp(-(x-mu)**2/(2*sigma**2))
    # plt.plot(x[1:],list(map(gaussian_fnc, x[1:])),'g--')
    #
    # # print the estimation of number of true positives
    # x.reverse()
    # while True:
    #     if x[-1] > 0:
    #         x.append(x[-1]-step)
    #     else:
    #         x = x[:-1]
    #         break
    #
    # estimate = 0.
    # for i in x:
    #     if i > mu-sigma:
    #         estimate += gaussian_fnc(i)
    #     else:
    #         break
    # print('One sigma position: %f, number of estimation: %f' % (mu-sigma, estimate))
    #
    # estimate = 0.
    # for i in x:
    #     if i > mu-2*sigma:
    #         estimate += gaussian_fnc(i)
    #     else:
    #         break
    # print( 'Two sigma position: %f, number of estimation: %f' % (mu-2*sigma, estimate))
    #
    # #if ccc_value:
    # #    ccc_value = float(ccc_value)
    # #    estimate = 0.
    # #    for i in x:
    # #        if i > ccc_value:
    # #            estimate += gaussian_fnc(i)
    # #        else:
    # #            break
    # #    print('CCC value position: %f, number of estimation: %f' % (ccc_value, estimate))
    #
    # if imageFile is None:
    #     pyplot.show()


    # generate the histogram
    matplotlib.rc('font', size=18)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(121)
    # ax.hist(x[1:],y,'ro-')
    y, x_hist, _ = ax.hist(scores, bins=num_bins, histtype='step') #, 'ro-')
    ax.set_xlabel('Score')
    ax.set_xlim(x_hist[0], x_hist[-1])
    ax.set_ylabel('Frequency')

    # define gaussian function with parameters to fit
    def gauss(x, mu, sigma, A):
        return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

    # integral of gaussian with certain sigma and A
    def gauss_integral(sigma, A):
        # mu does not influence the integral
        return A * np.abs(sigma) * np.sqrt(2 * np.pi)

    # exponential function
    def exponential(x, a, b, c):
        return a * np.exp(-b * x) + c

    # a cumulative distribution function
    def cumulative_dist(x, Lambda):
        return 1 - np.exp(- Lambda * x)

    # define bimodal function of two gaussians to fit both populations
    def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
        return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

    try:
        # adjust x to center of each bin so len(x)==len(y)
        x = (x_hist[1:] + x_hist[:-1]) / 2

        # expected values
        # left gaussian expectation: mu < x[0] and A > y[0]
        # right gaussian expectation: mu ~ x[half] and A ~ y[half]
        # initial guess for mu_1 should be 0 for better convergence
        expected = (0, .1, y[0], x[peak_index], .1, y[peak_index])
        # force peak of particle population to be at peak index
        if force_peak:
            bounds = ([0, 0, 0, x[peak_index] - 0.01, 0, 0],
                      [.1, 0.2, np.inf, x[peak_index] + 0.01, 0.2, 200])
        else:
            bounds = (0, [.1, 0.2, np.inf, 1.0, 0.2, 200])  # todo 200 is not a good upper bound for second gaussian
        # parameter names for output
        params_names = ['mu_1', 'sigma_1', 'A_1', 'mu_2', 'sigma_2', 'A_2']
        # params_names = ['a', 'b', 'c', 'mu_2', 'sigma_2', 'A_2']
        params, cov = curve_fit(bimodal, x, y, p0=expected, bounds=bounds)  # max iterations argument: maxfev=2000)
        # give sigma of fit for each parameter
        sigma = np.sqrt(np.diag(cov))

        # print information about fit of the model
        print('\nfit of the bimodal model:')
        print('\testimated\t\tsigma')
        for n, p, s in zip(params_names, params, sigma):
            print(f'{n}\t{p:.3f}\t\t{s:.3f}')
        print('\n')

        # plot bimodal model and the gaussian particle population
        ax.plot(x, bimodal(x, *params), color='blue', lw=2, label='bimodal model')
        population, noise = (params[:3], params[3:6]) if params[0] > params[3] else (params[3:6], params[:3])
        # population = params[2:5]
        ax.plot(x, gauss(x, *population), color='red', lw=2, label='particle population')
        ax.set_ylim(0, 3 * population[2])
        ax.legend()

        # todo add proposed cutoff based on function overlap?

        # print the estimation of number of true positives
        # x = np.flip(x)
        # gaussian_fnc = lambda x: gauss(x, *population)
        # mu, sigma = population[0], abs(population[1])
        # estimate = 0.
        # for i in x:
        #     if i > mu - sigma:
        #         estimate += gaussian_fnc(i)
        #     else:
        #         break
        # print('One sigma position: %f, particle number estimation: %f' % (mu - sigma, estimate))
        #
        # estimate = 0.
        # for i in x:
        #     if i > mu - 2 * sigma:
        #         estimate += gaussian_fnc(i)
        #     else:
        #         break
        # print( 'Two sigma position: %f, particle number estimation: %f' % (mu - 2 * sigma, estimate))

        # if ccc_value:
        #     ccc_value = float(ccc_value)
        #     estimate = 0.
        #     for i in x:
        #         if i > ccc_value:
        #             estimate += gaussian_fnc(i)
        #         else:
        #             break
        #     print('CCC value position: %f, number of estimation: %f' % (ccc_value, estimate))

        # Generate a ROC curve
        roc_steps = 200
        x_roc = np.flip(np.linspace(x[0], x[-1], roc_steps))
        # find ratio of hist step vs roc step
        hist_step = (x_hist[-1] - x_hist[0]) / num_bins
        roc_step = (x[-1] - x[0]) / roc_steps
        delta = hist_step / roc_step  # can be used to divide true/false positives by per roc step
        # variable for total number of tp and fp
        n_true_positives = .0
        n_false_positives = .0
        # list for storing probability of true positives and false positives for each cutoff
        recall = []  # recall = TP / (TP + FN); TP + FN is the full area under the Gaussian curve
        fdr = []  # false discovery rate = FP / (TP + FP); == 1 - TP / (TP + FP)

        # find integral of gaussian particle population; NEED TO DIVIDE BY HISTOGRAM BIN STEP
        population_integral = gauss_integral(population[1], population[2]) / hist_step
        print('estimation total number of true positives: ', population_integral)

        # should use CDF (cumulative distribution function) of Gaussian, gives probability from -infinity to x
        CDF = lambda x: 0.5 * (1 + erf((x - population[0])/(np.sqrt(2) * population[1])))
        gauss_noise = lambda x: gauss(x, *noise)
        # alternative (less precise) way of determining true positives
        # gauss_pop = lambda x: gauss(x, *population)

        for x_i in x_roc:
            # calculate probability of true positives x_i
            # n_true_positives += gauss_pop(x_i) / delta
            n_true_positives = (1 - CDF(x_i)) * population_integral

            # determine false positives up to this point, could also use CDF
            n_false_positives += gauss_noise(x_i) / delta

            # add probability
            recall.append(n_true_positives / population_integral)
            fdr.append(n_false_positives / (n_true_positives + n_false_positives))

        # ============= attempt to fit ROC curve by cumulative distribution function of exponential distribution =======
        # # ROC curve characterized by f(x) = 1 - e^(-lx)  (l=lambda)
        # params2_names = ['lambda']
        # params2, cov2 = curve_fit(cumulative_dist, [f for f,r in zip(fdr, recall) if f >0.05],
        #                           [r for f, r in zip(fdr, recall) if f > 0.05], p0=[1])
        # # give sigma of fit for each parameter
        # sigma2 = np.sqrt(np.diag(cov2))
        #
        # # print information about fit of the model
        # # print('\nfit of the cumulative distribution:')
        # # print('\testimated\t\tsigma')
        # for n, p, s in zip(params2_names, params2, sigma2):
        #     print(f'{n}\t{p:.3f}\t\t{s:.3f}')
        # print('\n')
        # x_cum = np.linspace(0,1,200)
        # y_cum = cumulative_dist(x_cum, *params2)

        # ============================= use splines to fit ROC curve ===================================================
        recall.append(1.)
        fdr.append(1.)
        recall.insert(0, .0)
        fdr.insert(0, .0)
        # use CubicSpline this fits the (sort of) cumulative dist best)
        spline = interpolate.CubicSpline(fdr, recall)
        # scipy CubicSpline can be easily integrated to find area under curve
        AUC = spline.integrate(0, 1)
        print('AUC: ', AUC)

        # points for plotting
        xs = np.linspace(0, 1, 200)
        # plot the fdr curve
        ax2 = fig.add_subplot(122)
        ax2.scatter(fdr, recall, facecolors='none', edgecolors='r', s=25)
        ax2.plot([0, 1], [0, 1], ls="--", c=".3", lw=1) # transform=ax.transAxes,
        ax2.plot(xs, spline(xs), label='spline')
        # ax2.plot(x_cum, y_cum, label='spline')
        ax2.plot([], [], ' ', label=f'AUC: {AUC:.2f}')
        ax2.legend()
        ax2.set_xlabel('False discovery rate')
        ax2.set_ylabel('Recall')
        ax2.set_xlim(0, 1)
        ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        ax2.set_ylim(0, 1)
        ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    except RuntimeError as e:
        # runtime error is because the model could not be fit, in that case print error and continue with execution
        print(e)

    if imageFile is None:
        plt.tight_layout()
        plt.show()

    # else:
    #     if not ('png' in imageFile or 'PNG' in imageFile):
    #         imageFile = imageFile + '.png'
    #
    #     plt.tight_layout()
    #     plt.savefig(imageFile)
