# This code is heavily based on a blog post by Andre Dietrich, currently
# found at http://www.aizac.info/simple-check-of-a-sample-against-80-distributions/

import scipy.stats
import warnings
from math import log
from os import listdir
import numpy as np
import sys

# just for surpressing warnings
warnings.simplefilter('ignore')

  
# list of all available distributions
cdfs = [
    "norm",            #Normal (Gaussian)
    "alpha",           #Alpha
    "anglit",          #Anglit
    "arcsine",         #Arcsine
    "beta",            #Beta
    "betaprime",       #Beta Prime
    "bradford",        #Bradford
    "burr",            #Burr
    "cauchy",          #Cauchy
    "chi",             #Chi
    "chi2",            #Chi-squared
    "cosine",          #Cosine
    "dgamma",          #Double Gamma
    "dweibull",        #Double Weibull
    "erlang",          #Erlang
    "expon",           #Exponential
    "exponweib",       #Exponentiated Weibull
    "exponpow",        #Exponential Power
    "fatiguelife",     #Fatigue Life (Birnbaum-Sanders)
    "foldcauchy",      #Folded Cauchy
    "f",               #F (Snecdor F)
    "fisk",            #Fisk
    "foldnorm",        #Folded Normal
    "frechet_r",       #Frechet Right Sided, Extreme Value Type II
    "frechet_l",       #Frechet Left Sided, Weibull_max
    "gamma",           #Gamma
    "gausshyper",      #Gauss Hypergeometric
    "genexpon",        #Generalized Exponential
    "genextreme",      #Generalized Extreme Value
    "gengamma",        #Generalized gamma
    "genlogistic",     #Generalized Logistic
    "genpareto",       #Generalized Pareto
    "genhalflogistic", #Generalized Half Logistic
    "gilbrat",         #Gilbrat
    "gompertz",        #Gompertz (Truncated Gumbel)
    "gumbel_l",        #Left Sided Gumbel, etc.
    "gumbel_r",        #Right Sided Gumbel
    "halfcauchy",      #Half Cauchy
    "halflogistic",    #Half Logistic
    "halfnorm",        #Half Normal
    "hypsecant",       #Hyperbolic Secant
    "invgamma",        #Inverse Gamma
#    "invnorm",         #Inverse Normal
    "invweibull",      #Inverse Weibull
    "johnsonsb",       #Johnson SB
    "johnsonsu",       #Johnson SU
    "laplace",         #Laplace
    "logistic",        #Logistic
    "loggamma",        #Log-Gamma
    "loglaplace",      #Log-Laplace (Log Double Exponential)
    "lognorm",         #Log-Normal
    "lomax",           #Lomax (Pareto of the second kind)
    "maxwell",         #Maxwell
    "mielke",          #Mielke's Beta-Kappa
    # "nakagami",        #Nakagami
    "ncx2",            #Non-central chi-squared
#    "ncf",             #Non-central F
    "nct",             #Non-central Student's T
    "pareto",          #Pareto
    "powerlaw",        #Power-function
    "powerlognorm",    #Power log normal
    "powernorm",       #Power normal
    "rdist",           #R distribution
    "reciprocal",      #Reciprocal
    "rayleigh",        #Rayleigh
    "rice",            #Rice
    "recipinvgauss",   #Reciprocal Inverse Gaussian
    "semicircular",    #Semicircular
    "t",               #Student's T
    "triang",          #Triangular
    "truncexpon",      #Truncated Exponential
    "truncnorm",       #Truncated Normal
    "tukeylambda",     #Tukey-Lambda
    "uniform",         #Uniform
    "vonmises",        #Von-Mises (Circular)
    "wald",            #Wald
    "weibull_min",     #Minimum Weibull (see Frechet)
    "weibull_max",     #Maximum Weibull (see Frechet)
    "wrapcauchy",      #Wrapped Cauchy
    "ksone",           #Kolmogorov-Smirnov one-sided (no stats)
    "kstwobign"]       #Kolmogorov-Smirnov two-sided test for Large N
    

t = sys.argv[1:]
dict = dict((cdf,[0,0,0]) for cdf in cdfs)
for filename in t:
    print "processing %s" % filename
    data = np.loadtxt(filename)
    new_cdfs=[]
    for cdf in cdfs:
        #fit our data set against every probability distribution
        parameters = eval("scipy.stats."+cdf+".fit(data)");
        #Applying the Kolmogorov-Smirnof one sided test
        D, p = scipy.stats.kstest(data, cdf, args=parameters);
        # Discard all distributions that are very poor fits
        if np.isnan(p):
            pass
        else:
            print cdf
            negative_log_likelihood = eval("scipy.stats."+cdf+".nnlf("+str(parameters)+",data)")
            if negative_log_likelihood<1e10 and p>0:
                new_cdfs.append(cdf)
                dict[cdf][0]-=log(p)
                dict[cdf][1]+=negative_log_likelihood
                dict[cdf][2]=len(parameters)
    cdfs=new_cdfs[:]

print "\n\nBest fits for selected data set. Contains name,\nsum of negative log p-values for KS test (lower values indicate better fit),\nsum of negative log likelihoods (lower is better),\nand number of parameters (degrees of freedom) \n"

print "name".ljust(12) + "neg log p KS-test".ljust(21) + "neg log likelihood".ljust(21)+"# parameters"

for cdf in new_cdfs:
    print cdf.ljust(14) + str(dict[cdf][0]).ljust(20) + str(dict[cdf][1]).ljust(24), str(dict[cdf][2])

