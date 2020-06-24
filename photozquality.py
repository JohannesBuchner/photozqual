"""
Copyright 2013-2019 Johannes Buchner

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Loads pdzs, and computes the best probability functions to incorporate systematic errors.

Usage: 
python photozquality.py ids...

Optimizes the parameters under the likelihood of the zspec values, 
giving the most likely description.

"""

import sys, os
import numpy
import scipy
import scipy.signal
import scipy.stats
import scipy.optimize
import matplotlib.pyplot as plt

def get_zphot(id):
	zphot1 = numpy.loadtxt('pdz/%s' % id)
	zphot1[:,1] /= zphot1[:,1].sum()
	return zphot1
def get_zspec(id):
	return float(numpy.loadtxt('specz/%s' % id))
def has_zspec(id):
	return os.path.exists('specz/%s' % id)


def weighted_median(y, p):
	items = sorted(zip(y,p), key=lambda x: x[0])
	s = 0.
	for i in items:
		s += i[1]
		if s >= 0.5: 
			return i[0]
def normalize(x, y):
	binwidth = x[1:] - x[:-1]
	total = (y[:-1] * binwidth).sum()
	assert total >= 0, (total, x, y)
	return y / total
data = {}

def create_pdz(zphot, w, woutliers):
	zphot_best = zphot[:,0][zphot[:,1].argmax()]
	zphot2a = zphot[:,1]
	zphot2 = scipy.signal.convolve(zphot2a, w, 'same')
	assert zphot2.sum() > 0, w
	zphot2 /= zphot2.sum()
	#zphot2 = zphot2**d # d disabled
	zphot3 = normalize(zphot[:,0], zphot2)
	zphot3 += woutliers
	#zphot3[0,0] = 0.
	zphot3[-1] = 0.
	zphot4 = normalize(zphot[:,0], zphot3)
	return zphot_best, zphot2, zphot3, zphot4


def calc_stats(k, wlogoutliers, plot_individual=False, write_individual=False, plot_corr=False):
	values = []
	woutliers = 10**(wlogoutliers-9)
	if not (k >= 0 and woutliers >= 0):
		print('ERROR: negative values passed')
		return {}
	distvalues = []
	outlier_list = []
	# k == 50 means sigma of 0.5 --> k * 100 gives width
	w = scipy.signal.gaussian(100, k / 100 * 100)
	if not (w.sum() > 0):
		print('ERROR: kernel width too small')
		return {}
	w /= w.max()
	
	#print(k, wlogoutliers)
	
	data_pdzs = dict([(i, create_pdz(zphot, w, woutliers)) for (i, (zphot, zspec)) in data.items()])
	
	for i, (zphot, zspec) in data.items():
		zphot_best, zphot2, zphot3, zphot4 = data_pdzs[i]
		
		#print 'norm:', numpy.trapz(x=zphot[:,0], y=zphot2)
		
		distvalues += list(zip([zspec] * len(zphot2), zphot[:,0], zphot2, numpy.abs(zspec - zphot[:,0]) / (1 + zspec), ((zspec - zphot[:,0]) / (1 + zspec) )**2))
		
		#print zphot
		zmed, zlo, zhi = numpy.interp(x=[0.5, 0.1, 0.9], xp=zphot2.cumsum(), fp=zphot[:,0])

		# gaussian kernel of 0.3 width

		#zphot_val = numpy.interp(zspec, zphot[:,0], zphot[:,1])
		#if zphot_val == 0:
		#	outlier_list.append(i)
		# max(1e-6, zphot_val)
		zphot_val = numpy.interp(zspec, zphot[:,0], zphot4)
		avgdist = numpy.average(numpy.abs(zphot[:,0] - zspec), weights = zphot2)
		
		#print zphot[:,1].sum(), zspec, zphot_val, max(1e-9, zphot_val)

		values.append( (zspec, zphot_best, zphot_val, zmed, zlo, zhi, avgdist) )
	
	if plot_individual or write_individual:
		for i in allids:
			zphot = get_zphot(i)
			zphot_best, zphot2, zphot3, zphot4 = create_pdz(zphot, w, woutliers)
			if plot_individual:
				plt.figure(figsize=(4,3))
				#plt.gca().set_yscale('log')
				plt.plot(zphot[:,0], zphot[:,1], label='photo-z', color='black', lw=2)
				#plt.plot(zphot[:,0], zphot2 / zphot2.max(), label='convolved')
				plt.plot(zphot[:,0], zphot3, label='$p_{D+Sys}(z)$', color='gray', ls='--')
				plt.fill_between(zphot[:,0], zphot3, color='gray', alpha=0.2)
				#plt.plot(zphot[:,0], zphot2/ zphot2.max(), label='zphot final')
				zeroes = zphot[zphot4 > max(zphot4) * 1e-2,0]
				plt.xlim(max(min(zeroes)-0.1, 0), min(max(zeroes)+0.1,7))
				plt.ylim(0, None)
				plt.yticks([])
				plt.vlines(zspec, plt.ylim()[0], plt.ylim()[1], linestyles='--')
				plt.ylabel('probability')
				plt.xlabel('redshift')
				plt.legend(loc='best', prop=dict(size=10), frameon=False)
				plt.savefig('pdz/' + i + ".png")
				plt.savefig('pdz/' + i + ".pdf", bbox_inches='tight')
				plt.close()
			if write_individual:
				numpy.savetxt('smoothened/' + i, numpy.vstack([zphot[:,0], zphot4]).transpose())
	
	v = numpy.array(values, dtype=[ ('zspec', 'f'), 
		('zphot_best', 'f'), # original best-fit photo-z
		('zphot_val', 'f'),  # pdz evaluated at zspec
		('zphot_med', 'f'), ('zphot_lo', 'f'), ('zphot_hi', 'f'), # limits
		('zdist', 'f') # weighted distance to zspec
		] )
	
	# compute usual statistics (for unconvolved)
	dist = numpy.abs(v['zspec'] - v['zphot_best']) / (1 + v['zspec'])

	outliers = dist > 0.15
	outlier_fraction = outliers.sum() * 1. / len(outliers)

	sigma_NMAD = 1.48 * numpy.median( dist )
	
	distvalues = numpy.array(distvalues, dtype=[ ('zspec', 'f'), ('zphot', 'f'), ('prob', 'f'), ('zabsdist', 'f'), ('zsqdist', 'f')])
	snad = (distvalues['zabsdist'] * distvalues['prob']).sum() / len(v)
	sigma_ad = (distvalues['zsqdist'] * distvalues['prob']).sum()**0.5 / len(v)
	if plot_corr:
		nmad = weighted_median(distvalues['zabsdist'], distvalues['prob'] / len(v))
	smooth_outlier_fraction = ((distvalues['zabsdist'] > 0.15) * distvalues['prob']).sum() / len(v)
	
	# compute likelihood
	loglike = numpy.sum(numpy.log(v['zphot_val']))
	information_loss = numpy.sum(v['zphot_val'] * numpy.log2(v['zphot_val']))
	# check the values where the likelihood is zero
	
	r = dict(n=len(outliers), outlier_fraction=outlier_fraction, 
		loglike=loglike, sigma_NMAD=sigma_NMAD, snad=snad, #nmad=nmad,
		outliers=outlier_list, sigma_ad=sigma_ad, smooth_outlier_fraction=smooth_outlier_fraction,
		smoothened_chi2 = -2 * loglike,
		information_loss = information_loss,
		)
	
	#print(r)

	if plot_corr:
		print()
		print("plotting...")
		plt.figure(figsize=(5,3))
		plt.plot(v['zspec'], v['zphot_best'], 'x ')
		plt.plot([0,plt.xlim()[1]], [0,plt.xlim()[1]], '--', color='grey')
		plt.xlabel('Spectroscopic redshift')
		plt.ylabel('Photometric redshift')
		plt.ylim(0, 4)
		plt.xlim(0, 4)
		plt.xticks(list(range(5)))
		plt.yticks(list(range(5)))
		plt.text(0.2, plt.ylim()[1]-0.2, "$\\eta = %.1f\\%% $\n$\\sigma_{NMAD} = %.4f$\n$n = %d$" % (outlier_fraction*100, sigma_NMAD, len(v)), va='top')
		plt.savefig('photozquality.pdf', bbox_inches='tight')
		plt.savefig('photozquality.png', bbox_inches='tight')
		plt.close()
		
		plt.figure()
		plt.errorbar(v['zspec'], v['zphot_med'], 
			yerr=[v['zphot_med'] - v['zphot_lo'], v['zphot_hi'] - v['zphot_med'] ], 
			marker='x', ls='')
		plt.plot([0,plt.xlim()[1]], [0,plt.xlim()[1]], '--', color='grey')
		plt.xlabel('Spectroscopic redshift')
		plt.ylabel('Photometric redshift')
		plt.text(0.2, plt.ylim()[1]-0.2, "$\\chi^2 = %.4f$\n$n = %d$" % (-2 * loglike, len(v)), va='top')
		plt.savefig('photozquality_err.pdf', bbox_inches='tight')
		plt.savefig('photozquality_err.png', bbox_inches='tight')
		plt.close()
		
		plt.figure(figsize=(3,2.5))
		plt.hexbin(x=distvalues['zspec'], y=distvalues['zphot'], 
			C=distvalues['prob'], gridsize=30, 
			reduce_C_function=lambda x: numpy.log(1+numpy.sum(x)),
			cmap = plt.cm.gray_r, vmax=numpy.log(1+len(v)/30.*4), extent=[0,4,0,4])
		plt.plot([0,plt.xlim()[1]], [0,plt.xlim()[1]], '--', color='grey')
		plt.ylim(0, 4)
		plt.xlim(0, 4)
		plt.xticks(list(range(5)))
		plt.yticks(list(range(5)))
		plt.xlabel('Spectroscopic redshift')
		plt.ylabel('Photometric redshift')
		plt.text(0.2, plt.ylim()[1]-0.2, "$\\widetilde{\\sigma}_{NMAD} = %.4f$\n$\\chi^2 = %.2f$\n$\\widetilde{\\eta} = %.1f$, $n = %d$\n$\\widetilde{\\sigma}_{tot} = %.3f$" % (nmad, -2 * loglike, smooth_outlier_fraction * 100, len(v), sigma_ad), va='top')
		plt.savefig('photozquality_dist.pdf', bbox_inches='tight')
		plt.savefig('photozquality_dist.png', bbox_inches='tight')
		plt.close()
		
	return r

print('loading data...')
allids = sys.argv[1:]
ids = [i for i in allids if has_zspec(i)]

for i in ids:
	if not has_zspec(i):
		print('skipped', i, os.path.exists('pdz/' + i), os.path.exists('specz/' + i))
		continue
	zphot = get_zphot(i)
	zspec = get_zspec(i)
	data[i] = [zphot, zspec]

print('data loaded')

def optfunc_smooth(args):
	return calc_stats(*args).get('smoothened_chi2', 1e300)

try:
	os.mkdir('smoothened')
except IOError: pass

print()
print('Without smoothing:')
#r = calc_stats(0.04, 1e-5, plot_corr=True)
r = calc_stats(0.04, 1)
print("\tclassic redshift error: sigma_NMAD = %.3f" % r['sigma_NMAD'])
print("\tclassic outlier fraction: eta = %.3f" % r['outlier_fraction'])
print()

x0 = [r['sigma_NMAD']*100, numpy.log10(r['outlier_fraction'])+9]
cons = [lambda x: x[0], lambda x: x[0]]

assert numpy.isfinite(optfunc_smooth(x0))

print('Finding systematic errors ...')
x1_smooth = scipy.optimize.fmin_cobyla(optfunc_smooth, x0=x0, cons=cons, rhobeg=[0.1,0.1], disp=0)
r_smooth = calc_stats(*x1_smooth)
print('Information loss of this method: %.2f bits (lower is better)' % r['information_loss'])
print('Chi^2: %.2f (lower is better)' % r['smoothened_chi2'])



k, wlogoutliers = x1_smooth

print()
print("with optimal systematic errors:")
print("\tsystematic z error: sigma_NMAD = %.3f" % (k/100.))
print("\tsystematic outlier fraction: eta = %.3f" % 10**(wlogoutliers-9))
print()

# plot out smoothened pdzs

calc_stats(*x1_smooth, plot_corr=True, write_individual=True)
