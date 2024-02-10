import numpy as np

def robust_mean(y, cut):
  """Computes a robust mean estimate in the presence of outliers.
  Args:
    y: 1D numpy array. Assumed to be normally distributed with outliers.
    cut: Points more than this number of standard deviations from the median are
      ignored.
  Returns:
    mean: A robust estimate of the mean of y.
    mean_stddev: The standard deviation of the mean.
    mask: Boolean array with the same length as y. Values corresponding to
        outliers in y are False. All other values are True.
  """
  # First, make a robust estimate of the standard deviation of y, assuming y is
  # normally distributed. The conversion factor of 1.4826 takes the median
  # absolute deviation to the standard deviation of a normal distribution.
  # See, e.g. https://www.mathworks.com/help/stats/mad.html.
  absdev = np.abs(y - np.median(y))
  sigma = 1.4826 * np.median(absdev)

  # If the previous estimate of the standard deviation using the median absolute
  # deviation is zero, fall back to a robust estimate using the mean absolute
  # deviation. This estimator has a different conversion factor of 1.253.
  # See, e.g. https://www.mathworks.com/help/stats/mad.html.
  if sigma < 1.0e-24:
    sigma = 1.253 * np.mean(absdev)

  # Identify outliers using our estimate of the standard deviation of y.
  mask = absdev <= cut * sigma

  # Now, recompute the standard deviation, using the sample standard deviation
  # of non-outlier points.
  sigma = np.std(y[mask])

  # Compensate the estimate of sigma due to trimming away outliers. The
  # following formula is an approximation, see
  # http://w.astro.berkeley.edu/~johnjohn/idlprocs/robust_mean.pro.
  sc = np.max([cut, 1.0])
  if sc <= 4.5:
    sigma /= (-0.15405 + 0.90723 * sc - 0.23584 * sc**2 + 0.020142 * sc**3)

  # Identify outliers using our second estimate of the standard deviation of y.
  mask = absdev <= cut * sigma

  # Now, recompute the standard deviation, using the sample standard deviation
  # with non-outlier points.
  sigma = np.std(y[mask])

  # Compensate the estimate of sigma due to trimming away outliers.
  sc = np.max([cut, 1.0])
  if sc <= 4.5:
    sigma /= (-0.15405 + 0.90723 * sc - 0.23584 * sc**2 + 0.020142 * sc**3)

  # Final estimate is the sample mean with outliers removed.
  mean = np.mean(y[mask])
  mean_stddev = sigma / np.sqrt(len(y) - 1.0)

  return mean, mean_stddev, mask

def lcbin(time2, lc, nbins, method = 'robustmean', sigmacutin = 3, errorin = None):
  time = time2 - np.min(time2)  
  binsize = float(np.max(time) - np.min(time))/float(nbins)
  
  fbin = np.zeros(nbins) + np.nan #replicate(!Values.F_NAN, nbins)
  ebin = np.copy(fbin) #replicate(!Values.F_NAN, nbins)
  tbin = (np.arange(nbins) + 0.5 ) * binsize
  tbin = tbin + np.min(time2)

  goodindmask = np.zeros(len(time2)) + np.nan
  
  for i in range(nbins):
    w = np.argwhere(np.logical_and(time >= i*binsize, time < (i+1)*binsize))
    if len(w) > 0:   
      if method == 'robustmean':
        thisfbin,thissigmafbin,thismask  = robust_mean(lc[w],sigmacutin)
        thismask = np.ravel(thismask)
        fbin[i] = thisfbin
        goodindmask[np.ravel(w)] = thismask
        goodind = np.argwhere(thismask)
        if len(w) > 1:
          if np.sum(thismask) > 1: # gt 1 then begin
            ebin[i] = np.std(lc[w[goodind]], ddof=1) / np.sqrt(np.sum(thismask))
            if errorin is not None: 
                fbin[i] = np.sum(lc[w[goodind]] / errorin[w[goodind]]**2) / np.sum(1.0 / errorin[w[goodind]]**2)
                ebin[i] = np.sqrt(1.0 / np.sum(1.0 / errorin[w[goodind]]**2))
          if np.sum(thismask) == 1: #(goodind) eq 1 then begin
            ebin[i] = np.inf
      
  return tbin,fbin, ebin, goodindmask
