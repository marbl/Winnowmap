
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "kmers.H"

namespace merylutil::inline kmers::v1 {


double
interpolate(double *h,     //  Array of values
            int32   pp,    //  Base interpolation on h[pp-bo] to h[pp+eo]
            int32   bo,    //   - low  range of interpolation = pp - bo
            int32   eo,    //   - high range of interpolation = pp + eo
            double  x) {   //  Coordinate we're interpolating for
  double y = 0;

  for (int32 ii=pp-bo; ii<=pp+eo; ii++) {
    double t = h[ii];
    for (int32 jj=pp-bo; jj<=pp+eo; jj++) {
      if (ii != jj)
        t *= ((x - jj) / (ii - jj));
    }
    y += t;
  }

  //fprintf(stderr, "Interpolate f(%f) = %f via values from %d to %d\n", x, y, pp-bo, pp+eo); 

  return y;
}


void
displayInterpolation(double *h, char tag) {
  char  N[256] = {0};

  for (int32 ii=6; ii<26; ii++) {
    sprintf(N, "%c%02d", tag, ii);

    FILE *O = fopen(N, "w");
    for (double x=ii-8; x<=ii+8.001; x += 0.025)
      fprintf(O, "  %8.4f %9.4f\n", x, interpolate(h, ii, 3, 3, x));
    fclose(O);
  }
}


//  Interpolate between points to find the 'true' min/max.
//
//  The interpolation uses +-range points around point pp,
//  but we only search for the min/max within +- 1.5 of pp.
//
//  For example, if pp=11 and range=5, we'll build the interpolation
//  so that it goes through points 6, 7, .., 16, but then find the
//  simple min/max between 9.5 and 12.5.  This was done so that the interpolation
//  itself is higher quality (allegedly) and because closely spaced peaks
//  can cause a false peak:
//
//  If we're trying to find the peak 'x',          f
//  but we scan too far, we'll find peak          * *  x
//  'f' instead.                                 *   ** **
//

double
findInterpolatedMinMax(double *h, int32 pp, double range, double step, bool wantmin) {
  double  pmin = pp + range;   double imin = DBL_MAX;
  double  pmax = pp - range;   double imax = DBL_MIN;

  double  pbgn = pp - 1.5;                //  Search, inclusive, between pp-range
  double  pend = pp + 1.5 + step * 1.1;   //  and pp+range; 1.1 to avoid rounding weirdness.

  for (double ptst=pbgn; ptst<pend; ptst += step) {
    double i = interpolate(h, pp, range+1, range+1, ptst);

    if (i < imin)   { pmin = ptst;  imin = i; }
    if (imax < i)   { pmax = ptst;  imax = i; }

    //fprintf(stderr, "interp at pp %d +- %f step %f - ptst %f -> i %f\n", pp, range, step, ptst, i);
  }

  return (wantmin) ? pmin : pmax;
}


//  Find slopes between points.                                          /-*-d0--*
//    d0 is the slope of the segment to the next point             /-*-d1
//    d1 is the slope of the segment to the previous point     *-d2
//    d2 is the slope of the segment before that             ii-2  ii-1  ii+0  ii+1
//
//  Then based on the sign of slopes d0 and d1
//  infer if we're at a min or max or in between.
//
//  Function returns the location of the idx'th peak:
//   - The first maximum is ignored; it is assumed to occur at x=0 and
//     represents noise.
//   - idx=0 is the first minimum, a distinguishing point between noise
//     and real data.
//   - idx=1 is the first maximum.
//   - idx=2 is the second minimum.
//   - idx=3 is the second maximum.
//
double
findExtrema(double *h, uint32 idx, uint32 range) {

  int32  ii=0;                         //  Find the first non-zero data, then
  while ((ii < 100) && (h[ii] == 0))   //  then advance two more so that d2,
    ii++;                              //  d1, d0 (as above) are for valid
  ii += 2;                             //  data points.

  for (int32 ii=3; ii<100; ii++) {
    double d0 = h[ii+1] - h[ii+0];    bool d0n = (d0 < 0);
    double d1 = h[ii+0] - h[ii-1];    bool d1n = (d1 < 0);
    double d2 = h[ii-1] - h[ii-2];    bool d2n = (d2 < 0);

    //  Still decreasing.
    if      ((d1n == true) && (d0n == true)) {
    }

    //  Change from decrease to increase - a min!
    else if ((d1n == true) && (d0n == false))  {
      double m = findInterpolatedMinMax(h, ii, range, 0.025, true);
      //fprintf(stderr, "MIN at value %3d with %f kmers - %7.3f\n", ii, h[ii], m);
      if (idx-- == 0)
        return m;
    }

    //  Change form increase to decrease - a max!
    else if ((d1n == false)  && (d0n == true)) {
      double m = findInterpolatedMinMax(h, ii, range, 0.025, false);
      //fprintf(stderr, "MAX at value %3d with %f kmers - %7.3f\n", ii, h[ii], m);
      if (idx-- == 0)
        return m;
    }

    //  Still increasing.
    else {
    }
  }

  return 0;
}



void
merylHistogram::dumpSmoothed(uint32 n, double *h, double *s) {
  FILE *sm = fopen("orig-vs-smooth.dat", "w");
  for (int32 ii=0; ii<n; ii++)
    fprintf(sm, "%d original %f smooth %f\n", ii, h[ii], s[ii]);
  fclose(sm);
}


void
merylHistogram::dumpDerivs(uint32 n, double *h) {
  int32  ii=0;                       //  Find the first non-zero data, then

  while ((ii < n) && (h[ii] == 0))   //  then advance two more so that d2,
    ii++;                            //  d1, d0 (as above) are for valid
  ii += 2;                           //  data points.

  FILE *sm = fopen("derivitives.dat", "w");
  for (int32 ii=3; ii<n; ii++) {
    double d0 = h[ii+1] - h[ii+0];    bool d0n = (d0 < 0);
    double d1 = h[ii+0] - h[ii-1];    bool d1n = (d1 < 0);
    double d2 = h[ii-1] - h[ii-2];    bool d2n = (d2 < 0);

    fprintf(sm, "%d  %f %f %f  %f %f\n", ii, d2, d1, d0, d1-d2, d0-d1);
  }
  fclose(sm);
}


void
merylHistogram::dumpMirror(uint32 n, double *h) {   //  This is, I think, attempting
  double *mirror = new double [n];                  //  to show the asymmetry around
  double  mp = _coveragePeaks[1];                   //  the first peak.

  for (int32 xx=0; xx<n; xx++)
    mirror[xx] = 3.33;

  for (int32 xx=(int32)floor(_coveragePeaks[0]*10); xx<(int32)floor(10*mp); xx++)
    mirror[xx] = (interpolate(h, (int32)round(       xx/10.0), 2, 2,        xx/10.0) -
                  interpolate(h, (int32)round(2*mp - xx/10.0), 2, 2, 2*mp - xx/10.0));

  FILE *sm = fopen("mirror.dat", "w");
  for (int32 xx=0; xx<n; xx++)
    fprintf(sm, "%.1f %f %f %f\n",
            xx/10.0,
            mirror[xx],
            interpolate(h, (int32)round(       xx/10.0), 2, 2,        xx/10.0),
            interpolate(h, (int32)round(2*mp - xx/10.0), 2, 2, 2*mp - xx/10.0));
  fclose(sm);

  delete [] mirror;

  double halfPeak = findExtrema(mirror, 1, 4);
  fprintf(stderr, "halfPeak %f\n", halfPeak);
}



void
merylHistogram::computePloidyPeaks(bool debugPloidy) {

  //
  //  Compute a Weierstrass transform -- convolve with a Gaussian kernel with
  //  specific parameters -- to smooth the counts.
  //
  //  The main parmeter below is 't'.
  //    t =   0.5  - 
  //    t =   1.00 - default
  //    t =   2.00 - low peaks shift +1; -12 to 12
  //    t =   5.00 - peaks exist, but the low ones shift minCov from ~8 to ~12; -20 to 20
  //    t = 100.0  - smoothed out all the peaks
  //

  double  t = 0.5;       //  Smoothing amount.
  double  ws[101];       //  Convolution kernel data.
  double *w = ws + 50;   //  Offset so w[-1] returns valid data.

  int32   iim = 9;       //  Width of the transform.

  double  t4s = 1.0 / sqrt(4 * M_PI * t);
  double  t4n = -4 * t;

  for (int32 ii=-iim; ii<=iim; ii++)
    w[ii] = t4s * exp(ii * ii / t4n);

  //
  //  Load the histogram data, converting from value[ii] and occurrences[ii]
  //  to a simple array indexed by value.
  //

  double  *h = new double [1024];
  double  *s = new double [1024];

  for (int32 ii=0; ii<1024; ii++)
    h[ii] = s[ii] = 0;

  for (int32 ii=0; ii<histogramLength(); ii++) {
    int64 hv = histogramValue(ii);
    int64 ho = histogramOccurrences(ii);

    if (hv < 1024)  //  The histogram is stored as value[ii] and occurrences[ii],
      h[hv] = ho;   //  where ii is just an index into the list.  The list does
    else            //  NOT contain empty cells (so we must clear h[] first) and
      break;        //  ii is NOT the value.
  }

  //
  //  Convolve to generate the smoothed histogram.  Fast enough, but could
  //  definitely be faster.
  //

  for (int32 hi=0; hi<1024; hi++) {
    for (int64 wi=-iim; wi<=iim; wi++) {  //  Index into Weierstrass values
      int64  whi = hi + wi;               //  Index into histogram occ for that Weierstrass value

      if      (whi < 1)     s[hi] += 0      * w[wi];
      else if (whi < 1024)  s[hi] += h[whi] * w[wi];
      else                  s[hi] += 0;
    }
  }

  if (debugPloidy)
    dumpSmoothed(1024, h, s);

  //
  //  Compute the noise/genomic threshold, and 1x, 2x, 3x and 4x peaks.
  //
  //  Use the original un-smoothed data to find the trough between noise and
  //  genomic - smoothed data shifts the trough to the right because the
  //  noise peak is so large - but use the smoothed data to find the ploidy
  //  peaks - unsmoothed data occasionally has local maxima that mess it up.
  //
  _coveragePloidy[0] = 0.0;
  _coveragePeaks[0]  = findExtrema(h, 0, 3);

  for (uint32 ii=1; ii<9; ii++) {
    _coveragePloidy[ii] = ii;
    _coveragePeaks[ii]  = findExtrema(s, 2 * ii - 1, 4);
  }

  //
  //  Some histograms confuse 1/2x and 1x.  Here we try to decide if peak ii=1
  //  represents 1/2x or 1x.
  //   - first find the max peak.
  //   - then, if that max peak is not at ploidy=1, divide all
  //     ploidies by 2 to make the mak peak be at ploidy 1.
  //

  uint32  maxPeak  = 1;
  double  maxPeakX = _coveragePeaks[1];
  double  maxPeakY = interpolate(h, (int32)round(maxPeakX), 3, 3, maxPeakX);

  for (uint32 ii=1; ii<9; ii++) {
    double mx = _coveragePeaks[ii];
    double my = interpolate(h, (int32)round(mx), 3, 3, mx);

    if (maxPeakY < my) {
      maxPeak  = ii;
      maxPeakX = mx;
      maxPeakY = my;
    }
  }

  if (debugPloidy)
    fprintf(stderr, "maxPeak #%d of %f at X=%f\n", maxPeak, maxPeakY, maxPeakX);

  while (--maxPeak > 0)             //  Divide poidy by 2 to normalize
    for (uint32 ii=0; ii<9; ii++)   //  it so that 1x coverage is max.
      _coveragePloidy[ii] /= 2.0;   //

  if (debugPloidy)
    dumpMirror(1024, h);
  if (debugPloidy)
    dumpDerivs(1024, h);

  //  Cleanup and quit.

  delete [] s;
  delete [] h;
}



}  //  namespace merylutil::kmers::v1
