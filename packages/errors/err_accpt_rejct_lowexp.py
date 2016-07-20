
from scipy.optimize import curve_fit
from scipy import stats
import numpy as np
from ..math_f import exp_function
from ..inp import input_params as g


def lin_func(x, a, b):
    '''Linear function.'''
    return a * np.asarray(x) + b


def predband(x, xd, yd, f_vars, conf=0.95):
    """
    Code adapted from Rodrigo Nemmen's post:
    http://astropython.blogspot.com.ar/2011/12/calculating-prediction-band-
    of-linear.html

    Calculates the prediction band of the linear regression model at the
    desired confidence level.

    Clarification of the difference between confidence and prediction bands:

    "The prediction bands are further from the best-fit line than the
    confidence bands, a lot further if you have many data points. The 95%
    prediction band is the area in which you expect 95% of all data points
    to fall. In contrast, the 95% confidence band is the area that has a
    95% chance of containing the true regression line."
    (from http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?
    reg_graphing_tips_linear_regressio.htm)

    Arguments:
    - x: array with x values to calculate the confidence band.
    - xd, yd: data arrays.
    - a, b, c: linear fit parameters.
    - conf: desired confidence level, by default 0.95 (2 sigma)

    References:
    1. http://www.JerryDallal.com/LHSP/slr.htm, Introduction to Simple Linear
    Regression, Gerard E. Dallal, Ph.D.
    """

    alpha = 1. - conf    # Significance
    N = xd.size          # data sample size
    var_n = len(f_vars)  # Number of variables used by the fitted function.

    # Quantile of Student's t distribution for p=(1 - alpha/2)
    q = stats.t.ppf(1. - alpha / 2., N - var_n)

    # Std. deviation of an individual measurement (Bevington, eq. 6.15)
    se = np.sqrt(1. / (N - var_n) * np.sum((yd - lin_func(xd, *f_vars)) ** 2))

    # Auxiliary definitions
    sx = (x - xd.mean()) ** 2
    sxd = np.sum((xd - xd.mean()) ** 2)

    # Predicted values (best-fit model)
    yp = lin_func(x, *f_vars)
    # Prediction band
    dy = q * se * np.sqrt(1. + (1. / N) + (sx / sxd))

    # Return only upper prediction band.
    upb = yp + dy

    return upb


def separate_stars(mag, e_mag, e_col1, be_m, popt_mag, popt_col1):
    '''
    Use the exponential curve obtained to accept/reject stars in the
    magnitude range beyond the (brightest star + be) limit.
    '''

    e_max, be_e = g.er_params[1], g.er_params[3]

    # Initialize empty lists.
    acpt_indx, rjct_indx = [], []

    # Iterate through all stars and accept or reject those beyond
    # the (brightest star + be mag) limit according to the curve
    # obtained for the errors in magnitude and color.
    for st_ind, st_mag in enumerate(mag):

        # Reject stars with at least one error >= e_max.
        if e_mag[st_ind] >= e_max or e_col1[st_ind] >= e_max:
            rjct_indx.append(st_ind)
        else:
            # For stars brighter than the bright end.
            if mag[st_ind] <= be_m:
                # For values in this range accept all stars with both errors
                # < be_e.
                if e_mag[st_ind] < be_e and e_col1[st_ind] < be_e:
                    # Accept star.
                    acpt_indx.append(st_ind)
                else:
                    # Reject star.
                    rjct_indx.append(st_ind)

            else:
                # For the reminder of stars, we check to see if they are
                # located above or below the exp envelope for both errors.
                # If they are above in either one, we reject them, otherwise
                # we accept them.

                # Compare with exponential curve.
                mag_rjct, col1_rjct = False, False
                if e_mag[st_ind] > exp_function.exp_2p(mag[st_ind], *popt_mag):
                    # Reject star.
                    mag_rjct = True

                if e_col1[st_ind] > exp_function.exp_2p(mag[st_ind],
                                                        *popt_col1):
                    # Reject star.
                    col1_rjct = True

                if mag_rjct or col1_rjct:
                    # Reject star.
                    rjct_indx.append(st_ind)
                else:
                    # Accept star.
                    acpt_indx.append(st_ind)

    return acpt_indx, rjct_indx


def main(mag, e_mag, e_col1, be_m):
    '''
    Find the exponential fit to the photometric errors in mag and color
    and reject stars beyond the N*sigma limit.
    '''

    C = g.er_params[4]
    C_val = C if 0. < C <= 1. else 0.95
    # Generate equi-spaced mag values.
    mag_l = np.linspace(mag.min(), mag.max(), 100)

    # Find best fit of data with linear function, previous use of the log()
    # function on the erros.
    for i, err_col in enumerate([np.asarray(e_mag), np.asarray(e_col1)]):
        # Replace bad values.
        err_col[err_col <= 0.] = 0.001
        # Apply natural logarithm to error values. Necessary since the
        # prediction band function expects a *linear* regression model.
        log_err = np.log(err_col)
        # Find best fit using a linear function of the form y=a*x+b.
        slope, intercept, d1, d2, d2 = stats.linregress(mag, log_err)
        # Call function to generate the C_val upper prediction band values.
        upb = predband(mag_l, mag, log_err, [slope, intercept], conf=C_val)
        # Obtain 2P exp function parameters that best fit the upper prediction
        # band. Exponentiate values because the prediction band is linear.
        popt_2p, dummy = curve_fit(exp_function.exp_2p, mag_l, np.exp(upb))
        if i == 0:
            popt_mag = popt_2p
        else:
            popt_col1 = popt_2p

    # Use the fitted curves to identify accepted/rejected stars and store
    # their indexes.
    acpt_indx, rjct_indx = separate_stars(mag, e_mag, e_col1, be_m, popt_mag,
                                          popt_col1)

    err_plot = [popt_mag, popt_col1]

    return acpt_indx, rjct_indx, err_plot
