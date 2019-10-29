
import numpy as np
from scipy.spatial import distance
from .cluster_area import circFrac


def main(clp, x, y, **kwargs):
    """
    Obtain the Radial Density Profile, given a fixed center value.

    Several methods were developed over time, but the standard circular
    rings method seems to be the best one.
    """

    # # Square rings
    # rdp_radii, rdp_points, poisson_error = squareRings(
    #     clp['bin_width'], clp['hist_2d'], clp['bin_cent'])
    # # Per cell distances
    # rdp_radii, rdp_points, poisson_error = radial_prof(
    #     clp['hist_2d'], clp['bin_cent'], clp['bin_width'])
    # # KDE based
    # rdp_radii, rdp_points, poisson_error = kdeRDP(x, y, clp['kde_cent'])

    # Circular rings
    rdp_radii, rdp_points, poisson_error = circRingsRDP(x, y, clp['kde_cent'])

    if rdp_points:
        print("Radial density profile (RDP) calculated")
    else:
        raise ValueError("ERROR: RDP is empty. Check the center coordinates")

    # Add data to dictionary.
    clp['rdp_radii'], clp['rdp_points'], clp['poisson_error'] =\
        rdp_radii, rdp_points, poisson_error
    return clp


def squareRings(bin_width, hist_2d, cent):
    """
    Calculate the density profile by counting the number of stars in the center
    bin first (r aprox width_bin/2 px), then moving to the 8 adjacent bins
    (r aprox width_bin + (width_bin/2) px), then the following 16 (r aprox
    2*width_bin + (width_bin/2) px), then the next 24 (r aprox 3*width_bin +
    (width_bin/2) px) and so forth. These "square rings" have consecutive
    multiples of 8 bins each (if no bin falls outside the range of the x,y
    chart) and an area of N*width_bin^2 [px^2], where N is the number of bins.
    Dividing the number of stars in each "square ring" by this area, we get all
    the "ring densities" for those approximate values of r.
    """
    radii, rdp_points, poisson_error = [], [], []

    # Estimate the central density averaging the 9 central bins (ie:
    # the center bin plus the surrounding 8 in x,y)
    x_c_b, y_c_b = cent
    stars_N, bin_N = 0., 0.
    for i in range(-1, 2):
        for j in range(-1, 2):
            try:
                stars_N += hist_2d[x_c_b + i][y_c_b + j]
                bin_N += 1
            except IndexError:
                pass
    rdp_points.append(stars_N / (bin_N * bin_width ** 2))
    radii.append(bin_width * 1.5)
    poisson_error.append(np.sqrt(stars_N) / (bin_N * bin_width ** 2))

    # Use max x,y length defined in the 2D histogram.
    rdp_length = max(len(hist_2d), len(hist_2d[0]))

    # Iterate through all the bins in the largest dimension.
    for i in range(2, rdp_length):

        # Initialize bin_count for this square ring.
        ring_count, bin_count = 0, 0

        # Iterate through bins in the x dimension for the 2D hist.
        for xindex, xitem in enumerate(hist_2d):
            # Iterate through bins in the y dimension for the 2D hist.
            for yindex, st_in_bin in enumerate(xitem):

                # Bin is in the top row
                if yindex == (y_c_b + i) and abs(xindex - x_c_b) <= i:
                    # Add stars in bin to corresponding ring.
                    ring_count += st_in_bin
                    # Add 1 more bin to the "square ring".
                    bin_count += 1
                # Bin is in the bottom row
                elif yindex == (y_c_b - i) and abs(xindex - x_c_b) <= i:
                    ring_count += st_in_bin
                    bin_count += 1
                # Bin is in the left column
                elif xindex == (x_c_b - i) and abs(yindex - y_c_b) <= (i - 1):
                    ring_count += st_in_bin
                    bin_count += 1
                # Bin is in the right column
                elif xindex == (x_c_b + i) and abs(yindex - y_c_b) <= (i - 1):
                    ring_count += st_in_bin
                    bin_count += 1

        # Break when no more bins are stored in this square ring. This means
        # we reached the border of the frame.
        if bin_count == 0:
            break

        # If no stars are inside this square ring, set value to 1 to avoid a
        # division by zero.
        bin_count = 1 if bin_count == 0 else bin_count
        # The number of bins times the area of each bin gives the area of
        # this square ring.
        area = bin_count * (bin_width ** 2)

        # Calculate density corresponding to "square ring" i
        rdp_points.append(ring_count / area)
        # Obtain the Poisson error bar for each value
        poisson_error.append(np.sqrt(ring_count) / area)

        # Store values for radii to go with the densities obtained above
        # and stored in 'rdp_points'
        radii.append(bin_width / 2. + (bin_width * i))

    return radii, rdp_points, poisson_error


def radial_prof(hist_2d, center, bw):
    """
    https://stackoverflow.com/a/42660542/1391441
    """
    # from .xy_density import cent_bin as center_bin
    # # 2D histogram of the coordinates in the entire frame.
    # hist_2d, xedges, yedges = np.histogram2d(x, y, bins=Nbins)
    # xl, yl = np.ptp(x) / Nbins, np.ptp(y) / Nbins
    # b_area = xl * yl
    # center = center_bin(xedges, yedges, kde_cent)

    xi, yi = np.indices((hist_2d.shape))
    r = np.sqrt((xi - center[0])**2 + (yi - center[1])**2)
    un_radii = np.unique(r)

    radii, rdp, poisson_error = [], [], []
    N_stars, area, rad = 0., 0., []
    for un in un_radii:
        bins_at_dist = hist_2d[r == un]
        N_stars += np.sum(hist_2d[r == un])
        area += len(bins_at_dist) * bw**2
        rad.append(un)
        if N_stars > 10:
            rdp.append(N_stars / area)
            # print(len(rad), min(rad) * bw, max(rad) *bw, np.mean(rad) * bw)
            radii.append(np.mean(rad))
            poisson_error.append(np.sqrt(N_stars) / area)
            # Reset
            N_stars = 0.
            area = 0.
            rad = []

    return np.array(radii) * bw, rdp, poisson_error


def kdeRDP(x, y, kde_cent):
    """
    This  is very dependent on the bandwidth value.
    """
    from scipy.stats import gaussian_kde

    # Distances of all stars to the center of the cluster.
    dist = distance.cdist([kde_cent], np.array([x, y]).T)[0]
    ds = np.sort(dist)

    x0, x1, y0, y1 = min(x), max(x), min(y), max(y)

    dx0, dx1 = abs(kde_cent[0] - x0), abs(kde_cent[0] - x1)
    dy0, dy1 = abs(kde_cent[1] - y0), abs(kde_cent[1] - y1)
    dxy = min(dx0, dx1, dy0, dy1)

    kde = gaussian_kde(ds)  # , bw_method=.01)
    # For the 16 clusters this auto value is in the range (.15, .25)
    print(kde.covariance_factor())
    N_stars = len(ds)

    radii = np.linspace(0., ds[-1], 55)[:-5]

    N_tot = 250000
    rdp_kde, poisson_error = [], []

    fr_area_l = 1.
    for l, h in zip(*[radii[:-1], radii[1:]]):

        if h < dxy:
            area = np.pi * (h**2 - l**2)
        else:
            # From this point on, the rings fall outside the frame.
            fr_area_h = circFrac(kde_cent, h, x0, x1, y0, y1, N_tot=N_tot)
            area = (np.pi * h**2 * fr_area_h) - (np.pi * l**2 * fr_area_l)
            fr_area_l = fr_area_h

        # area = np.pi * (h**2 - l**2)
        N_in_box = kde.integrate_box(l, h) * N_stars
        # print(l, h, N_in_box, area)
        rdp_kde.append(N_in_box / area)
        poisson_error.append(np.sqrt(N_in_box) / area)

    radii = (radii[:-1] + radii[1:]) * .5

    return radii, rdp_kde, poisson_error


def circRingsRDP(x, y, kde_cent, N_tot=250000):
    """
    This is the standard circular rings method, combined with a simple Monte
    Carlo algorithm to estimate the area once the rings overflow the frame.
    """

    # Distances of all stars to the center of the cluster.
    dist = distance.cdist([kde_cent], np.array([x, y]).T)[0]
    # Sorted min to max
    ds = np.sort(dist)

    # Frame limits
    x0, x1, y0, y1 = min(x), max(x), min(y), max(y)
    # Minimum distance to frame border
    dx0, dx1 = abs(kde_cent[0] - x0), abs(kde_cent[0] - x1)
    dy0, dy1 = abs(kde_cent[1] - y0), abs(kde_cent[1] - y1)
    dxy = min(dx0, dx1, dy0, dy1)

    # Discard the last 15 steps. This are the ones with the largest MC variance
    radii = np.linspace(0., ds[-1], 55)[:-15]

    rdp_radii, rdp_points, poisson_error = [], [], []
    rad, stars_in_ring, area_tot, fr_area_l = [], 0., 0., 1.
    for l, h in zip(*[radii[:-1], radii[1:]]):

        # Count stars inside the (r_min=l, r_max=h) ring.
        stars_in_ring += sum((ds >= l) & (ds < h))

        # Area of ring.
        if h < dxy:
            area = np.pi * (h**2 - l**2)
        else:
            # From this point on, the rings fall outside the frame.
            fr_area_h = circFrac(kde_cent, h, x0, x1, y0, y1, N_tot=N_tot)
            area = (np.pi * h**2 * fr_area_h) - (np.pi * l**2 * fr_area_l)
            fr_area_l = fr_area_h

        rad.append(h)
        area_tot += area
        # Use a minimum of 10 stars, else store all values and go to the
        # next radius value.
        if stars_in_ring > 10:
            rdp_radii.append(np.mean(rad))
            rdp_points.append(stars_in_ring / area_tot)
            poisson_error.append(np.sqrt(stars_in_ring) / area_tot)
            # Reset
            rad, stars_in_ring, area_tot = [], 0., 0.

    return rdp_radii, rdp_points, poisson_error
