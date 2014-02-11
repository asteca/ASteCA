"""
@author: gabriel
"""

def get_clust_rad(backg_value, radii, ring_density, width_bins):
    """
    Obtain the value of the cluster's radius by counting the number of points
    that fall within a given interval of the background or lower. If this number
    is equal to the number of points left behind (to the left) then assign the
    radius as the first of those values in the string that shows a "stable"
    behaviour around the background value.
    Iterate increasing the interval around the background until n_left points
    are found or the interval reaches 25% of the total delta or no more
    points are available or a radius of 500 px was reached. In any of these
    last three cases, assign r=500 px.
    """
    
    # Initiate empty list of as many items as different bin widths were used.
    clust_rad = [0]*len(width_bins)
    
    # Use the density profile and the background value obtained with the
    # smallest bin width.
    rd_index, rd_item = 0, ring_density[0]
        
    # Difference between max density value and the background value.
    delta_total = (max(rd_item) - backg_value[rd_index])
    
    # If the difference between the max density value and the background is 
    # less than 3 times the value of the background, raise a flag.
    flag_delta_total = False
    if delta_total < 3*backg_value[rd_index]:
        flag_delta_total = True
    
    # Start i value to cap the number of iterations to 5.
    i = 1
    # Initialize condition to break out of 'while'.
    flag_not_stable = True
    # Iterate until a stable condition is attained or for a maximum of 5 values
    # of delta_backg.
    while flag_not_stable and i < 6:
    
        # i*5% of difference between max density value and background.
        delta_backg = i*5*delta_total/100.
        # Increase value of i for next iteration (if it happens)
        i += 1
        
        # Initialize density values counter (points inside the range determined
        # by the delta around the background) and assign a value to the
        # number of points that determine the 'stabilized' condition.
        in_delta_val, n_left = 0, 4
        # This value stores the number of points that fall outside the delta
        # limits around the background, for the value of 'delta_backg' used.
        outside_val = 0
        
        # Initialize index_rad value.
        index_rad = 0

        # Iterate through all values of star density in this "square ring".
        for index, item in enumerate(rd_item):
            
            # Condition to iterate until at least n_left points close to the
            # background value are found.
            if in_delta_val < n_left:
                
                # If the density value is closer than 'delta_backg' to the
                # background value or lower --> add it.
                # The condition index!=0 is there to avoid counting the first
                # point which sometimes presents a very low density value due
                # to the small width of the bin used.
                if (item - backg_value[rd_index]) <= delta_backg and \
                index !=0:
                    # Augment value of counter.
                    in_delta_val += 1
                    # Store first radius value that falls below the upper delta
                    # limit. This will be the assigned radius if the stable
                    # condition is attained.
                    if in_delta_val == 1:
                        index_rad = index
                else:
                    # If the point falls above the upper delta limit around the
                    # background --> increase value.
                    outside_val += 1
            else:
                # Exit 'for' and 'while' loops if n_left consecutive values were
                # found == "stable" condition.
                flag_not_stable = False
                break
    
    # If radius is > 500px --> re-assign it as 500 and raise a flag.
    flag_rad_500 = False
    # The first condition is there in case that the stable condition was reached
    # with the last item.
    if (in_delta_val < n_left) and flag_not_stable == True:
        # Assign maximum radius value of 500 px.
        clust_rad[rd_index] = 500
    else:
        # Stable condition was attained, assign radius value as the first point
        # that fell inside the delta limits around the background.
        clust_rad[rd_index] = radii[rd_index][index_rad]
        if clust_rad[rd_index] > 500:
            clust_rad[rd_index] = 500
            flag_rad_500 = True
    
    delta_percentage = (i-1)*5
    # Raise a flag if the delta used to get the stable condition is greater
    # than 10%.
    flag_delta = False
    if delta_percentage > 10:
        flag_delta = True
        
    # Raise a flag if the number of points that fall outside the delta limits
    # is higher than the number of points to the left of the one used as the
    # cluster's radius. This indicates a possible variable background.
    # outside_val > (index_rad+1) means that there are at least 2 points that
    # fall outside the delta limits between index_rad and (index_rad+3)
    flag_delta_points = False
    if outside_val > (index_rad+1):
        flag_delta_points = True
    

    return clust_rad, delta_backg, delta_percentage, \
    flag_delta_total, flag_not_stable, flag_rad_500, flag_delta, \
    flag_delta_points