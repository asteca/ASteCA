

# DEPRECATED 23/03/22

def main(flag_tf, tf_range, data, x_col, y_col):
    """
    Remove stars outside of the (xmin, xmax, ymin, ymax) range.
    """

    if flag_tf:
        xmin, xmax, ymin, ymax = tf_range
        msk = (data[x_col] >= xmin) & (data[x_col] <= xmax) &\
            (data[y_col] >= ymin) & (data[y_col] <= ymax)
        if msk.sum() < 10:
            raise ValueError("ERROR: <10 stars left after frame trimming.")
        data = data[msk]
    else:
        pass

    return data
