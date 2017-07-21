
import numpy as np
import matplotlib.pyplot as plt
import display_frame


def main(cld, run_mode, coords, **kwargs):
    '''
    Trim frame according to given values of new center and side lengths.
    '''
    if run_mode == 'manual':
        # Unpack dictionary.
        ids, x, y, mags, em, cols, ec = cld['ids'], cld['x'], cld['y'],\
            cld['mags'], cld['em'], cld['cols'], cld['ec']

        # Show full frame plot. Use main magnitude.
        display_frame.main(x, y, mags[0], coords)
        plt.show()

        # Ask to trim frame.
        while True:
            temp_cent, temp_side = [], []
            answer_fra = raw_input('Trim frame? (y/n) ')
            if answer_fra == 'n':
                break
            elif answer_fra == 'y':
                try:
                    print('Input center of new frame (in px).')
                    temp_cent.append(float(raw_input('x: ')))
                    temp_cent.append(float(raw_input('y: ')))
                    print('Input full side length for new frame (in px).')
                    temp_side.append(float(raw_input('x_side: ')))
                    temp_side.append(float(raw_input('y_side: ')))

                    # Empty new lists.
                    ids2, x2, y2 = [], [], []
                    mags2, em2, cols2, ec2 = [[] for _ in mags],\
                        [[] for _ in mags], [[] for _ in cols],\
                        [[] for _ in cols]

                    # Iterate through all stars.
                    for i, star in enumerate(ids):

                        # Check if star is inside new frame boundaries.
                        if abs(temp_cent[0] - x[i]) < temp_side[0] / 2. and \
                                abs(temp_cent[1] - y[i]) < temp_side[1] / 2.:

                            ids2.append(star)
                            x2.append(x[i])
                            y2.append(y[i])
                            for j, m in enumerate(mags):
                                mags2[j].append(m[i])
                                em2[j].append(em[j][i])
                            for j, c in enumerate(cols):
                                cols2[j].append(c[i])
                                ec2[j].append(ec[j][i])

                    # Re-create dictionary.
                    cld['ids'], cld['x'], cld['y'], cld['mags'], cld['em'],\
                        cld['cols'], cld['ec'] = np.array(ids2),\
                        np.array(x2), np.array(y2), np.array(mags2),\
                        np.array(em2), np.array(cols2), np.array(ec2)
                    break
                except Exception:
                    print("Sorry, input is not valid. Try again.")
            else:
                print("Sorry, input is not valid. Try again.")

    return cld
