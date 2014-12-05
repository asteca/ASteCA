"""
@author: gabriel
"""


def spiral():
    '''
    Create a list containing lists of two items. These two items are the
    x and y coordinates of a spiral array starting from [0,0] and moving
    right and bottom to obtain finally: [[0,0], [1,0], [1,-1], [0,1], ...]
    '''

    # Total number of spiral rings to be generated
    max_rings = 170
    # Initialize the list that will contain the spiral coordinates.
    spir_lst = []
    # Append center bin
    spir_lst.append([0, 0])

    for ring in range(1, max_rings):

        # Populate the right column first, then the bottom row, then the
        # left column and finally the top row.

        # Right column
        xi = ring
        for yi in range((ring - 1), -ring, -1):
            spir_lst.append([xi, yi])

        # Bottom row
        yi = -ring
        for xi in range(ring, -(ring + 1), -1):
            spir_lst.append([xi, yi])

        # Left column
        xi = -ring
        for yi in range(-(ring - 1), ring, 1):
            spir_lst.append([xi, yi])

        # Top row
        yi = ring
        for xi in range(-ring, (ring + 1), 1):
            spir_lst.append([xi, yi])

    return spir_lst