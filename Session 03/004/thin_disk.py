"""
Thind disk model

@author: AshCat
"""


def thin_disk(r, R_min , R_max):
    if r>R_min and r<R_max:
        return 1.
    else:
        return 0.
