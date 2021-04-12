import numpy as np

def list_to_array(list):
    """
    :param list:
    :return:
    """

    arr = np.zeros(len(list))
    for i in range(len(list)):
        arr[i] = list[i]

    return arr