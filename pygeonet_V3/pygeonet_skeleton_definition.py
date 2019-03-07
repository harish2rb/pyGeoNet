import numpy as np


def compute_skeleton_by_single_threshold(inputArray, threshold):
    """
    compute the skeleton based on a single threshold
    Skeleton by thresholding one grid measure e.g. flow or curvature

    :param inputArray: the input array
    :param threshold: the threshold value to be used on the input array
    :return: a binary numpy array called the skeletonarray
    """
    skeletonArray = np.zeros((inputArray.shape))
    skeletonArray[np.where(inputArray > threshold)] = 1
    return skeletonArray


def compute_skeleton_by_dual_threshold(inputArray1, inputArray2, threshold1, threshold2):
    """
    Compute Skeleton by using thresholds on two grid measures e.g. flow and curvature

    :param inputArray1: the input array 1
    :param inputArray2: the input array 2
    :param threshold1: the threshold value for input array 1
    :param threshold2: the threshold value for input array 2
    :return: a binary numpy array called the skeletonarray
    """
    skeletonArray = np.zeros((inputArray1.shape))
    mask1 = np.where(inputArray1 > threshold1, 1, False)
    mask2 = np.where(inputArray2 > threshold2, 1, False)
    skeletonArray = mask1 * mask2
    return skeletonArray