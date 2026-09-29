import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import prepare_pygeonet_defaults as defaults


def raster_plot(Array, title):
    """
    Plot raster image

    :param Array: The input Array
    :param title: The raster plot title
    :return: A matplotlib figure
    """
    if not hasattr(defaults, 'figureNumber'):
        defaults.figureNumber = 0
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(Array, cmap=cm.coolwarm)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(title)
    plt.colorbar()
    if defaults.doPlot == 1:
        plt.show()


def raster_point_plot(Array, PointsList, title, color=cm.coolwarm, point_style='go'):
    """
     Plot raster image with points on top of it.

    :param Array: The input raster to plot
    :param PointsList: The points list to plot on top of raster
    :param title: The plot title
    :param color: the color of the raster image
    :param point_style: the point style
    :return: A matplotlib figure
    """
    if not hasattr(defaults, 'figureNumber'):
        defaults.figureNumber = 0
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(Array, cmap=color)
    plt.plot(PointsList[1], PointsList[0], point_style)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(title)
    plt.colorbar()
    if defaults.doPlot == 1:
        plt.show()


def geodesic_contour_plot(geodesicDistanceArray, title):
    """
    Plot the geodesic distance contout plot

    :param geodesicDistanceArray: The input geodesic Distance Array
    :param title: The plot title
    :return: A matplotlib figure
    """
    if not hasattr(defaults, 'figureNumber'):
        defaults.figureNumber = 0
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(np.log10(geodesicDistanceArray), cmap=cm.coolwarm)
    plt.contour(geodesicDistanceArray, 140, cmap=cm.coolwarm)
    plt.title(title)
    plt.colorbar()
    if defaults.doPlot == 1:
        plt.show()


def channel_plot(flowDirectionsArray, geodesicPathsCellList,
                 xx, yy, title, color=cm.coolwarm,
                 point_style='go', line_style='k-'):
    """
    plot the channels

    :param flowDirectionsArray: The flow directions array
    :param geodesicPathsCellList: The geodesic paths cell list
    :param xx: the x locations of the points
    :param yy: the y locations of the points
    :param title: the plot title
    :param color: the color scale
    :param point_style: the point style
    :param line_style: the line style
    :return: A matplotlib figure
    """
    if not hasattr(defaults, 'figureNumber'):
        defaults.figureNumber = 0
    defaults.figureNumber = defaults.figureNumber + 1
    plt.figure(defaults.figureNumber)
    plt.imshow(flowDirectionsArray, cmap=color)
    for pp in range(0, len(geodesicPathsCellList)):
        plt.plot(geodesicPathsCellList[pp][1, :], geodesicPathsCellList[pp][0, :], line_style)
    plt.plot(xx, yy, point_style)
    plt.title(title)
    plt.colorbar()
    if defaults.doPlot == 1:
        plt.show()
