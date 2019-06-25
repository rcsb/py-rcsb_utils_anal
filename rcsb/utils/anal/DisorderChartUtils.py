##
# File: DisorderChartUtils.py
#

import logging
import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, ScalarFormatter  # pylint: disable=ungrouped-imports

matplotlib.use("Agg")


logger = logging.getLogger(__name__)


class DisorderChartUtils(object):
    def __init__(self, **kwargs):
        pass

    def __binBy(self, xIn, yIn, nbins=30):
        """
        Bin y by x.
        Returns the binned "y" values and the left edges of the bins
        """
        xV = np.asarray(xIn)
        yV = np.asarray(yIn)
        bins = np.linspace(xV.min(), xV.max(), nbins + 1)
        # To avoid extra bin for the max value
        bins[-1] += 0.0001

        indicies = np.digitize(xV, bins)

        output = []
        for i in range(1, len(bins)):
            output.append(yV[indicies == i])

        # Just return the left edges of the bins
        # bins = bins[:-1]

        return output, bins

    def __binByDiscrete(self, xIn, yIn):
        """
        Bin y by x.
        Returns the binned "y" values and the left edges of the bins
        """
        imin = min(xIn)
        imax = max(xIn)
        bins = range(imin, imax + 2)
        indicies = np.digitize(xIn, bins)

        # x = np.asarray(xIn)
        yV = np.asarray(yIn)

        # bins = np.linspace(x.min(), x.max(), nbins + 1)
        # To avoid extra bin for the max value
        # bins[-1] += .0001
        # indicies = np.digitize(x, bins)

        output = []
        for i in range(1, len(bins)):
            output.append(yV[indicies == i])

        # Just return the left edges of the bins
        # bins = bins[:-1]

        return output, bins

    def doHistogramChart(self, valList, **kwargs):
        logger.info("Data set length %d", len(valList))
        # plotLabel = kwargs.get("plotLabel", "Figure #")
        xLabel = kwargs.get("xPlotLabel", "X Label")
        yLabel = kwargs.get("yPlotLabel", "Y Label")
        title = kwargs.get("plotTitle", "Title")
        #
        plotPath = kwargs.get("plotPath", "myFigure.png")
        plotFormat = kwargs.get("plotFormat", "png")
        yPlotScale = kwargs.get("yPlotScale", None)
        plotColor = kwargs.get("plotColor", "#6894bf")
        plotGridY = kwargs.get("plotGridY", True)
        plotNumBins = kwargs.get("plotNumBins", 100)
        #
        logScaleOpt = True if yPlotScale == "log" else False
        #
        #
        fig = plt.figure(frameon=False)
        ax = fig.gca()
        fig.set_size_inches(5, 5)
        #
        xPlotMax = kwargs.get("xPlotMax", None)
        xPlotMin = kwargs.get("xPlotMin", None)
        xPlotIncr = kwargs.get("xPlotIncr", None)
        binBoundaries = None
        if xPlotMin is not None and xPlotMax is not None and xPlotIncr is not None:
            plt.xlim(xPlotMin, xPlotMax)
            ax.xaxis.set_ticks(np.arange(xPlotMin, xPlotMax, xPlotIncr))
            if plotNumBins is not None:
                binBoundaries = np.linspace(xPlotMin, xPlotMax, plotNumBins)
        #
        # plt.hist(valList, bins=20)
        if binBoundaries is not None:
            logger.info("Using (%d) bins", plotNumBins)
            num, bins, _ = plt.hist(x=valList, bins=binBoundaries, color=plotColor, alpha=0.7, rwidth=0.85, log=logScaleOpt)
        else:
            num, bins, _ = plt.hist(x=valList, bins="auto", color=plotColor, alpha=0.7, rwidth=0.85, log=logScaleOpt)
        #
        if plotGridY:
            plt.grid(axis="y", alpha=0.75)

        pltOpt = False
        if pltOpt:
            ax = fig.gca()
            yPlotMax = kwargs.get("yPlotMax", num.max())
            yPlotMin = kwargs.get("yPlotMin", num.min())
            logger.info("yPlotMin %r yPlotMax = %r", yPlotMin, yPlotMax)
            xPlotMax = kwargs.get("xPlotMax", bins.min())
            xPlotMin = kwargs.get("xPlotMin", bins.max())

            logger.info("xPlotMin %r xPlotMax = %r", xPlotMin, xPlotMax)

            plt.ylim(yPlotMin, yPlotMax)
            plt.xlim(xPlotMin, xPlotMax)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            for axis in [ax.xaxis]:
                axis.set_major_formatter(ScalarFormatter())
        #
        yPlotMax = kwargs.get("yPlotMax", num.max())
        yPlotMin = kwargs.get("yPlotMin", num.min())
        logger.info("yPlotMin %r yPlotMax %r", yPlotMin, yPlotMax)
        plt.ylim(yPlotMin, yPlotMax)
        #
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # plt.legend()
        # plt.text(23, 45, r'$\mu=15, b=3$')
        plt.xlabel(xLabel)
        plt.ylabel(yLabel)
        plt.title(title)
        #
        plt.tight_layout()
        plt.savefig(plotPath, format=plotFormat, dpi=200)

    def doIntegerBarChart(self, iValList, **kwargs):
        plotLabel = kwargs.get("plotLabel", "Figure #")
        xLabel = kwargs.get("xPlotLabel", "X Label")
        yLabel = kwargs.get("yPlotLabel", "Y Label")
        title = kwargs.get("plotTitle", "Title")
        #
        plotPath = kwargs.get("plotPath", "myFigure.png")
        plotFormat = kwargs.get("plotFormat", "png")
        logger.info("Integer data set length %d", len(iValList))
        yPlotScale = kwargs.get("yPlotScale", None)
        plotColor = kwargs.get("plotColor", "#6894bf")
        plotGridY = kwargs.get("plotGridY", False)
        #
        cD = {}
        for iVal in iValList:
            cV = int(iVal)
            if cV in cD:
                cD[cV] += 1
            else:
                cD[cV] = 1
        xL = []
        yL = []
        for k in sorted(cD.keys()):
            # if k == 0:
            #    continue
            xL.append(k)
            yL.append(cD[k])
        #
        logger.info("Plot data %r", list(zip(xL, yL)))
        if yPlotScale == "log":
            yL = [math.log10(y) for y in yL]
        yMaxVal = max(yL)
        yMinVal = min(yL)
        #
        xMinVal = min(xL)
        xMaxVal = min(xL)
        #
        fig = plt.figure(frameon=False)
        fig.set_size_inches(5, 5)
        ax = fig.gca()

        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        plt.bar(xL, yL, label=plotLabel, color=plotColor)
        if plotGridY:
            plt.grid(axis="y", alpha=0.75)
        #
        yPlotMax = kwargs.get("yPlotMax", yMaxVal)
        yPlotMin = kwargs.get("yPlotMin", yMinVal)
        logger.info("yPlotMin %r yPlotMax = %r", yPlotMin, yPlotMax)

        xPlotMax = kwargs.get("xPlotMax", xMaxVal)
        xPlotMin = kwargs.get("xPlotMin", xMinVal)
        logger.info("xPlotMin %r xPlotMax = %r", xPlotMin, xPlotMax)

        plt.ylim(yPlotMin, yPlotMax)
        plt.xlim(xPlotMin, xPlotMax)
        # if yPlotScale:
        #    plt.yscale(yPlotScale)
        # plt.xscale('log')
        # plt.legend()
        plt.xlabel(xLabel)
        plt.ylabel(yLabel)
        plt.title(title)

        for axis in [ax.xaxis]:
            axis.set_major_formatter(ScalarFormatter())
        plt.tight_layout()
        plt.savefig(plotPath, format=plotFormat, dpi=200)
