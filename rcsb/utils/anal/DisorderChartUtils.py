##
# File: DisorderChartUtils.py
#

import logging
import math

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import numpy as np


logger = logging.getLogger(__name__)


class DisorderChartUtils(object):

    def __init__(self, **kwargs):
        pass

    def __bin_by(self, xIn, yIn, nbins=30):
        """
        Bin y by x.
        Returns the binned "y" values and the left edges of the bins
        """
        x = np.asarray(xIn)
        y = np.asarray(yIn)
        bins = np.linspace(x.min(), x.max(), nbins + 1)
        # To avoid extra bin for the max value
        bins[-1] += .0001

        indicies = np.digitize(x, bins)

        output = []
        for i in range(1, len(bins)):
            output.append(y[indicies == i])

        # Just return the left edges of the bins
        # bins = bins[:-1]

        return output, bins

    def __bin_by_discrete(self, xIn, yIn):
        """
        Bin y by x.
        Returns the binned "y" values and the left edges of the bins
        """
        imin = min(xIn)
        imax = max(xIn)
        bins = range(imin, imax + 2)
        indicies = np.digitize(xIn, bins)

        # x = np.asarray(xIn)
        y = np.asarray(yIn)

        # bins = np.linspace(x.min(), x.max(), nbins + 1)
        # To avoid extra bin for the max value
        # bins[-1] += .0001
        # indicies = np.digitize(x, bins)

        output = []
        for i in range(1, len(bins)):
            output.append(y[indicies == i])

        # Just return the left edges of the bins
        # bins = bins[:-1]

        return output, bins

    def doFigOneTest(self):
        idC = self.__deserializeJson('../results/structures-per-nme.json')
        logger.info("Data set one %d: %r" % (len(idC), idC))
        cD = {}
        for c in idC:
            if c in cD:
                cD[c] += 1
            else:
                cD[c] = 1
        xL = []
        yL = []
        for k in sorted(cD.keys()):
            if k == 0:
                continue
            xL.append(k)
            yL.append(cD[k])
        logger.info("Plot data %r" % list(zip(xL, yL)))

        fig = plt.figure(frameon=False)
        fig.set_size_inches(5, 5)
        ax = fig.gca()
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.bar(xL, yL, label="Figure 1", color='#6894bf')
        plt.ylim(0, 18)
        plt.xscale('log')
        # plt.legend()
        plt.xlabel('Number PDB Structures')
        plt.ylabel('NMEs')
        plt.title('Related PDB Structures by NME')

        for axis in [ax.xaxis]:
            axis.set_major_formatter(ScalarFormatter())
        plt.tight_layout()
        plt.savefig('fig1.png', format='png', dpi=200)

    def doFigTwoTest(self):
        tmL = self.__deserializeJson('../results/deposit-vs-approval-data.json')
        logger.info("Data set two %d: %r" % (len(tmL), tmL))
        cD = {}
        for tm in tmL:
            c = int(tm)
            if c in cD:
                cD[c] += 1
            else:
                cD[c] = 1
        xL = []
        yL = []
        for k in sorted(cD.keys()):
            if k == 0:
                continue
            xL.append(k)
            yL.append(cD[k])
        logger.info("Plot data %r" % list(zip(xL, yL)))

        fig = plt.figure(frameon=False)
        fig.set_size_inches(5, 5)
        ax = fig.gca()
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        plt.bar(xL, yL, label="Figure 1", color='#6894bf')
        plt.ylim(0, 15)
        # plt.xscale('log')
        # plt.legend()
        plt.xlabel('Years before FDA Approval')
        plt.ylabel('NMEs')
        plt.title('Deposition to Approval Interval by NME')

        for axis in [ax.xaxis]:
            axis.set_major_formatter(ScalarFormatter())
        plt.tight_layout()
        plt.savefig('fig2.png', format='png', dpi=200)

    def doHistogramChart(self, valList, **kwargs):
        logger.info("Data set length %d" % (len(valList)))
        plotLabel = kwargs.get('plotLabel', "Figure #")
        xLabel = kwargs.get('xPlotLabel', "X Label")
        yLabel = kwargs.get('yPlotLabel', "Y Label")
        title = kwargs.get('plotTitle', "Title")
        #
        plotPath = kwargs.get('plotPath', "myFigure.png")
        plotFormat = kwargs.get('plotFormat', "png")
        yPlotScale = kwargs.get('yPlotScale', None)
        plotColor = kwargs.get('plotColor', '#6894bf')
        plotGridY = kwargs.get("plotGridY", True)
        #
        logScaleOpt = True if yPlotScale == 'log' else False
        #
        fig = plt.figure(frameon=False)
        fig.set_size_inches(5, 5)

        # plt.hist(valList, bins=20)
        n, bins, patches = plt.hist(x=valList, bins='auto', color=plotColor, alpha=0.7, rwidth=0.85, log=logScaleOpt)
        if plotGridY:
            plt.grid(axis='y', alpha=0.75)
        #
        ax = fig.gca()
        if False:
            ax = fig.gca()
            yPlotMax = kwargs.get('yPlotMax', n.max())
            yPlotMin = kwargs.get('yPlotMin', n.min())
            logger.info("yPlotMin %r yPlotMax = %r" % (yPlotMin, yPlotMax))
            xPlotMax = kwargs.get('xPlotMax', bins.min())
            xPlotMin = kwargs.get('xPlotMin', bins.max())

            logger.info("xPlotMin %r xPlotMax = %r" % (xPlotMin, xPlotMax))

            plt.ylim(yPlotMin, yPlotMax)
            plt.xlim(xPlotMin, xPlotMax)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            for axis in [ax.xaxis]:
                axis.set_major_formatter(ScalarFormatter())
        #
        xPlotMax = kwargs.get('xPlotMax', None)
        xPlotMin = kwargs.get('xPlotMin', None)
        xPlotIncr = kwargs.get('xPlotIncr', None)
        if xPlotMin is not None and xPlotMax is not None and xPlotIncr is not None:
            plt.xlim(xPlotMin, xPlotMax)
            ax.xaxis.set_ticks(np.arange(xPlotMin, xPlotMax, xPlotIncr))

        yPlotMax = kwargs.get('yPlotMax', n.max())
        yPlotMin = kwargs.get('yPlotMin', n.min())
        logger.info("yPlotMin %r yPlotMax %r" % (yPlotMin, yPlotMax))
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
        plotLabel = kwargs.get('plotLabel', "Figure #")
        xLabel = kwargs.get('xPlotLabel', "X Label")
        yLabel = kwargs.get('yPlotLabel', "Y Label")
        title = kwargs.get('plotTitle', "Title")
        #
        plotPath = kwargs.get('plotPath', "myFigure.png")
        plotFormat = kwargs.get('plotFormat', "png")
        logger.info("Integer data set length %d" % (len(iValList)))
        yPlotScale = kwargs.get('yPlotScale', None)
        plotColor = kwargs.get('plotColor', '#6894bf')
        plotGridY = kwargs.get("plotGridY", False)
        #
        cD = {}
        for iVal in iValList:
            c = int(iVal)
            if c in cD:
                cD[c] += 1
            else:
                cD[c] = 1
        xL = []
        yL = []
        for k in sorted(cD.keys()):
            # if k == 0:
            #    continue
            xL.append(k)
            yL.append(cD[k])
        #
        logger.info("Plot data %r" % list(zip(xL, yL)))
        if yPlotScale == 'log':
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
            plt.grid(axis='y', alpha=0.75)
        #
        yPlotMax = kwargs.get('yPlotMax', yMaxVal)
        yPlotMin = kwargs.get('yPlotMin', yMinVal)
        logger.info("yPlotMin %r yPlotMax = %r" % (yPlotMin, yPlotMax))

        xPlotMax = kwargs.get('xPlotMax', xMaxVal)
        xPlotMin = kwargs.get('xPlotMin', xMinVal)
        logger.info("xPlotMin %r xPlotMax = %r" % (xPlotMin, xPlotMax))

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
