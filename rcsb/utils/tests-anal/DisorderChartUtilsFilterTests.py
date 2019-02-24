##
# File:    DisorderChartUtilsFilterTests.py
# Author:  J. Westbrook
# Date:   23-Feb-2019
#
# Updates:
#
##
"""
Tests for plotting tools for disorder data analysis

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"


import logging
import os
import time
import unittest


from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.anal.DisorderChartUtils import DisorderChartUtils
from rcsb.utils.anal.TaxonomyUtils import TaxonomyUtils

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class DisorderChartUtilsFilterTests(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super(DisorderChartUtilsFilterTests, self).__init__(methodName)
        self.__verbose = True

    def setUp(self):
        #
        #
        self.__mockTopPath = os.path.join(TOPDIR, 'rcsb', 'mock-data')
        configPath = os.path.join(TOPDIR, 'rcsb', 'mock-data', 'config', 'dbload-setup-example.yml')
        configName = 'site_info'
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        #
        self.__workPath = os.path.join(HERE, 'test-output')
        self.__instanceSavePath = os.path.join(HERE, 'test-data', 'instance-data-dictionary.pic')
        self.__entryResDataSavePath = os.path.join(HERE, 'test-data', 'entry-data-dictionary-with-res.pic')
        #
        self.__plotGapCount = os.path.join(HERE, 'test-output', 'gap-count-fig-all.png')
        self.__plotGapLength = os.path.join(HERE, 'test-output', 'gap-length-fig-all.png')
        self.__plotOwabSegmentCount = os.path.join(HERE, 'test-output', 'owab-segment-count-fig-all.png')
        self.__plotOwabSegmentLength = os.path.join(HERE, 'test-output', 'owag-segment-length-fig-all.png')
        self.__plotCoverageRefDb = os.path.join(HERE, 'test-output', 'coverage-refdb-fig-all.png')
        self.__plotCoverageSample1 = os.path.join(HERE, 'test-output', 'coverage-sample-fig-all-1.png')
        self.__plotCoverageSample2 = os.path.join(HERE, 'test-output', 'coverage-sample-fig-all-2.png')
        #
        self.__mU = MarshalUtil()
        #
        self.__entryResD = self.__mU.doImport(self.__entryResDataSavePath, format="pickle")
        #
        self.__tU = TaxonomyUtils()
        self.__startTime = time.time()
        logger.debug("Starting %s at %s" % (self.id(),
                                            time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)\n" % (self.id(),
                                                              time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                              endTime - self.__startTime))

    def __writeLegend(self, fp, message):
        try:
            fn, ext = os.path.splitext(fp)
            lfp = fn + '.txt'
            #
            with open(lfp, 'w') as ofh:
                ofh.write("%s" % message)
            return True
        except Exception as e:
            logger.exception("Failing with %s %s" % (fp, str(e)))
        return False

    def __filterEntry(self, entryId, **kwargs):
        """
            'homomeric protein' 'Single protein entity'
            'heteromeric protein' 'Multiple protein entities'
            'DNA' 'DNA entity/entities only'
            'RNA' 'RNA entity/entities only'
            'NA-hybrid' 'DNA/RNA hybrid entity/entities only'
            'protein/NA' 'Both protein and nucleic acid polymer entities'
            'DNA/RNA' 'Both DNA and RNA polymer entities'
            'oligosaccharide' 'One of more oligosaccharide entities'
            'protein/oligosaccharide' 'Both protein and oligosaccharide entities'
            'NA/oligosaccharide' 'Both NA and oligosaccharide entities'

        """
        resLimit = kwargs.get("resLimit", None)
        polymerComposition = kwargs.get("polymerComposition", None)
        #
        #
        ok1 = True
        if resLimit is not None:
            ok1 = self.__entryResD[entryId]['ls_d_res_high'] <= resLimit

        ok2 = True
        if polymerComposition is not None:
            ok2 = polymerComposition == self.__entryResD[entryId]['polymer_composition']

        return ok1 and ok2

    def __filterEntity(self, entityD, **kwargs):
        taxIdTarget = kwargs.get('taxId', None)
        taxClassTarget = kwargs.get('taxType', None)
        #
        ok1 = False
        if taxIdTarget is not None:
            taxId = entityD['ncbi_taxonomy_id'] if 'ncbi_taxonomy_id' in entityD else None
            ok1 = True if taxId == taxIdTarget else False
        #
        ok2 = False
        if taxClassTarget is not None:
            taxId = entityD['ncbi_taxonomy_id'] if 'ncbi_taxonomy_id' in entityD else None
            if taxId is None:
                ok2 = False
            elif taxClassTarget == 'Eukaryota':
                ok2 = self.__tU.isEukaryota(taxId)
            elif taxClassTarget == 'Bacteria':
                ok2 = self.__tU.isBacteria(taxId)
            elif taxClassTarget == 'Virus':
                ok2 = self.__tU.isVirus(taxId)
            elif taxClassTarget == 'Archaea':
                ok2 = self.__tU.isArchaea(taxId)
            elif taxClassTarget == 'Unclassified':
                ok2 = self.__tU.isUnclassified(taxId)
            elif taxClassTarget == 'Other':
                ok2 = self.__tU.isOther(taxId)
            else:
                ok2 = False
        #
        return ok1 or ok2

    def testViewData(self):
        """ Test case - extract entity instance data -

        """
        try:
            entryD = self.__mU.doImport(self.__instanceSavePath, format="pickle")
            for entryId in entryD:
                for entityId, eD in entryD[entryId]['selected_polymer_entities'].items():

                    analD = eD['anal_instances'] if 'anal_instances' in eD else {}
                    for asymId, aD in analD.items():
                        logger.info("entryId %s entityId %s asymId %s analD: %r" % (entryId, entityId, asymId, aD))

        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testViewEntryData(self):
        """ Test case - extract entity instance data -

        """
        try:
            entryD = self.__mU.doImport(self.__instanceSavePath, format="pickle")
            for entryId, enD in entryD.items():
                logger.debug("Entry keys %r" % list(self.__entryResD[entryId].keys()))
                if self.__filterEntry(entryId, resLimit=1.2, polymerComposition='heteromeric protein'):
                    logger.info("%s resolution %r pc %r " % (entryId, self.__entryResD[entryId]['ls_d_res_high'], self.__entryResD[entryId]['polymer_composition']))
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testViewEntityData(self):
        """ Test case - extract entity instance data -

        """
        try:
            iCountH = 0
            iCountB = 0
            entryD = self.__mU.doImport(self.__instanceSavePath, format="pickle")
            for entryId, enD in entryD.items():
                if not self.__filterEntry(entryId, resLimit=1.2):
                    continue
                resV = self.__entryResD[entryId]['ls_d_res_high']
                for entityId, eD in entryD[entryId]['selected_polymer_entities'].items():
                    if self.__filterEntity(eD, taxId=9606):
                        taxId = eD['ncbi_taxonomy_id'] if 'ncbi_taxonomy_id' in eD else None
                        logger.debug("%s %s selected %.2f taxId %r" % (entryId, entityId, resV, taxId))
                        iCountH += 1
                    if self.__filterEntity(eD, taxType='Bacteria'):
                        sn = eD['ncbi_scientific_name'] if 'ncbi_scientific_name' in eD else None
                        logger.info("%s %s selected %.2f sn %r" % (entryId, entityId, resV, sn))
                        iCountB += 1
            #
            logger.info("Homo sapiens %d bacteria %d" % (iCountH, iCountB))
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    ###
    ###
    def __getPlotFileName(self, fnBase, resLimit=3.5, polymerComposition=None, taxIdTarget=None, taxType=None, format="png"):
        try:
            fnQ = fnBase + "res-%.1f-" % resLimit
            if polymerComposition is not None:
                tS = polymerComposition.replace(' ', '-')
                fnQ += tS + '-'
            if taxIdTarget is not None:
                fnQ += 'taxid-%d' % taxIdTarget
            if taxType is not None:
                fnQ += 'tax-%s' % taxType.lower()
            fnQ += "." + format
        except Exception as e:
            logger.exception("Failing with % s" % str(e))
            fnQ = None
        #
        return fnQ

    def __getLegend(self, legendBase, numEntries=0, resLimit=3.5, polymerComposition=None, taxIdTarget=None, taxType=None):
        """
         t= "Gap count statistics for all (%d) protein sequences (X-ray resolution limit < 3.5 Angstoms) "
        """
        legend = None
        try:
            tS = legendBase
            if polymerComposition is not None:
                tS += " (%s)" % polymerComposition

            if taxIdTarget is not None:
                tS += " with TaxId=%d" % taxIdTarget

            if taxType is not None:
                tS += " with taxonomy %s" % taxType

            tS += " in %d X-ray structures with resolution limit < %.1f Angstroms" % (numEntries, resLimit)
            legend = tS
        except Exception as e:
            logger.exception("Failing with % s" % str(e))
            legend = None
        #
        return legend

    def testPlotFilteredGapData(self):
        try:
            for taxType in ['Eukaryota', 'Bacteria', 'Virus', 'Archaea']:
                self.__filterViewGapData(resLimit=3.5, polymerComposition=None, taxIdTarget=None, taxType=taxType)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def __filterViewGapData(self, resLimit=3.5, polymerComposition=None, taxIdTarget=None, taxType=None):
        """ Test case - protein sequence instance gaps and widths all cases
        """
        try:
            entryD = self.__mU.doImport(self.__instanceSavePath, format="pickle")
            gapCountList = []
            gapLengthList = []
            #
            numEntries = 0
            for entryId, enD in entryD.items():
                if not self.__filterEntry(entryId, resLimit=resLimit, polymerComposition=polymerComposition):
                    continue
                numEntries += 1
                logger.debug("%s resolution %r pc %r " % (entryId, self.__entryResD[entryId]['ls_d_res_high'], self.__entryResD[entryId]['polymer_composition']))

                for entityId, eD in entryD[entryId]['selected_polymer_entities'].items():
                    if not self.__filterEntity(eD, taxId=taxIdTarget, taxType=taxType):
                        continue

                    #
                    analD = eD['anal_instances'] if 'anal_instances' in eD else {}

                    for asymId, aD in analD.items():

                        gapCount = len(aD['gapD'])
                        tL = list(aD['gapD'].values())
                        tL = [t if t > 0 else 0 for t in tL]
                        gapL = tL if len(tL) > 0 else [0]
                        gapCountList.append(gapCount)
                        gapLengthList.extend(gapL)
            #
            logger.info("gaps %d gap lengths %d" % (len(gapCountList), len(gapLengthList)))
            #
            fnBase = 'gap-count-fig-'
            fn = self.__getPlotFileName(fnBase, resLimit=resLimit, polymerComposition=polymerComposition, taxIdTarget=taxIdTarget, taxType=taxType)
            fp = os.path.join(self.__workPath, fn)
            legendBase = "Gap counts for (%d) protein sequences" % len(gapCountList)
            legend = self.__getLegend(legendBase, numEntries=numEntries, resLimit=resLimit, polymerComposition=polymerComposition, taxIdTarget=taxIdTarget, taxType=taxType)

            cu = DisorderChartUtils()
            # cu.doIntegerBarChart(gapCountList, plotPath=self.__plotGapCount, yPlotScale=None, yPlotMax=300000)
            cu.doIntegerBarChart(gapCountList, plotPath=fp, yPlotScale='log',
                                 yPlotMax=6, xPlotMax=30, xPlotLabel="Gap Count", yPlotLabel="Protein Instances (log)",
                                 plotTitle="Protein Instance Gap Count")
            #
            self.__writeLegend(fp, legend)
            cu.doIntegerBarChart(gapLengthList, plotPath=self.__plotGapLength, yPlotScale='log',
                                 yPlotMax=6, xPlotMax=150, xPlotLabel="Gap width (residues)", yPlotLabel="Gap Instances (log)",
                                 plotTitle="Protein Instance Gap Widths")

            #
            legendBase = "Gap width statistics for (%d) protein sequences" % len(gapCountList)
            legend = self.__getLegend(legendBase, numEntries=numEntries, resLimit=resLimit, polymerComposition=polymerComposition, taxIdTarget=taxIdTarget, taxType=taxType)
            fnBase = 'gap-width-fig-'
            fn = self.__getPlotFileName(fnBase, resLimit=resLimit, polymerComposition=polymerComposition, taxIdTarget=taxIdTarget, taxType=taxType)
            fp = os.path.join(self.__workPath, fn)
            self.__writeLegend(fp, legend)
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testViewOccData(self):
        """ Test case - protein sequence instance occupancy segment counts and segment widths all cases
        """
        try:
            entryD = self.__mU.doImport(self.__instanceSavePath, format="pickle")
            segmentCountList = []
            segmentLengthList = []
            for entryId in entryD:
                for entityId, eD in entryD[entryId]['selected_polymer_entities'].items():

                    analD = eD['anal_instances'] if 'anal_instances' in eD else {}

                    for asymId, aD in analD.items():
                        segmentCount = len(aD['owabRegiond'])
                        segmentLengths = [d['length'] for sId, d in aD['owabRegiond'].items()]

                        segmentCountList.append(segmentCount)
                        segmentLengthList.extend(segmentLengths)
            #
            logger.info("gaps %d gap lengths %d" % (len(segmentCountList), len(segmentLengthList)))
            #
            cu = DisorderChartUtils()
            cu.doIntegerBarChart(segmentCountList, plotPath=self.__plotOwabSegmentCount, yPlotScale='log',
                                 yPlotMax=6, xPlotMax=100, xPlotLabel="Segment Count", yPlotLabel="Protein Instances (log)",
                                 plotTitle="Segment counts (OWAB > 2 * mean OWAB)")
            self.__writeLegend(
                self.__plotOwabSegmentCount,
                "Segment counts for all (%d) protein sequences (OWAB > 2 * mean OWAB and X-ray resolution limit < 3.5 Angstoms) " %
                len(segmentCountList))
            cu.doIntegerBarChart(segmentLengthList, plotPath=self.__plotOwabSegmentLength, yPlotScale='log',
                                 yPlotMax=6, xPlotMax=100, xPlotLabel="Segment width (residues)", yPlotLabel="Segmnt Instances (log)",
                                 plotTitle="Segment widths (OWAB > 2 * mean OWAB)")
            self.__writeLegend(
                self.__plotOwabSegmentLength,
                "Segment widths for all (%d) protein sequences (OWAB > 2 * mean OWAB and X-ray resolution limit < 3.5 Angstoms) " %
                len(segmentLengthList))
        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()

    def testViewCoverageData(self):
        """ Test case - protein sequence reference db and sample sequence coverage all cases

        """
        try:
            entryD = self.__mU.doImport(self.__instanceSavePath, format="pickle")
            covRefDbList = []
            covSampleList = []
            for entryId in entryD:
                for entityId, eD in entryD[entryId]['selected_polymer_entities'].items():

                    analD = eD['anal_instances'] if 'anal_instances' in eD else {}

                    for asymId, aD in analD.items():
                        covRefDb = aD['coverage_inst_refdb']
                        covSample = aD['coverage_inst_entity']
                        if covRefDb is not None:
                            covRefDb = 0.0 if covRefDb < 0.0 else covRefDb
                            covRefDb = 1.0 if covRefDb > 1.0 else covRefDb
                            covRefDbList.append(covRefDb)
                        if covSample is not None:
                            covSample = 0.0 if covSample < 0.0 else covSample
                            covSample = 1.0 if covSample > 1.0 else covSample
                            covSampleList.append(covSample)
            #
            logger.info("covRefDbList %d covSampleList %d" % (len(covRefDbList), len(covSampleList)))
            #
            cu = DisorderChartUtils()
            cu.doHistogramChart(covRefDbList, plotPath=self.__plotCoverageRefDb,
                                yPlotScale='log',
                                yPlotMax=100000, yPlotMin=1000,
                                # yPlotMax=100000,
                                xPlotLabel="Coverage Fraction",
                                yPlotLabel="Protein Instances",
                                plotTitle="Reference Sequence Coverage")
            self.__writeLegend(
                self.__plotCoverageRefDb,
                "UniProt reference sequence coverage for all (%d) protein sequences (X-ray resolution limit < 3.5 Angstoms) " %
                len(covRefDbList))
            cu.doHistogramChart(covSampleList, plotPath=self.__plotCoverageSample1,
                                yPlotScale='log',
                                yPlotMax=100000, yPlotMin=1000,
                                # yPlotMax=100000,
                                xPlotLabel="Coverage Fraction",
                                yPlotLabel="Protein Instances",
                                plotTitle="Sample Sequence Coverage")
            self.__writeLegend(self.__plotCoverageSample1, "Sample sequence coverage for all (%d) protein sequences (X-ray resolution limit < 3.5 Angstoms) " % len(covSampleList))
            #
            cu.doHistogramChart(covSampleList, plotPath=self.__plotCoverageSample2,
                                yPlotScale='log',
                                yPlotMax=100000, yPlotMin=1000,
                                xPlotMin=0.8, xPlotMax=1.001, xPlotIncr=0.1,
                                # yPlotMax=100000,
                                xPlotLabel="Coverage Fraction",
                                yPlotLabel="Protein Instances",
                                plotTitle="Sample Sequence Coverage")
            self.__writeLegend(self.__plotCoverageSample1, "Sample sequence coverage for all (%d) protein sequences (X-ray resolution limit < 3.5 Angstoms) " % len(covSampleList))

        except Exception as e:
            logger.exception("Failing with %s" % str(e))
            self.fail()


def filterSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DisorderChartUtilsFilterTests("testPlotFilteredGapData"))
    return suiteSelect


if __name__ == '__main__':
    if (True):
        mySuite = filterSuite()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
