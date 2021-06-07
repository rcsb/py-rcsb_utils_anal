##
# File:    testChemRefStats.py
# Author:  J. Westbrook
# Date:    14-Apr-2021
#
# Updates:
#
##
"""
Chemical reference data statistics extracted from the ExDB collections.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import copy
import logging
import os
import re
import time
import unittest

from rcsb.exdb.utils.ObjectExtractor import ObjectExtractor
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.TimeUtil import TimeUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class ChemRefStatsTests(unittest.TestCase):
    def __init__(self, methodName="runTest"):
        super(ChemRefStatsTests, self).__init__(methodName)
        self.__verbose = True

    def setUp(self):
        #
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        configPath = os.path.join(TOPDIR, "rcsb", "mock-data", "config", "dbload-setup-example.yml")
        #
        configName = "site_info_remote_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        #
        self.__resourceName = "MONGO_DB"
        # self.__workPath = os.path.join(TOPDIR, "CACHE", "exdb")
        self.__workPath = os.path.join(HERE, "test-output")
        #
        self.__testEntryCacheKwargs = {"fmt": "json", "indent": 3}
        self.__objectLimitTest = None
        #
        self.__totalMolecules = 33603 + 696
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testExtractCcTotals(self):
        """Test case - extract total released components"""
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="bird_chem_comp_core",
                collectionName="bird_chem_comp_core",
                cacheFilePath=os.path.join(self.__workPath, "cc-all-cache.json"),
                useCache=False,
                cacheKwargs=self.__testEntryCacheKwargs,
                keyAttribute="chem_comp",
                uniqueAttributes=["rcsb_id"],
                selectionQuery={
                    "rcsb_id": {"$not": re.compile("^PRD_")},
                    "rcsb_chem_comp_info.release_status": {"$eq": "REL"},
                },
                selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers"],
            )
            eCount = obEx.getCount()
            logger.info("CC total count (%d)", eCount)
            self.assertGreaterEqual(eCount, 3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExtractBirdTotals(self):
        """Test case - extract total released components"""
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="bird_chem_comp_core",
                collectionName="bird_chem_comp_core",
                cacheFilePath=os.path.join(self.__workPath, "bird-all-cache.json"),
                useCache=False,
                cacheKwargs=self.__testEntryCacheKwargs,
                keyAttribute="chem_comp",
                uniqueAttributes=["rcsb_id"],
                selectionQuery={
                    "rcsb_id": {"$regex": re.compile("^PRD_")},
                    "rcsb_chem_comp_info.release_status": {"$eq": "REL"},
                },
                selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers"],
            )
            eCount = obEx.getCount()
            logger.info("Bird total count (%d)", eCount)
            self.assertGreaterEqual(eCount, 3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExtractDrugbankMapping(self):
        """Test case - extract Drugbank mapping"""
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="bird_chem_comp_core",
                collectionName="bird_chem_comp_core",
                cacheFilePath=os.path.join(self.__workPath, "drugbank-mapping-cache.json"),
                useCache=False,
                cacheKwargs=self.__testEntryCacheKwargs,
                keyAttribute="chem_comp",
                uniqueAttributes=["rcsb_id"],
                selectionQuery={
                    "rcsb_chem_comp_container_identifiers.drugbank_id": {"$exists": True},
                },
                selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers", "rcsb_chem_comp_related"],
            )
            eCount = obEx.getCount()
            logger.info("DrugBank mapping count (%d)", eCount)
            self.assertGreaterEqual(eCount, 3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExtractMapping(self):
        """Test case - extract related chemref mappings"""
        try:
            refD = {}
            for resourceName in ["PubChem", "ChEMBL", "ChEBI", "DrugBank", "Pharos"]:
                obEx = ObjectExtractor(
                    self.__cfgOb,
                    databaseName="bird_chem_comp_core",
                    collectionName="bird_chem_comp_core",
                    cacheFilePath=os.path.join(self.__workPath, "chemref-%s-mapping-cache.json" % resourceName),
                    useCache=False,
                    cacheKwargs=self.__testEntryCacheKwargs,
                    keyAttribute="chem_comp",
                    uniqueAttributes=["rcsb_id"],
                    selectionQuery={"rcsb_chem_comp_related.resource_name": {"$eq": resourceName}},
                    selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers", "rcsb_chem_comp_related"],
                )
                eCount = obEx.getCount()
                logger.info("%s mapping count is (%d)", resourceName, eCount)
                self.assertGreaterEqual(eCount, 3)
                refD[resourceName] = eCount
            logger.info("Resource mapping counts %r", refD)
            refD = {k: round(float(v) / float(self.__totalMolecules), 3) for k, v in refD.items()}
            logger.info("Resource mapping coverage pc %r", refD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testChangeCcByYear(self):
        """Test case - new and modified chemical components by year"""
        try:

            tU = TimeUtil()
            yrTotalD = {}
            yearList = [2016, 2017, 2018, 2019, 2020, 2021]
            for year in yearList:
                tsEnd = "%s-12-31" % str(year)
                tdEnd = tU.getDateTimeObj(tsEnd)
                obEx = ObjectExtractor(
                    self.__cfgOb,
                    databaseName="bird_chem_comp_core",
                    collectionName="bird_chem_comp_core",
                    cacheFilePath=os.path.join(self.__workPath, "chemref-cc-total-%d-cache.json" % year),
                    useCache=False,
                    cacheKwargs=self.__testEntryCacheKwargs,
                    keyAttribute="chem_comp",
                    uniqueAttributes=["rcsb_id"],
                    selectionQuery={
                        "rcsb_id": {"$not": re.compile("^PRD_")},
                        "rcsb_chem_comp_info.release_status": {"$eq": "REL"},
                        "rcsb_chem_comp_info.initial_release_date": {"$lte": tdEnd},
                    },
                    selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers", "chem_comp"],
                )
                eCount = obEx.getCount()
                logger.info("Total components in %d (%d)", year, eCount)
                yrTotalD[year] = eCount
            logger.info("Total components %r", yrTotalD)
            yrTotalPcD = {}
            for year in yearList[1:]:
                yrTotalPcD[year] = round(float(yrTotalD[year] - yrTotalD[year - 1]) / (yrTotalD[year]), 3)
            logger.info("Total pc increase by year %r", yrTotalPcD)
            #
            #
            for chngType in ["initial_release_date", "revision_date"]:
                yrD = {}
                yearList = [2017, 2018, 2019, 2020, 2021]
                for year in yearList:
                    tsStart = "%s-01-01" % str(year)
                    tdStart = tU.getDateTimeObj(tsStart)
                    tsEnd = "%s-12-31" % str(year)
                    tdEnd = tU.getDateTimeObj(tsEnd)
                    obEx = ObjectExtractor(
                        self.__cfgOb,
                        databaseName="bird_chem_comp_core",
                        collectionName="bird_chem_comp_core",
                        cacheFilePath=os.path.join(self.__workPath, "chemref-cc-%d-cache.json" % year),
                        useCache=False,
                        cacheKwargs=self.__testEntryCacheKwargs,
                        keyAttribute="chem_comp",
                        uniqueAttributes=["rcsb_id"],
                        selectionQuery={
                            "rcsb_id": {"$not": re.compile("^PRD_")},
                            "rcsb_chem_comp_info.release_status": {"$eq": "REL"},
                            "rcsb_chem_comp_info.%s" % chngType: {"$gte": tdStart, "$lte": tdEnd},
                        },
                        selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers", "chem_comp"],
                    )
                    eCount = obEx.getCount()
                    logger.info("%s components in %d (%d)", chngType, year, eCount)
                    yrD[year] = eCount
                    self.assertGreaterEqual(eCount, 3)

                logger.info("%s components by year %r", chngType, yrD)
                yrPcD = {}
                for year in yearList:
                    yrPcD[year] = round(float(yrD[year]) / float(yrTotalD[year]), 3)
                logger.info("%s pc by year %r", chngType, yrPcD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testChangeBirdByYear(self):
        """Test case - new and modified chemical components by year"""
        try:
            tU = TimeUtil()
            yrTotalD = {}
            yearList = [2016, 2017, 2018, 2019, 2020, 2021]
            for year in yearList:
                tsEnd = "%s-12-31" % str(year)
                tdEnd = tU.getDateTimeObj(tsEnd)
                obEx = ObjectExtractor(
                    self.__cfgOb,
                    databaseName="bird_chem_comp_core",
                    collectionName="bird_chem_comp_core",
                    cacheFilePath=os.path.join(self.__workPath, "chemref-cc-total-%d-cache.json" % year),
                    useCache=False,
                    cacheKwargs=self.__testEntryCacheKwargs,
                    keyAttribute="chem_comp",
                    uniqueAttributes=["rcsb_id"],
                    selectionQuery={
                        # "rcsb_id": {"$regex": re.compile("^PRD_")},
                        "rcsb_chem_comp_container_identifiers.prd_id": {"$exists": True},
                        "rcsb_chem_comp_info.initial_release_date": {"$exists": True},
                        # "rcsb_chem_comp_info.release_status": {"$eq": "REL"},
                        # "rcsb_chem_comp_info.initial_release_date": {"$lte": tdEnd},
                    },
                    selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers", "chem_comp"],
                )
                eCount = obEx.getCount()
                logger.info("Total birds in %d (%d)", year, eCount)
                yrTotalD[year] = eCount
            logger.info("Total birds %r", yrTotalD)
            yrTotalPcD = {}
            for year in yearList[1:]:
                yrTotalPcD[year] = round(float(yrTotalD[year] - yrTotalD[year - 1]) / (yrTotalD[year]), 3)
            logger.info("Total pc bird increase by year %r", yrTotalPcD)
            #
            for chngType in ["initial_release_date", "revision_date"]:
                yrD = {}
                yearList = [2017, 2018, 2019, 2020, 2021]
                for year in yearList:
                    tsStart = "%s-01-01" % str(year)
                    tdStart = tU.getDateTimeObj(tsStart)
                    tsEnd = "%s-12-31" % str(year)
                    tdEnd = tU.getDateTimeObj(tsEnd)
                    obEx = ObjectExtractor(
                        self.__cfgOb,
                        databaseName="bird_chem_comp_core",
                        collectionName="bird_chem_comp_core",
                        cacheFilePath=os.path.join(self.__workPath, "chemref-cc-%d-cache.json" % year),
                        useCache=False,
                        cacheKwargs=self.__testEntryCacheKwargs,
                        keyAttribute="chem_comp",
                        uniqueAttributes=["rcsb_id"],
                        selectionQuery={
                            "rcsb_chem_comp_container_identifiers.prd_id": {"$exists": True},
                            "rcsb_chem_comp_info.release_status": {"$eq": "REL"},
                            "rcsb_chem_comp_info.%s" % chngType: {"$gte": tdStart, "$lte": tdEnd},
                        },
                        selectionList=["rcsb_id", "rcsb_chem_comp_container_identifiers", "chem_comp"],
                    )
                    eCount = obEx.getCount()
                    logger.info("%s components in %d (%d)", chngType, year, eCount)
                    yrD[year] = eCount
                    self.assertGreaterEqual(eCount, 3)
                logger.info("%s components by year %r", chngType, yrD)
                yrPcD = {}
                for year in yearList:
                    yrPcD[year] = round(float(yrD[year]) / float(yrTotalD[year]), 3)
                logger.info("%s pc by year %r", chngType, yrPcD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def chemRefStatsSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemRefStatsTests("testExtractDrugbankMapping"))
    suiteSelect.addTest(ChemRefStatsTests("testExtractNewMolecules"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = chemRefStatsSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
