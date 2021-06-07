##
# File:    testCitationMeshTerms.py
# Author:  J. Westbrook
# Date:    18-Nov-2020
#
# Updates:
#
##
"""
Fetch mesh terms for primary citations.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest

from rcsb.db.mongo.Connection import Connection
from rcsb.db.mongo.MongoDbUtil import MongoDbUtil
from rcsb.utils.config.ConfigUtil import ConfigUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class CountEnumsTests(unittest.TestCase):
    def __init__(self, methodName="runTest"):
        super(CountEnumsTests, self).__init__(methodName)
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
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testCountEnumsTerms(self):
        """Test case - extract mesh terms"""
        try:
            pass
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testCountPolymerEnums(self):
        """Test case - count selected enums """
        try:
            dbName = "pdbx_core"
            collectionName = "pdbx_core_entry"
            with Connection(cfgOb=self.__cfgOb, resourceName=self.__resourceName) as client:
                mg = MongoDbUtil(client)
                ii = mg.count(dbName, collectionName)
                logger.info("%s %s collection length %d", dbName, collectionName, ii)
                #
                atNameL = ["rcsb_entry_info.polymer_composition", "rcsb_entry_info.selected_polymer_entity_types", "rcsb_entry_info.experimental_method"]
                for atName in atNameL:
                    vL0 = mg.distinct(dbName, collectionName, atName)
                    logger.info("enum values %r", vL0)
                    for v in vL0:
                        num = mg.count(dbName, collectionName, countFilter={atName: v})
                        logger.info("%s value %s (%d)", atName, v, num)
                #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testCountFileTypeEnums(self):
        """Test case - count selected enums """
        try:
            dbName = "repository_holdings"
            collectionName = "repository_holdings_current_entry"
            with Connection(cfgOb=self.__cfgOb, resourceName=self.__resourceName) as client:
                mg = MongoDbUtil(client)
                ii = mg.count(dbName, collectionName)
                logger.info("%s %s collection length %d", dbName, collectionName, ii)
                #
                atNameL = ["rcsb_repository_holdings_current.repository_content_types"]
                for atName in atNameL:
                    vL0 = mg.distinct(dbName, collectionName, atName)
                    logger.info("enum values %r", vL0)
                    for v in vL0:
                        num = mg.count(dbName, collectionName, countFilter={atName: {"$in": [v]}})
                        logger.info("%s value %s (%d)", atName, v, num)
                #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def enumSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CountEnumsTests("testCountEnums"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = enumSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
