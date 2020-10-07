##
# File:    testExtractHumanProteinData.py
# Author:  J. Westbrook
# Date:    7-Oct-2020
#
# Updates:
#
##
"""
Larger tests/examples of extractor functions selecting values from collections.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"


import logging
import os

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


class ExampleExtractorTests(unittest.TestCase):
    def __init__(self, methodName="runTest"):
        super(ExampleExtractorTests, self).__init__(methodName)
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
        self.__workPath = os.path.join(TOPDIR, "CACHE", "exdb")
        #
        self.__testEntryCacheKwargs = {"fmt": "json", "indent": 3}
        self.__objectLimitTest = None
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testExtractEntryDetails(self):
        """Test case - extract selected entry details"""
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="pdbx_core",
                collectionName="pdbx_core_entry",
                cacheFilePath=os.path.join(self.__workPath, "entry-data-details-example.json"),
                useCache=False,
                keyAttribute="entry",
                uniqueAttributes=["rcsb_id"],
                cacheKwargs=self.__testEntryCacheKwargs,
                objectLimit=self.__objectLimitTest,
                selectionList=["rcsb_id", "rcsb_accession_info", "struct", "exptl"],
            )
            eCount = obEx.getCount()
            logger.info("Entry count is %d", eCount)

            objD = obEx.getObjects()
            for _, obj in objD.items():
                rcsbId = obj["rcsb_id"]
                logger.debug("%s rcsb_accession_info %s", rcsbId, obj["rcsb_accession_info"]["initial_release_date"])

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExtractEntityTaxonomyAndRefDetails(self):
        """Test case - extract unique entity source and host taxonomies"""
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="pdbx_core",
                collectionName="pdbx_core_polymer_entity",
                cacheFilePath=os.path.join(self.__workPath, "entity-taxonomy-ref-example.json"),
                useCache=False,
                keyAttribute="entity",
                uniqueAttributes=["rcsb_id"],
                cacheKwargs=self.__testEntryCacheKwargs,
                objectLimit=None,
                selectionQuery={"entity_poly.rcsb_entity_polymer_type": "Protein", "rcsb_entity_source_organism.ncbi_taxonomy_id": 9606},
                selectionList=[
                    "rcsb_id",
                    "rcsb_entity_source_organism.ncbi_taxonomy_id",
                    "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers",
                    "rcsb_polymer_entity",
                ],
            )
            eCount = obEx.getCount()
            logger.info("Polymer entity count is %d", eCount)
            taxIdS = set()
            objD = obEx.getObjects()
            for _, eD in objD.items():
                try:
                    for tD in eD["rcsb_entity_source_organism"]:
                        taxIdS.add(tD["ncbi_taxonomy_id"])
                except Exception:
                    pass
            logger.info("Unique taxons %d", len(taxIdS))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def __getEntryDetails(self):
        mU = MarshalUtil(workPath=self.__workPath)
        entryFilePath = os.path.join(self.__workPath, "entry-data-details-example.json")
        tD = mU.doImport(entryFilePath, fmt="json")
        return tD["entry"]

    def testFindFirst(self):
        """Find first reference sequence in PDB

        Returns:
            [type]: [description]

        Example entity object:

            "1ERU_1": {
            "rcsb_entity_source_organism": [
                {
                "ncbi_taxonomy_id": 9606
                }
            ],
            "rcsb_polymer_entity": {
                "details": "ACTIVE SITE CYSTEINES 32 AND 35 IN THE OXIDIZED FORM",
                "formula_weight": 11.75,
                "src_method": "man",
                "rcsb_multiple_source_flag": "N",
                "rcsb_source_part_count": 1,
                "rcsb_source_taxonomy_count": 1,
                "pdbx_description": "THIOREDOXIN",
                "pdbx_number_of_molecules": 1,
                "rcsb_macromolecular_names_combined": [
                {
                    "name": "THIOREDOXIN",
                    "provenance_code": "ECO:0000304",
                    "provenance_source": "PDB Preferred Name"
                }
                ]
            },
            "rcsb_polymer_entity_container_identifiers": {
                "reference_sequence_identifiers": [
                {
                    "database_name": "UniProt",
                    "database_accession": "P10599",
                    "provenance_source": "SIFTS"
                }
                ]
            },
            "rcsb_id": "1ERU_1"
            },
        """
        tU = TimeUtil()
        entryD = self.__getEntryDetails()

        # tD = tU.getDateTimeObj(tS)

        seqRefD = {}
        multiPartCount = 0
        missCount = 0
        mU = MarshalUtil(workPath=self.__workPath)
        entityFilePath = os.path.join(self.__workPath, "entity-taxonomy-ref-example.json")
        entityD = mU.doImport(entityFilePath, fmt="json")
        logger.info("entityD.keys() %s", list(entityD.keys()))
        for entityId, eObj in entityD["entity"].items():
            # logger.info("eObj %r", eObj)
            # entityId = eObj["rcsb_id"]
            entryId = entityId[:4]
            releaseTs = entryD[entryId]["rcsb_accession_info"]["initial_release_date"]
            releaseDt = tU.getDateTimeObj(releaseTs)
            partCount = eObj["rcsb_polymer_entity"]["rcsb_source_part_count"]
            if partCount > 1:
                multiPartCount += 1
                continue
            ok = False
            try:
                rD = eObj["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"][0]
                if rD["database_name"] == "UniProt":
                    uniProtId = rD["database_accession"]
                    seqRefD.setdefault(uniProtId, []).append((entityId, releaseDt))
                    ok = True
            except Exception:
                pass
            if not ok:
                missCount += 1
                logger.info("Missing Uniprot reference for %r", entityId)

        logger.info("Length seqRefD %d", len(seqRefD))
        logger.info("mulitPartCount %d", multiPartCount)
        logger.info("missCount %d", missCount)


def objectExtractorExampleSuite():

    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ExampleExtractorTests("testExtractEntryDetails"))
    suiteSelect.addTest(ExampleExtractorTests("testExtractEntityTaxonomyAndRefDetails"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = objectExtractorExampleSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
