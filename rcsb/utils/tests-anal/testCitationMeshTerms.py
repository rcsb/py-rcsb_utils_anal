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

from collections import defaultdict
import copy
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


class CitationMeshTermsTests(unittest.TestCase):
    def __init__(self, methodName="runTest"):
        super(CitationMeshTermsTests, self).__init__(methodName)
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

    def testExportMesh(self):
        """Export MESH terms for PDB primary citations"""
        ommitL = [
            "Models, Molecular",
            "Protein Conformation",
            "Amino Acid Sequence",
            "Binding Sites",
            "Molecular Sequence Data",
            "Protein Binding",
            "Animals",
            "Protein Structure, Tertiary",
            "Protein Structure, Secondary",
            "Structure-Activity Relationship",
            "Catalytic Domain",
            "Sequence Homology, Amino Acid",
            "Substrate Specificity",
            "Mutation",
            "Sequence Alignment",
            "Recombinant Proteins",
            "Kinetics",
            "Crystallization",
            "Ligands",
            "Molecular Structure",
            "Protein Structure, Quaternary",
            "Protein Folding",
            "Protein Multimerization",
            "Nucleic Acid Conformation",
            "Enzyme Inhibitors",
            "Base Sequence",
            "Amino Acid Substitution",
            "Protein Domains",
            "Amino Acid Motifs",
            "Hydrogen-Ion Concentration",
            "Cell Line",
            "Models, Chemical",
            "Models, Biological",
            "Macromolecular Substances",
        ]
        mU = MarshalUtil()
        enD = mU.doImport(os.path.join(self.__workPath, "entry-citation-details-example.json"), fmt="json")
        pmD = mU.doImport(os.path.join(self.__workPath, "entry-citation-mesh-example.json"), fmt="json")
        #
        logger.info("pmD (%d)", len(pmD["entry"]))
        mCount = defaultdict(int)
        for pmId, pD in pmD["entry"].items():
            if "rcsb_pubmed_mesh_descriptors" in pD:
                meshL = pD["rcsb_pubmed_mesh_descriptors"]
                for mesh in meshL:
                    mCount[mesh] += 1
        mCountS = sorted(mCount.items(), key=lambda kv: kv[1], reverse=True)
        tL = []
        logger.info("mCountS (%d)", len(mCount))
        for m, num in mCountS:
            if num > 500:
                tL.append((m, num))
        #
        logger.info("%r", tL)
        #
        pubMedD = pmD["entry"]
        #
        rowList = []
        #
        for entryId, obj in enD["entry"].items():
            dD = {}
            doi = ""
            pubMed = None
            if "rcsb_accession_info" in obj:
                relYr = int(obj["rcsb_accession_info"]["initial_release_date"][:4])
            if "rcsb_primary_citation" in obj and "pdbx_database_id_PubMed" in obj["rcsb_primary_citation"]:
                pubMed = str(obj["rcsb_primary_citation"]["pdbx_database_id_PubMed"])
            if "rcsb_primary_citation" in obj and "pdbx_database_id_DOI" in obj["rcsb_primary_citation"]:
                doi = obj["rcsb_primary_citation"]["pdbx_database_id_DOI"]
            #
            if pubMed and pubMed in pubMedD and "rcsb_pubmed_mesh_descriptors" in pubMedD[pubMed]:
                meshL = [m for m in pubMedD[pubMed]["rcsb_pubmed_mesh_descriptors"] if m not in ommitL]
                if "rcsb_pubmed_doi" in pubMedD[pubMed]:
                    doi = pubMedD[pubMed]["rcsb_pubmed_doi"]
                dD = {"entry": entryId, "release_year": relYr, "pubmed": pubMed, "doi": doi, "mesh_terms": ";".join(meshL)}

                rowList.append(dD)
        #
        logger.info("rowList length %d", len(rowList))
        ok = mU.doExport(os.path.join(self.__workPath, "pdb-citation-mesh-terms.json"), rowList, fmt="json", indent=3)
        ok = mU.doExport(os.path.join(self.__workPath, "pdb-citation-mesh-terms.csv"), rowList, fmt="csv")

    def testExtractCitationEntryDetails(self):
        """Test case - extract selected entry details

        "pdbx_database_id_PubMed" : 24305054,
        "pdbx_database_id_DOI" : "10.1038/nature12725",
        """
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="dw",
                collectionName="core_entry",
                cacheFilePath=os.path.join(self.__workPath, "entry-citation-details-example.json"),
                useCache=False,
                keyAttribute="entry",
                uniqueAttributes=["rcsb_id"],
                cacheKwargs=self.__testEntryCacheKwargs,
                objectLimit=self.__objectLimitTest,
                selectionList=[
                    "rcsb_id",
                    "rcsb_primary_citation.pdbx_database_id_PubMed",
                    "rcsb_primary_citation.pdbx_database_id_DOI",
                    "rcsb_accession_info.initial_release_date",
                ],
            )
            eCount = obEx.getCount()
            logger.info("Entry count is %d", eCount)

            objD = obEx.getObjects()
            for _, obj in objD.items():
                rcsbId = obj["rcsb_id"]
                if "rcsb_primary_citation" in obj and "pdbx_database_id_PubMed" in obj["rcsb_primary_citation"]:
                    logger.debug("%s PubMed %s", rcsbId, obj["rcsb_primary_citation"]["pdbx_database_id_PubMed"])

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExtractMeshTerms(self):
        """Test case - extract mesh terms"""
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="dw",
                collectionName="core_pubmed",
                cacheFilePath=os.path.join(self.__workPath, "entry-citation-mesh-example.json"),
                useCache=False,
                keyAttribute="entry",
                uniqueAttributes=["rcsb_id"],
                cacheKwargs=self.__testEntryCacheKwargs,
                objectLimit=self.__objectLimitTest,
                selectionList=["rcsb_id", "rcsb_pubmed_mesh_descriptors", "rcsb_pubmed_doi"],
            )
            eCount = obEx.getCount()
            logger.info("PubMed count is %d", eCount)

            objD = obEx.getObjects()
            for _, obj in objD.items():
                rcsbId = obj["rcsb_id"]
                if "rcsb_pubmed_mesh_descriptors" in obj:
                    logger.debug("%s Mesh %s", rcsbId, obj["rcsb_pubmed_mesh_descriptors"])

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def citationSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CitationMeshTermsTests("testExtractCitationEntryDetails"))
    suiteSelect.addTest(CitationMeshTermsTests("testExtractMeshTerms"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = citationSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
