# File:    NeighborUtilsTests.py
# Author:  J. Westbrook
# Date:    12-Feb-2021
# Version: 0.001
#
# Update:
#
##
"""
Tests for finding ligand macromolecular neighbors.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest
from collections import namedtuple
from operator import itemgetter
import numpy as np
from scipy import spatial

from rcsb.utils.dictionary.DictMethodCommonUtils import DictMethodCommonUtils
from rcsb.utils.repository.RepositoryProvider import RepositoryProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


LigandTargetFields = (
    "ligandCompId",
    "ligandAtomName",
    "connectType",
    "partnerEntityType",
    "partnerEntityId",
    "partnerCompId",
    "partnerAsymId",
    "partnerSeqId",
    "partnerAtomName",
    "distance",
)
LigandTargetInstance = namedtuple("LigandTargetInstance", LigandTargetFields, defaults=(None,) * len(LigandTargetFields))

ReferenceFields = ("entityId", "entityType", "asymId", "compId", "seqId", "atomName")
ReferenceInstance = namedtuple("ReferenceInstance", ReferenceFields, defaults=(None,) * len(ReferenceFields))


class NeighborUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__numProc = 2
        self.__fileLimit = 200
        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_configuration"
        self.__configName = configName
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=mockTopPath)
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        self.__rpP = RepositoryProvider(cfgOb=self.__cfgOb, numProc=self.__numProc, fileLimit=self.__fileLimit, cachePath=self.__cachePath)
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def __getContentLocators(self, contentType="pdbx_core", mergeContent=None):
        """Read and process test fixture data files from the input content type."""
        try:
            mergeContent = mergeContent if mergeContent else []
            locatorObjList = self.__rpP.getLocatorObjList(contentType=contentType, mergeContentTypes=mergeContent)
            # containerList = self.__rpP.getContainerList(locatorObjList)
            return locatorObjList
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()
        return []

    def testGetNeighbors(self):
        """Test get neighbors."""
        commonU = DictMethodCommonUtils()
        ltiD = {}
        distLimit = 5.0
        locatorObjList = self.__getContentLocators()
        for locatorObj in locatorObjList:
            dataContainerList = self.__rpP.getContainerList([locatorObj])
            for dataContainer in dataContainerList:
                entryId = dataContainer.getName()
                # instD = self.__getNonpolymerInstanceNeighbors(dataContainer, distLimit=distLimit)
                instD = commonU.getNonpolymerInstanceNeighbors(dataContainer, distLimit=distLimit)
                ltiD[entryId] = instD
        for entryId, nptD in ltiD.items():
            for asymId, nL in nptD.items():
                logger.info("%s (%s) %r", entryId, asymId, nL)
        #
        logger.info("Returning %d", len(ltiD))

    def __getNonpolymerInstanceNeighbors(self, dataContainer, distLimit=5.0):
        """Get bound and unbound neighbors for each non-polymer instance.

        Args:
            dataContainer (obj): DataContainer object
            distLimit (float, optional): neighbor distance limit (Angstroms). Defaults to 5.0.

        Returns:
              (dict): {asymId: [LigandTargetInstance()]}
        """
        try:
            startTime = time.time()
            ligandTargetInstanceD = {}
            commonU = DictMethodCommonUtils()
            instanceTypeD = commonU.getInstanceTypes(dataContainer)
            if "non-polymer" not in instanceTypeD.values():
                return ligandTargetInstanceD
            #
            entryId = dataContainer.getName()
            logger.info("Starting entry %s", entryId)
            nonPolymerBoundD = commonU.getBoundNonpolymersByInstance(dataContainer)
            #
            instanceTypeD = commonU.getInstanceTypes(dataContainer)
            instancePolymerTypeD = commonU.getInstancePolymerTypes(dataContainer)
            instanceEntityD = commonU.getInstanceEntityMap(dataContainer)
            # -----
            targetXyzL = []
            targetRefL = []
            ligandXyzD = {}
            ligandRefD = {}
            # partition the cooordinates between ligands and candidate targets
            aObj = dataContainer.getObj("atom_site")
            for ii in range(aObj.getRowCount()):
                selectType = None
                asymId = aObj.getValue("label_asym_id", ii)
                instanceType = instanceTypeD[asymId]
                polymerType = instancePolymerTypeD[asymId] if asymId in instancePolymerTypeD else None
                if (instanceType == "polymer" and polymerType in ["Protein", "DNA", "RNA", "NA-hybrid"]) or instanceType == "branched":
                    selectType = "target"
                elif instanceType == "non-polymer":
                    selectType = "ligand"
                if selectType not in ["target", "ligand"]:
                    continue
                #
                atomName = aObj.getValue("label_atom_id", ii)
                seqId = aObj.getValue("label_seq_id", ii)
                compId = aObj.getValue("label_comp_id", ii)
                xC = aObj.getValue("Cartn_x", ii)
                yC = aObj.getValue("Cartn_y", ii)
                zC = aObj.getValue("Cartn_z", ii)
                entityId = instanceEntityD[asymId]
                #
                if selectType == "target":
                    targetXyzL.append((float(xC), float(yC), float(zC)))
                    targetRefL.append(ReferenceInstance(entityId, instanceType, asymId, compId, int(seqId) if seqId not in [".", "?"] else None, atomName))
                elif selectType == "ligand":
                    ligandXyzD.setdefault(asymId, []).append((float(xC), float(yC), float(zC)))
                    ligandRefD.setdefault(asymId, []).append(ReferenceInstance(entityId, instanceType, asymId, compId, None, atomName))
            #
            # ------
            logger.debug("%s targetXyzL (%d) targetRef (%d) ligandXyzD (%d) ", entryId, len(targetXyzL), len(targetXyzL), len(ligandXyzD))
            if not targetXyzL:
                return ligandTargetInstanceD
            tArr = np.array(targetXyzL, order="F")
            logger.debug("targetXyzL[0] %r tArr.shape %r tArr[0] %r", targetXyzL[0], tArr.shape, tArr[0])
            tree = spatial.cKDTree(tArr)
            #
            for asymId, ligXyzL in ligandXyzD.items():
                #
                if asymId in nonPolymerBoundD:
                    # Process bound ligands
                    for tup in nonPolymerBoundD[asymId]:
                        if tup.partnerEntityType not in ["non-polymer", "water"]:
                            ligandTargetInstanceD.setdefault(asymId, []).append(
                                LigandTargetInstance(
                                    tup.targetCompId,
                                    tup.targetAtomId,
                                    tup.connectType,
                                    tup.partnerEntityType,
                                    tup.partnerEntityId,
                                    tup.partnerCompId,
                                    tup.partnerAsymId,
                                    tup.partnerSeqId,
                                    tup.partnerAtomId,
                                    float(tup.bondDistance) if tup.bondDistance else None,
                                )
                            )
                else:
                    # Calculate ligand - target interactions
                    lArr = np.array(ligXyzL, order="F")
                    distance, index = tree.query(lArr, distance_upper_bound=distLimit)
                    logger.debug("%s lig asymId %s distance %r  index %r", entryId, asymId, distance, index)
                    for ligIndex, (dist, ind) in enumerate(zip(distance, index)):
                        if dist == np.inf:
                            continue
                        # ----
                        ligandTargetInstanceD.setdefault(asymId, []).append(
                            LigandTargetInstance(
                                ligandRefD[asymId][ligIndex].compId,
                                ligandRefD[asymId][ligIndex].atomName,
                                "non-bonded",
                                targetRefL[ind].entityType,
                                targetRefL[ind].entityId,
                                targetRefL[ind].compId,
                                targetRefL[ind].asymId,
                                targetRefL[ind].seqId,
                                targetRefL[ind].atomName,
                                dist,
                            )
                        )
                        # ----
                    if not (asymId in ligandTargetInstanceD and len(ligandTargetInstanceD[asymId])):
                        logger.debug("%s no neighbors for ligand asymId %s within %.2f", entryId, asymId, distLimit)
                        continue
            #
            # re-sort by distance -
            for asymId in ligandTargetInstanceD:
                ligandTargetInstanceD[asymId] = sorted(ligandTargetInstanceD[asymId], key=itemgetter(-1))
            # ------
        except Exception as e:
            logger.exception("Failing for %r with %r", dataContainer.getName() if dataContainer else None, str(e))
        #
        logger.info("Completed %s at %s (%.4f seconds)", dataContainer.getName(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return ligandTargetInstanceD


def getNeighborsSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(NeighborUtilsTests("testGetNeighbors"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = getNeighborsSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
