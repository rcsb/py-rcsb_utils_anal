##
# File: EntityClusterDataPrep.py
# Date: 8-Oct-2020
#
##
import os
import logging

from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))


class EntityClusterDataPrep(object):
    """
    Using mmseq2 entity clusters -

    clusters-by-entity-100.txt
    clusters-by-entity-30.txt
    clusters-by-entity-40.txt
    clusters-by-entity-50.txt
    clusters-by-entity-70.txt
    clusters-by-entity-90.txt
    clusters-by-entity-95.txt

    """

    def __init__(self, clusterPath=None, workPath=None):
        self.__workPath = workPath if workPath else "."
        self.__clusterPath = clusterPath if clusterPath else "mmseqs"
        self.__mU = MarshalUtil(workPath=self.__workPath)
        #
        # clusters-by-entity-
        self.__mmseqsData = [
            ("100", "clusters-by-entity-100.txt"),
            ("30", "clusters-by-entity-30.txt"),
            ("40", "clusters-by-entity-40.txt"),
            ("50", "clusters-by-entity-50.txt"),
            ("70", "clusters-by-entity-70.txt"),
            ("90", "clusters-by-entity-90.txt"),
            ("95", "clusters-by-entity-95.txt"),
        ]

    def __serialize(self, filePath, iD, fmt="json"):
        try:
            ok = self.__mU.doExport(filePath, iD, fmt=fmt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __makeIdEntityDict(self, clusterD):
        entityIdD = {}
        for clusterId, entityIdL in clusterD.items():
            for idEntity in entityIdL:
                entityIdD[idEntity] = str(clusterId)
        return entityIdD

    def __build(self):
        entityD = {}
        clusterD = {}
        # levelList = ["100", "95", "90", "70", "50", "30"]
        for level, fileName in self.__mmseqsData:
            pth = os.path.join(self.__clusterPath, fileName)
            cL = self.__mU.doImport(pth, fmt="list")
            cD = {str(ii): line.split() for ii, line in enumerate(cL, 1)}
            logger.info("Cluster level %s distinct members %d", level, len(cD))
            entityD[level] = self.__makeIdEntityDict(cD)
            clusterD[level] = cD
        return {"entityD": entityD, "clusterD": clusterD}
        #

    def fetch(self):
        return self.__build()

    def export(self, filePath, fmt="json"):
        rD = self.__build()
        ok = self.__mU.doExport(filePath, rD, fmt=fmt)
        return ok


if __name__ == "__main__":
    clusterTopPath = os.path.join(TOPDIR, "mock-data", "cluster_data", "mmseqs_clusters_current")
    exportPath = os.path.join(".", "entity-clust-combined.json")
    fbc = EntityClusterDataPrep(clusterPath=clusterTopPath)
    fbc.export(exportPath, fmt="json")
    #
