##
# File: EntityClusterDataPrep.py
# Date: 21-Mar-2017
#
# Update:
#    18-May-2017  jdw adapt for os config
#    19-Apr-2018  jdw adding wos package
#    15-Jul-2018. jdw migrated to latest entity sequence clusters.
#    28-Jul-2018  jdw generalize path to mmseqs, add as input to constructor
##
import os
import pickle
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("root")

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))


class FormatEntityClusters(object):
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

    def __init__(self, clusterPath):
        self.__clusterPath = clusterPath if clusterPath else "../mmseqs"
        self.__pickleProtocol = 0
        # clusters-by-entity-
        self.__bcData = [
            ("100", "clusters-by-entity-100.txt"),
            ("30", "clusters-by-entity-30.txt"),
            ("40", "clusters-by-entity-40.txt"),
            ("50", "clusters-by-entity-50.txt"),
            ("70", "clusters-by-entity-70.txt"),
            ("90", "clusters-by-entity-90.txt"),
            ("95", "clusters-by-entity-95.txt"),
        ]

    def __serialize(self, filePath, iD):
        try:
            with open(filePath, "wb") as fb:
                pickle.dump(iD, fb, self.__pickleProtocol)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def __readFile(self, fPath):
        ifh = open(fPath, "r")
        iClust = 1
        clusterD = {}
        for line in ifh:
            fields = str(line[:-1]).split()
            clusterD[iClust] = fields
            iClust += 1
        ifh.close()
        #
        clusterD[iClust] = ["1D01_2"]
        return clusterD

    def __makeIdChainDict(self, clusterD):
        idChD = {}
        for k, vL in clusterD.items():
            for idChain in vL:
                idChD[idChain] = k
        return idChD

    def __makeIdEntityDict(self, clusterD):
        idEntityD = {}
        for k, vL in clusterD.items():
            for idEntity in vL:
                idEntityD[idEntity] = k
        return idEntityD

    def __build(self):
        mD = {}
        sL = ["100", "95", "90", "70", "50", "30"]
        for simV, fn in self.__bcData:
            clusterD = self.__readFile(os.path.join(self.__clusterPath, fn))
            mD[simV] = self.__makeIdEntityDict(clusterD)
        kL = sorted(mD["100"].keys())
        #
        rD = {}
        for k in kL:
            tL = []
            for tS in sL:
                if k in mD[tS]:
                    tL.append(str(mD[tS][k]))
                else:
                    logger.info("Missing value for sim %s case %s\n", tS, k)
                    tL.append("na")

            rD[k] = tuple(tL)
        return rD

    def export(self, fPath, fType="pic"):
        rD = self.__build()
        if fType == "pic":
            self.__serialize(os.path.join(fPath), rD)
        else:
            ofh = open(os.path.join(fPath), "w")
            kys = sorted(rD.keys())
            ofh.write("PDB_ID,ENTITY_ID,100,95,90,70,50,30\n")
            for k in kys:
                v = rD[k]
                ff = k.split("_")
                ff.extend(v)
                ofh.write("%s\n" % ",".join(ff))
            ofh.close()


if __name__ == "__main__":
    clusterTopPath = "../mmseqs"
    clustPath = os.path.join(clusterTopPath, "entity-clust-combined.pic")
    fbc = FormatEntityClusters(clusterPath=clusterTopPath)
    fbc.export(clustPath, fType="pic")
    #
