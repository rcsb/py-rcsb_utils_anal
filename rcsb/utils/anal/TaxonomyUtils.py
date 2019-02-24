##
# File: TaxonomyUtils.py
#
import logging

from ete3 import NCBITaxa

logger = logging.getLogger(__name__)


class TaxonomyUtils(object):

    def __init__(self, **kwargs):
        self.__ncbi = NCBITaxa()

    def reload(self):

        self.__ncbi.update_taxonomy_database()

    def getLineage(self, taxId):
        tL = []
        try:
            tL = self.__ncbi.get_lineage(taxId)
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return tL

    def getLineageNames(self, taxId):
        nmL = []
        try:
            lineage = self.__ncbi.get_lineage(taxId)
            names = self.__ncbi.get_taxid_translator(lineage)
            nmL = [names[taxid] for taxid in lineage]
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return nmL

    def isBacteria(self, taxId):
        try:
            lineage = self.__ncbi.get_lineage(taxId)
            return True if 2 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isEukaryota(self, taxId):
        try:
            lineage = self.__ncbi.get_lineage(taxId)
            return True if 2759 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isVirus(self, taxId):
        try:
            lineage = self.__ncbi.get_lineage(taxId)
            return True if 10239 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isArchaea(self, taxId):
        try:
            lineage = self.__ncbi.get_lineage(taxId)
            return True if 2157 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isOther(self, taxId):
        """ other/synthetic
        """
        try:
            lineage = self.__ncbi.get_lineage(taxId)
            return True if 28384 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False

    def isUnclassified(self, taxId):
        try:
            lineage = self.__ncbi.get_lineage(taxId)
            return True if 12908 in lineage else False
        except Exception as e:
            logger.exception("Failing for taxId %r with %s" % (taxId, str(e)))
        return False
