##
# File:    testLeadingHumanProtein.py
# Author:  J. Westbrook
# Date:    13-Oct-2020
#
# Updates:
#
##
"""
Human protein with only 9606 taxonomy.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import copy
import logging
import os
import time
import unittest

from rcsb.exdb.utils.ObjectExtractor import ObjectExtractor
from rcsb.utils.anal.EntityClusterDataPrep import EntityClusterDataPrep
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.TimeUtil import TimeUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class LeadingHumanProteinTests(unittest.TestCase):
    def __init__(self, methodName="runTest"):
        super(LeadingHumanProteinTests, self).__init__(methodName)
        self.__verbose = True

    def setUp(self):
        #
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        configPath = os.path.join(TOPDIR, "rcsb", "mock-data", "config", "dbload-setup-example.yml")
        self.__clusterTopPath = os.path.join(TOPDIR, "rcsb", "mock-data", "cluster_data", "mmseqs_clusters_current")
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

    def testExtractUniProtDetails(self):
        """Test case - extract selected UniProt details"""
        try:
            obEx = ObjectExtractor(
                self.__cfgOb,
                databaseName="uniprot_exdb",
                collectionName="reference_entry",
                cacheFilePath=os.path.join(self.__workPath, "uniprot-data-details-example.json"),
                useCache=False,
                keyAttribute="uniprot",
                uniqueAttributes=["rcsb_id"],
                cacheKwargs=self.__testEntryCacheKwargs,
                objectLimit=self.__objectLimitTest,
                selectionList=["rcsb_id", "taxonomy_id", "host_taxonomy_id", "names", "gene", "source_scientific", "sequence"],
            )
            eCount = obEx.getCount()
            logger.info("Reference sequence count is %d", eCount)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

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

    def __extractEntityTaxonomyAndRefDetails(self, targetTaxIdL):
        """Test case - extract unique entity source and host taxonomies"""
        try:
            # targetTaxIdL = [9606]
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
                #
                selectionQuery={"entity_poly.rcsb_entity_polymer_type": "Protein", "rcsb_entity_source_organism.ncbi_taxonomy_id": {"$in": targetTaxIdL}},
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
                rcsbId = eD["rcsb_id"]
                try:
                    for tD in eD["rcsb_entity_source_organism"]:
                        taxIdS.add(tD["ncbi_taxonomy_id"])
                    if not set(targetTaxIdL) & taxIdS:
                        logger.error("%s missing any target taxId,  has only %r", rcsbId, taxIdS)
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

    def __getUniProtDetails(self):
        mU = MarshalUtil(workPath=self.__workPath)
        unpFilePath = os.path.join(self.__workPath, "uniprot-data-details-example.json")
        tD = mU.doImport(unpFilePath, fmt="json")
        return tD["uniprot"]

    def __getHostedTaxonomies(self, targetHostTaxId, unpD):
        hostedTaxIdL = []
        for _, uD in unpD.items():
            if "host_taxonomy_id" in uD and uD["host_taxonomy_id"] == targetHostTaxId:
                hostedTaxIdL.append(uD["taxonomy_id"])
        return list(set(hostedTaxIdL))

    def __getEntityDetails(self):
        mU = MarshalUtil(workPath=self.__workPath)
        entityFilePath = os.path.join(self.__workPath, "entity-taxonomy-ref-example.json")
        tD = mU.doImport(entityFilePath, fmt="json")
        return tD["entity"]

    def __getClusterDetails(self):
        ecdp = EntityClusterDataPrep(clusterPath=self.__clusterTopPath)
        return ecdp.fetch()

    def testfindFirst(self):
        for level in ["100", "95", "90"]:
            self.__testFindFirst(level)

    def __testFindFirst(self, level="95"):
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
        # level = "95"
        targetTaxId = 9606
        tU = TimeUtil()
        uniProtD = self.__getUniProtDetails()
        # hostedTaxIdL = self.__getHostedTaxonomies(targetTaxId, uniProtD)
        # logger.info("For %d host taxa count %d", targetTaxId, len(hostedTaxIdL))
        #
        hostedTaxIdL = []
        self.__extractEntityTaxonomyAndRefDetails([targetTaxId])
        entityD = self.__getEntityDetails()
        entryD = self.__getEntryDetails()

        uniProtTaxIdD = {k: v["taxonomy_id"] for k, v in uniProtD.items()}

        # tD = tU.getDateTimeObj(tS)
        # entities assigned to each ref sequence [uniprot]=[(entityId, releaseDT), ]
        refSeqD = {}
        # ref sequences assigned for each entity  [entityId] = [uniprot,...]
        entityRefSeqD = {}
        entityTaxCountD = {}
        entityReleaseDateD = {}
        #
        multiTaxCount = 0
        multiTaxRefSeqCount = 0
        missCount = 0
        badAssignCount = 0
        #
        entriesPerYearMT = {}
        entriesPerYear = {}
        #
        for entityId, eObj in entityD.items():
            # logger.info("eObj %r", eObj)
            # entityId = eObj["rcsb_id"]
            entryId = entityId[:4]
            releaseTs = entryD[entryId]["rcsb_accession_info"]["initial_release_date"]
            releaseDt = tU.getDateTimeObj(releaseTs)
            releaseYear = int(releaseDt.strftime("%Y"))
            # partCount = eObj["rcsb_polymer_entity"]["rcsb_source_part_count"]
            taxCount = eObj["rcsb_polymer_entity"]["rcsb_source_taxonomy_count"]
            entityTaxCountD[entityId] = taxCount
            #
            entriesPerYearMT.setdefault(releaseYear, set()).add(entryId)
            #
            if taxCount > 1:
                multiTaxCount += 1
                continue
            #
            entriesPerYear.setdefault(releaseYear, set()).add(entryId)
            ok = False
            try:
                rDL = eObj["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"]
                for rD in rDL:
                    if rD["database_name"] == "UniProt":
                        uniProtId = rD["database_accession"]
                        if uniProtId in uniProtTaxIdD and uniProtTaxIdD[uniProtId] in [targetTaxId] + hostedTaxIdL:
                            refSeqD.setdefault(uniProtId, []).append((entityId, releaseDt))
                            entityRefSeqD.setdefault(entityId, []).append(uniProtId)
                            ok = True
                        else:
                            badAssignCount += 1
                            uTax = uniProtTaxIdD[uniProtId] if uniProtId in uniProtTaxIdD else None
                            logger.debug("%s has bad uniprot assignment %s utaxId %r", entityId, uniProtId, uTax)
                if ok and taxCount > 1:
                    multiTaxRefSeqCount += 1
            except Exception:
                pass
            if not ok:
                missCount += 1
                logger.debug("Missing Uniprot reference for %r", entityId)
            # only for items in the taxonomy scope
            entityReleaseDateD[entityId] = releaseDt
        #
        logger.info("Level %s", level)
        logger.info("Total polymer entities %d", len(entityD))
        logger.info("Length reference sequences assigned %d", len(refSeqD))
        logger.info("multi-taxonomy entities %d", multiTaxCount)
        logger.info("multi-taxonomy entities with assigned references sequences  %d", multiTaxRefSeqCount)
        logger.info("Unassigned polymer entities %d", missCount)
        logger.info("Bad assignments (switch taxa) %d", badAssignCount)
        #
        logger.debug("entriesPerYearMT %r", entriesPerYearMT)
        logger.debug("entriesPerYear   %r", entriesPerYear)
        entriesPerYearMTCount = {}
        sumEntryMT = 0
        for yr in sorted(entriesPerYearMT):
            nn = len(entriesPerYearMT[yr])
            entriesPerYearMTCount[yr] = nn
            sumEntryMT += nn
        #
        entriesPerYearCount = {}
        entriesPerYearCountCum = {}
        sumEntry = 0
        for yr in sorted(entriesPerYear):
            nn = len(entriesPerYear[yr])
            entriesPerYearCount[yr] = nn
            sumEntry += nn
            entriesPerYearCountCum[yr] = sumEntry
        #
        logger.info("entriesPerYearCount   (%d)  %r", sumEntry, entriesPerYearCount)
        logger.info("entriesPerYearCount   (%d)  %r", sumEntry, entriesPerYearCountCum)
        logger.info("entriesPerYearMTCount (%d)  %r", sumEntryMT, entriesPerYearMTCount)
        #
        entriesPerYearCountRowL = []
        for yr, count in entriesPerYearCount.items():
            entriesPerYearCountRowL.append({"Year": yr, "Entries_Containing_Human_Proteins": count})
        #
        for uniProtId in refSeqD:
            refSeqD[uniProtId] = sorted(refSeqD[uniProtId], key=lambda item: item[1])
        for entityId in entityRefSeqD:
            entityRefSeqD[entityId] = list(set(entityRefSeqD[entityId]))
        #
        retRowFullL, retRowL, firstEntriesPerYearCountRowL = self.clusterFirstRep(level, entityRefSeqD, entityTaxCountD, entityReleaseDateD, entryD, entityD, uniProtD)
        #
        mU = MarshalUtil(workPath=self.__workPath)
        mU.doExport(os.path.join(self.__workPath, "%s-first-human-entity-full.csv" % level), retRowFullL, fmt="csv")
        mU.doExport(os.path.join(self.__workPath, "%s-first-human-entity-abbrev.csv" % level), retRowL, fmt="csv")
        #
        pth = os.path.join(self.__workPath, "human-containing-entries-by-year.csv")
        mU.doExport(pth, entriesPerYearCountRowL, fmt="csv")
        pth = os.path.join(self.__workPath, "human-containing-entries-by-year.md")
        self.__formatTable(pth, entriesPerYearCountRowL)
        #
        pth = os.path.join(self.__workPath, "%s-leading-human-containing-entries-by-year.csv" % level)
        mU.doExport(os.path.join(self.__workPath, "%s-leading-human-containing-entries-by-year.csv" % level), firstEntriesPerYearCountRowL, fmt="csv")
        pth = os.path.join(self.__workPath, "%s-leading-human-containing-entries-by-year.md" % level)
        self.__formatTable(pth, firstEntriesPerYearCountRowL)

    def __formatTable(self, pth, rowDL):
        tsL = []
        rowKeys = list(rowDL[0].keys())
        sepL = [":-------:" for n in range(len(rowKeys))]
        ln = "| " + " | ".join(rowKeys) + " |"
        tsL.append(ln)
        ln = "| " + " | ".join(sepL) + " |"
        tsL.append(ln)
        for rowD in rowDL:
            ln = "| " + " | ".join([str(v) for v in rowD.values()]) + " |"
            tsL.append(ln)
        with open(pth, "w") as ofh:
            for ts in tsL:
                ofh.write("%s\n" % ts)

    def __getEntityAnnotations(self, entityId, entityD):
        descr = ""
        try:
            descr = entityD[entityId]["rcsb_polymer_entity"]["pdbx_description"]
        except Exception:
            pass
        return descr

    def __getEntryAnnotations(self, entryId, entryD):
        title = ""
        descr = ""
        try:
            title = entryD[entryId]["struct"]["title"]
            descr = entryD[entryId]["struct"]["pdbx_descriptor"]
        except Exception:
            pass
        return title, descr

    def __getUniProtAnnotations(self, upIdL, uniProtD, jChar="|"):
        geneL = []
        nameL = []
        for upId in upIdL:
            try:
                for nD in uniProtD[upId]["names"]:
                    if nD["nameType"] == "recommendedName":
                        nameL.append(nD["name"])
                        break
            except Exception:
                pass
            try:
                for gD in uniProtD[upId]["gene"]:
                    if gD["type"] == "primary":
                        geneL.append(gD["name"])
                        break
            except Exception:
                pass
        return jChar.join(nameL), jChar.join(geneL)

    def clusterFirstRep(
        self,
        level,
        entityRefSeqD,
        entityTaxCountD,
        entityReleaseDateD,
        entryD,
        entityD,
        uniProtD,
    ):
        retRowL = []
        retRowFullL = []
        clD = self.__getClusterDetails()
        clusterD = clD["clusterD"]
        firstEntriesPerYear = {}
        #
        # list of reference sequences assigned in this cluster
        clusterRefD = {}
        clusterMembersD = {}
        for clusterId, entityIdL in clusterD[level].items():
            for entityId in entityIdL:
                if entityId not in entityReleaseDateD:
                    # skip if not in taxonomy scope
                    continue
                if entityId in entityRefSeqD:  # ref assigned
                    rL = entityRefSeqD[entityId]
                    clusterRefD.setdefault(clusterId, []).extend(rL)
                else:  # ref unassigned
                    rL = []
                clusterMembersD.setdefault(clusterId, []).append((entityId, entityReleaseDateD[entityId], rL))
        #
        uclusterRefD = {}
        for clusterId, rL in clusterRefD.items():
            uclusterRefD[clusterId] = list(set(rL))

        logger.info("Populated clusters at %s - %d", level, len(clusterMembersD))
        for clusterId, memberList in clusterMembersD.items():
            clusterRefSeqCount = len(uclusterRefD[clusterId]) if clusterId in uclusterRefD else 0
            clusterRefSeqL = uclusterRefD[clusterId] if clusterId in uclusterRefD else []
            memberListS = sorted(memberList, key=lambda item: item[1])
            (firstEntityId, firstRelDT, firstRefL) = memberListS[0]
            firstReleaseYear = int(firstRelDT.strftime("%Y"))
            #
            if clusterRefSeqCount:
                logger.debug("Cluster %-8s members %-5d (of %-5d) reference sequences %2d", clusterId, len(memberList), len(clusterD[level][clusterId]), clusterRefSeqCount)
                # count the instances in the
                #
                refCount = len([True for m in memberList if firstRefL == m[2]])
                logger.debug(">>>  %5s (%d/%d) %s %r", firstEntityId, entityTaxCountD[firstEntityId], refCount, firstRelDT.strftime("%Y-%m-%d"), firstRefL)

                entryTitle, entryDescr = self.__getEntryAnnotations(firstEntityId[:4], entryD)
                entityDescr = self.__getEntityAnnotations(firstEntityId, entityD)
                uName, uGene = self.__getUniProtAnnotations(firstRefL, uniProtD)

                logger.debug("-------- %s ", entryTitle)
                logger.debug("-------- %s ", entryDescr)
                logger.debug("-------- %s ", entityDescr)
                logger.debug("-------- %s ", uName)
                logger.debug("-------- %s ", uGene)
                #
                recD = {
                    "Cluster_ID": clusterId,
                    "Cluster_Members_Total": len(clusterD[level][clusterId]),
                    "Cluster_Members_Human": len(memberList),
                    "Cluster_UniProt_Count": clusterRefSeqCount,
                    "PDB_Entity_ID": firstEntityId,
                    "PDB_Release_Year": int(firstRelDT.strftime("%Y")),
                    "UniProt_IDs": ",".join(firstRefL),
                    "Assign_Count": refCount,
                    "PDB_Struct_title": entryTitle,
                    "PDB_Struct_Descr": entryDescr,
                    "PDB_Entity_Descr": entityDescr,
                    "UniProt_Name": uName,
                    "Uniprot_Gene": uGene,
                }
                retRowL.append(recD)
                retRowFullL.append(recD)
                firstEntriesPerYear.setdefault(firstReleaseYear, set()).add(firstEntityId[:4])
                #
                tL = copy.copy(clusterRefSeqL)
                try:
                    for upId in firstRefL:
                        tL.remove(upId)
                except ValueError:
                    pass
                for (entityId, relDT, refL) in memberListS[1:]:
                    for upId in refL:
                        if upId in tL:
                            refCount = len([True for m in memberList if m[2] == refL])
                            logger.debug(">>> >>>  %5s (%d/%d) %s %r", entityId, entityTaxCountD[entityId], refCount, relDT.strftime("%Y-%m-%d"), refL)
                            entryTitle, entryDescr = self.__getEntryAnnotations(entityId[:4], entryD)
                            entityDescr = self.__getEntityAnnotations(entityId, entityD)
                            uName, uGene = self.__getUniProtAnnotations(refL, uniProtD)
                            #
                            logger.debug("------------------- %s ", entryTitle)
                            logger.debug("------------------- %s ", entryDescr)
                            logger.debug("------------------- %s ", entityDescr)
                            logger.debug("------------------- %s ", uName)
                            logger.debug("------------------- %s ", uGene)
                            # relYear = int(relDT.strftime("%Y"))
                            # firstEntriesPerYear.setdefault(relYear, set()).add(entityId[:4])
                            recD = {
                                "Cluster_ID": clusterId,
                                "Cluster_Members_Total": len(clusterD[level][clusterId]),
                                "Cluster_Members_Human": len(memberList),
                                "Cluster_UniProt_Count": clusterRefSeqCount,
                                "PDB_Entity_ID": entityId,
                                "PDB_Release_Year": int(relDT.strftime("%Y")),
                                "UniProt_IDs": ",".join(refL),
                                "Assign_Count": refCount,
                                "PDB_Struct_title": entryTitle,
                                "PDB_Struct_Descr": entryDescr,
                                "PDB_Entity_Descr": entityDescr,
                                "UniProt_Name": uName,
                                "Uniprot_Gene": uGene,
                            }
                            retRowFullL.append(recD)
                            try:
                                tL.remove(upId)
                            except ValueError:
                                pass
            else:
                # No reference sequeuces -
                logger.debug("cluster %-8s members %-5d (of %-5d) no reference sequences", clusterId, len(memberList), len(clusterD[level][clusterId]))
                logger.debug(">>>  %5s (%d) %s", firstEntityId, entityTaxCountD[firstEntityId], firstRelDT.strftime("%Y-%m-%d"))
                entryTitle, entryDescr = self.__getEntryAnnotations(firstEntityId[:4], entryD)
                entityDescr = self.__getEntityAnnotations(firstEntityId, entityD)
                logger.debug("-------- %s ", entryTitle)
                logger.debug("-------- %s ", entryDescr)
                logger.debug("-------- %s ", entityDescr)
                firstEntriesPerYear.setdefault(firstReleaseYear, set()).add(firstEntityId[:4])
                recD = {
                    "Cluster_ID": clusterId,
                    "Cluster_Members_Total": len(clusterD[level][clusterId]),
                    "Cluster_Members_Human": len(memberList),
                    "Cluster_UniProt_Count": clusterRefSeqCount,
                    "PDB_Entity_ID": firstEntityId,
                    "PDB_Release_Year": int(firstRelDT.strftime("%Y")),
                    "UniProt_IDs": "",
                    "Assign_Count": 0,
                    "PDB_Struct_title": entryTitle,
                    "PDB_Struct_Descr": entryDescr,
                    "PDB_Entity_Descr": entityDescr,
                    "UniProt_Name": "",
                    "Uniprot_Gene": "",
                }
                retRowL.append(recD)
                retRowFullL.append(recD)
        #
        #
        firstEntriesPerYearCount = {}
        sumEntry = 0
        for yr in sorted(firstEntriesPerYear):
            nn = len(firstEntriesPerYear[yr])
            firstEntriesPerYearCount[yr] = nn
            sumEntry += nn
        #
        logger.info("%s firstEntriesPerYearCount   (%d)  %r", level, sumEntry, firstEntriesPerYearCount)
        firstEntriesPerYearCountRowL = []
        for yr, count in firstEntriesPerYearCount.items():
            firstEntriesPerYearCountRowL.append({"Year": yr, "Entries_With_Leading_Human_Protein": count})
        # NOTE - ONLY THE LEADING CLUSTER MEMBERS in firstEntriesPerYearCountRowL...
        return retRowFullL, retRowL, firstEntriesPerYearCountRowL


def leadingHumanProteinSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(LeadingHumanProteinTests("testExtractUniProtDetails"))
    suiteSelect.addTest(LeadingHumanProteinTests("testExtractEntryDetails"))
    suiteSelect.addTest(LeadingHumanProteinTests("testfindFirst"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = leadingHumanProteinSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
