##
# File:    testGeneCoverage.py
# Author:  J. Westbrook
# Date:    10-May-2021
#
# Updates:
#
##
"""
Evaluate selected gene coverage in PDB for human and model organism
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

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.FastaUtil import FastaUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.TimeUtil import TimeUtil
from rcsb.utils.seq.UniProtUtils import UniProtUtils
from rcsb.utils.seqalign.MMseqsUtils import MMseqsUtils
from rcsb.workflow.targets.ProteinTargetSequenceWorkflow import ProteinTargetSequenceWorkflow

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class GeneCoverageTests(unittest.TestCase):
    def __init__(self, methodName="runTest"):
        super(GeneCoverageTests, self).__init__(methodName)
        self.__verbose = True

    def setUp(self):
        #
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        configPath = os.path.join(TOPDIR, "rcsb", "mock-data", "config", "dbload-setup-example.yml")
        configName = "site_info_remote_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        #
        self.__resourceName = "MONGO_DB"
        self.__workPath = os.path.join(TOPDIR, "CACHE", "gene-coverage")
        #
        self.__testEntryCacheKwargs = {"fmt": "json", "indent": 3}
        self.__objectLimitTest = None
        #
        self.__cdxGeneList = [
            "ABL1",
            "ACVR1B",
            "AKT1",
            "AKT2",
            "AKT3",
            "ALK",
            "ALOX12B",
            "AMER1",
            "APC",
            "AR",
            "ARAF",
            "ARFRP1",
            "ARID1A",
            "ASXL1",
            "ATM",
            "ATR",
            "ATRX",
            "AURKA",
            "AURKB",
            "AXIN1",
            "AXL",
            "BAP1",
            "BARD1",
            "BCL2",
            "BCL2L1",
            "BCL2L2",
            "BCL6",
            "BCOR",
            "BCORL1",
            "BRAF",
            "BRCA1",
            "BRCA2",
            "BRD4",
            "BRIP1",
            "BTG1",
            "BTG2",
            "BTK",
            "C11orf30",
            "CALR",
            "CARD11",
            "CASP8",
            "CBFB",
            "CBL",
            "CCND1",
            "CCND2",
            "CCND3",
            "CCNE1",
            "CD22",
            "CD274",
            "CD70",
            "CD79A",
            "CD79B",
            "CDC73",
            "CDH1",
            "CDK12",
            "CDK4",
            "CDK6",
            "CDK8",
            "CDKN1A",
            "CDKN1B",
            "CDKN2A",
            "CDKN2B",
            "CDKN2C",
            "CEBPA",
            "CHEK1",
            "CHEK2",
            "CIC",
            "CREBBP",
            "CRKL",
            "CSF1R",
            "CSF3R",
            "CTCF",
            "CTNNA1",
            "CTNNB1",
            "CUL3",
            "CUL4A",
            "CXCR4",
            "CYP17A1",
            "DAXX",
            "DDR1",
            "DDR2",
            "DIS3",
            "DNMT3A",
            "DOT1L",
            "EED",
            "EGFR",
            "EP300",
            "EPHA3",
            "EPHB1",
            "EPHB4",
            "ERBB2",
            "ERBB3",
            "ERBB4",
            "ERCC4",
            "ERG",
            "ERRFI1",
            "ESR1",
            "EZH2",
            "FAM46C",
            "FANCA",
            "FANCC",
            "FANCG",
            "FANCL",
            "FAS",
            "FBXW7",
            "FGF10",
            "FGF12",
            "FGF14",
            "FGF19",
            "FGF23",
            "FGF3",
            "FGF4",
            "FGF6",
            "FGFR1",
            "FGFR2",
            "FGFR3",
            "FGFR4",
            "FH",
            "FLCN",
            "FLT1",
            "FLT3",
            "FOXL2",
            "FUBP1",
            "GABRA6",
            "GATA3",
            "GATA4",
            "GATA6",
            "GID4",
            "GNA11",
            "GNA13",
            "GNAQ",
            "GNAS",
            "GRM3",
            "GSK3B",
            "H3F3A",
            "HDAC1",
            "HGF",
            "HNF1A",
            "HRAS",
            "HSD3B1",
            "ID3",
            "IDH1",
            "IDH2",
            "IGF1R",
            "IKBKE",
            "IKZF1",
            "INPP4B",
            "IRF2",
            "IRF4",
            "IRS2",
            "JAK1",
            "JAK2",
            "JAK3",
            "JUN",
            "KDM5A",
            "KDM5C",
            "KDM6A",
            "KDR",
            "KEAP1",
            "KEL",
            "KIT",
            "KLHL6",
            "KMT2A",
            "KMT2D",
            "KRAS",
            "LTK",
            "LYN",
            "MAF",
            "MAP2K1",
            "MAP2K2",
            "MAP2K4",
            "MAP3K1",
            "MAP3K13",
            "MAPK1",
            "MCL1",
            "MDM2",
            "MDM4",
            "MED12",
            "MEF2B",
            "MEN1",
            "MERTK",
            "MET",
            "MITF",
            "MKNK1",
            "MLH1",
            "MPL",
            "MRE11A",
            "MSH2",
            "MSH3",
            "MSH6",
            "MST1R",
            "MTAP",
            "MTOR",
            "MUTYH",
            "MYC",
            "MYCL",
            "MYCN",
            "MYD88",
            "NBN",
            "NF1",
            "NF2",
            "NFE2L2",
            "NFKBIA",
            "NKX2-1",
            "NOTCH1",
            "NOTCH2",
            "NOTCH3",
            "NPM1",
            "NRAS",
            "NT5C2",
            "NTRK1",
            "NTRK2",
            "NTRK3",
            "P2RY8",
            "PALB2",
            "PARK2",
            "PARP1",
            "PARP2",
            "PARP3",
            "PAX5",
            "PBRM1",
            "PDCD1",
            "PDCD1LG2",
            "PDGFRA",
            "PDGFRB",
            "PDK1",
            "PIK3C2B",
            "PIK3C2G",
            "PIK3CA",
            "PIK3CB",
            "PIK3R1",
            "PIM1",
            "PMS2",
            "POLD1",
            "POLE",
            "PPARG",
            "PPP2R1A",
            "PPP2R2A",
            "PRDM1",
            "PRKAR1A",
            "PRKCI",
            "PTCH1",
            "PTEN",
            "PTPN11",
            "PTPRO",
            "QKI",
            "RAC1",
            "RAD21",
            "RAD51",
            "RAD51B",
            "RAD51C",
            "RAD51D",
            "RAD52",
            "RAD54L",
            "RAF1",
            "RARA",
            "RB1",
            "RBM10",
            "REL",
            "RET",
            "RICTOR",
            "RNF43",
            "ROS1",
            "RPTOR",
            "SDHA",
            "SDHB",
            "SDHC",
            "SDHD",
            "SETD2",
            "SF3B1",
            "SGK1",
            "SMAD2",
            "SMAD4",
            "SMARCA4",
            "SMARCB1",
            "SMO",
            "SNCAIP",
            "SOCS1",
            "SOX2",
            "SOX9",
            "SPEN",
            "SPOP",
            "SRC",
            "STAG2",
            "STAT3",
            "STK11",
            "SUFU",
            "SYK",
            "TBX3",
            "TEK",
            "TET2",
            "TGFBR2",
            "TIPARP",
            "TNFAIP3",
            "TNFRSF14",
            "TP53",
            "TSC1",
            "TSC2",
            "TYRO3",
            "U2AF1",
            "VEGFA",
            "VHL",
            "WHSC1",
            "WHSC1L1",
            "WT1",
            "XPO1",
            "XRCC2",
            "ZNF217",
            "ZNF703",
        ]
        self.__hemeGeneList = [
            "ABL1",
            "ACTB",
            "AKT1",
            "AKT2",
            "AKT3",
            "ALK",
            "AMER1",
            "APC",
            "APH1A",
            "AR",
            "ARAF",
            "ARFRP1",
            "ARHGAP26",
            "ARID1A",
            "ARID2",
            "ASMTL",
            "ASXL1",
            "ATM",
            "ATR",
            "ATRX",
            "AURKA",
            "AURKB",
            "AXIN1",
            "AXL",
            "B2M",
            "BAP1",
            "BARD1",
            "BCL10",
            "BCL11B",
            "BCL2",
            "BCL2L2",
            "BCL6",
            "BCL7A",
            "BCOR",
            "BCORL1",
            "BIRC3",
            "BLM",
            "BRAF",
            "BRCA1",
            "BRCA2",
            "BRD4",
            "BRIP1",
            "BRSK1",
            "BTG2",
            "BTK",
            "BTLA",
            "C11orf30",
            "CAD",
            "CALR",
            "CARD11",
            "CBFB",
            "CBL",
            "CCND1",
            "CCND2",
            "CCND3",
            "CCNE1",
            "CCT6B",
            "CD22",
            "CD274",
            "CD36",
            "CD58",
            "CD70",
            "CD79A",
            "CD79B",
            "CDC73",
            "CDH1",
            "CDK12",
            "CDK4",
            "CDK6",
            "CDK8",
            "CDKN1B",
            "CDKN2A",
            "CDKN2B",
            "CDKN2C",
            "CEBPA",
            "CHD2",
            "CHEK1",
            "CHEK2",
            "CIC",
            "CIITA",
            "CKS1B",
            "CPS1",
            "CREBBP",
            "CRKL",
            "CRLF2",
            "CSF1R",
            "CSF3R",
            "CTCF",
            "CTNNA1",
            "CTNNB1",
            "CUX1",
            "CXCR4",
            "DAXX",
            "DDR2",
            "DDX3X",
            "DNM2",
            "DNMT3A",
            "DOT1L",
            "DTX1",
            "DUSP2",
            "DUSP9",
            "EBF1",
            "ECT2L",
            "EED",
            "EGFR",
            "ELP2",
            "EP300",
            "EPHA3",
            "EPHA5",
            "EPHA7",
            "EPHB1",
            "ERBB2",
            "ERBB3",
            "ERBB4",
            "ERG",
            "ESR1",
            "ETS1",
            "ETV6",
            "EXOSC6",
            "EZH2",
            "FAF1",
            "FAM46C",
            "FANCA",
            "FANCC",
            "FANCD2",
            "FANCE",
            "FANCF",
            "FANCG",
            "FANCL",
            "FAS",
            "FBXO11",
            "FBXO31",
            "FBXW7",
            "FGF10",
            "FGF14",
            "FGF19",
            "FGF23",
            "FGF3",
            "FGF4",
            "FGF6",
            "FGFR1",
            "FGFR2",
            "FGFR3",
            "FGFR4",
            "FHIT",
            "FLCN",
            "FLT1",
            "FLT3",
            "FLT4",
            "FLYWCH1",
            "FOXL2",
            "FOXO1",
            "FOXO3",
            "FOXP1",
            "FRS2",
            "GADD45B",
            "GATA1",
            "GATA2",
            "GATA3",
            "GID4",
            "GNA11",
            "GNA12",
            "GNA13",
            "GNAQ",
            "GNAS",
            "GPR124",
            "GRIN2A",
            "GSK3B",
            "GTSE1",
            "HDAC1",
            "HDAC4",
            "HDAC7",
            "HGF",
            "HIST1H1C",
            "HIST1H1D",
            "HIST1H1E",
            "HIST1H2AC",
            "HIST1H2AG",
            "HIST1H2AL",
            "HIST1H2AM",
            "HIST1H2BC",
            "HIST1H2BJ",
            "HIST1H2BK",
            "HIST1H2BO",
            "HIST1H3B",
            "HNF1A",
            "HRAS",
            "HSP90AA1",
            "ICK",
            "ID3",
            "IDH1",
            "IDH2",
            "IGF1R",
            "IKBKE",
            "IKZF1",
            "IKZF2",
            "IKZF3",
            "IL7R",
            "INHBA",
            "INPP4B",
            "INPP5D",
            "IRF1",
            "IRF4",
            "IRF8",
            "IRS2",
            "JAK1",
            "JAK2",
            "JAK3",
            "JARID2",
            "JUN",
            "KAT6A",
            "KDM2B",
            "KDM4C",
            "KDM5A",
            "KDM5C",
            "KDM6A",
            "KDR",
            "KEAP1",
            "KIT",
            "KLHL6",
            "KMT2A",
            "KMT2C",
            "KMT2D",
            "KRAS",
            "LEF1",
            "LRP1B",
            "LRRK2",
            "MAF",
            "MAFB",
            "MAGED1",
            "MALT1",
            "MAP2K1",
            "MAP2K2",
            "MAP2K4",
            "MAP3K1",
            "MAP3K14",
            "MAP3K6",
            "MAP3K7",
            "MAPK1",
            "MCL1",
            "MDM2",
            "MDM4",
            "MED12",
            "MEF2B",
            "MEF2C",
            "MEN1",
            "MET",
            "MIB1",
            "MITF",
            "MKI67",
            "MLH1",
            "MPL",
            "MRE11A",
            "MSH2",
            "MSH3",
            "MSH6",
            "MTOR",
            "MUTYH",
            "MYC",
            "MYCL",
            "MYCN",
            "MYD88",
            "MYO18A",
            "NCOR2",
            "NCSTN",
            "NF1",
            "NF2",
            "NFE2L2",
            "NFKBIA",
            "NKX2-1",
            "NOD1",
            "NOTCH1",
            "NOTCH2",
            "NPM1",
            "NRAS",
            "NSD1",
            "NT5C2",
            "NTRK1",
            "NTRK2",
            "NTRK3",
            "NUP93",
            "NUP98",
            "P2RY8",
            "PAG1",
            "PAK3",
            "PALB2",
            "PASK",
            "PAX5",
            "PBRM1",
            "PC",
            "PCBP1",
            "PCLO",
            "PDCD1",
            "PDCD11",
            "PDCD1LG2",
            "PDGFRA",
            "PDGFRB",
            "PDK1",
            "PHF6",
            "PIK3CA",
            "PIK3CG",
            "PIK3R1",
            "PIK3R2",
            "PIM1",
            "PLCG2",
            "POT1",
            "PPP2R1A",
            "PRDM1",
            "PRKAR1A",
            "PRKDC",
            "PRSS8",
            "PTCH1",
            "PTEN",
            "PTPN11",
            "PTPN2",
            "PTPN6",
            "PTPRO",
            "RAD21",
            "RAD50",
            "RAD51",
            "RAF1",
            "RARA",
            "RASGEF1A",
            "RB1",
            "RELN",
            "RET",
            "RHOA",
            "RICTOR",
            "RNF43",
            "ROS1",
            "RPTOR",
            "RUNX1",
            "S1PR2",
            "SDHA",
            "SDHB",
            "SDHC",
            "SDHD",
            "SERP2",
            "SETBP1",
            "SETD2",
            "SF3B1",
            "SGK1",
            "SMAD2",
            "SMAD4",
            "SMARCA1",
            "SMARCA4",
            "SMARCB1",
            "SMC1A",
            "SMC3",
            "SMO",
            "SOCS1",
            "SOCS2",
            "SOCS3",
            "SOX10",
            "SOX2",
            "SPEN",
            "SPOP",
            "SRC",
            "SRSF2",
            "STAG2",
            "STAT3",
            "STAT4",
            "STAT5A",
            "STAT5B",
            "STAT6",
            "STK11",
            "SUFU",
            "SUZ12",
            "TAF1",
            "TBL1XR1",
            "TCF3",
            "TCL1A",
            "TET2",
            "TGFBR2",
            "TLL2",
            "TMEM30A",
            "TMSB4X",
            "TNFAIP3",
            "TNFRSF11A",
            "TNFRSF14",
            "TNFRSF17",
            "TOP1",
            "TP53",
            "TP63",
            "TRAF2",
            "TRAF3",
            "TRAF5",
            "TSC1",
            "TSC2",
            "TSHR",
            "TUSC3",
            "TYK2",
            "U2AF1",
            "U2AF2",
            "VHL",
            "WDR90",
            "WHSC1",
            "WISP3",
            "WT1",
            "XBP1",
            "XPO1",
            "YY1AP1",
            "ZMYM3",
            "ZNF217",
            "ZNF24",
            "ZNF703",
            "ZRSR2",
        ]
        #
        self.__timeOut = 3600
        self.__identityCutoff = 0.90
        self.__seqDbTopPath = os.path.join(self.__workPath, "db")
        self.__seqDataTupL = [
            (os.path.join(self.__workPath, "FASTA", "pdb-protein-targets.fa"), os.path.join(self.__workPath, "FASTA", "pdb-protein-targets-taxon.tdd"), "pdbpr", 0),
            (os.path.join(self.__workPath, "special-sequences.fa"), None, "special", 100),
        ]

        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testAnalMatchResults(self):
        #
        # get gene correspondences
        mU = MarshalUtil()
        geneUniProtMappingPath = os.path.join(self.__workPath, "gene-uniprot-mapping-reviewed.json")
        gD = mU.doImport(geneUniProtMappingPath, fmt="json")
        idD = {}
        for gene, unpIdList in gD.items():
            for unpId in unpIdList:
                idD.setdefault(unpId, []).append(gene)
        #
        for level in ["95", "90", "70", "50", "40"]:
            # Get sequence matching data
            resultDirPath = os.path.join(self.__workPath, "search-results-db", "special-results-%s.json" % level)
            mDL = mU.doImport(resultDirPath, fmt="json")
            logger.info("%s Total match length (%d)", level, len(mDL))
            #
            cdxMatchHumanD = {}
            hemeMatchHumanD = {}
            cdxMatchOtherD = {}
            hemeMatchOtherD = {}
            cdxMatchMouseD = {}
            hemeMatchMouseD = {}
            for mD in mDL:
                unpIdMatch = mD["query"].split("|")[0].strip()
                if unpIdMatch not in idD:
                    logger.info("Missing id %s", unpIdMatch)
                    continue
                #
                if len(idD[unpIdMatch]) < 1:
                    logger.info("Missing genes for %s", unpIdMatch)
                for gene in idD[unpIdMatch]:
                    if gene not in self.__hemeGeneList and gene not in self.__cdxGeneList:
                        logger.info("Missing gene %r for %r", gene, unpIdMatch)

                pdbEntityId = mD["target"].split("|")[0]
                taxId = mD["targetTaxId"]
                taxName = mD["targetTaxName"]
                seqIdentity = mD["sequenceIdentity"]
                if taxId == 9606:
                    for gene in idD[unpIdMatch]:
                        if gene in self.__hemeGeneList:
                            hemeMatchHumanD.setdefault(gene, []).append((pdbEntityId, taxId, taxName, unpIdMatch, seqIdentity))
                        if gene in self.__cdxGeneList:
                            cdxMatchHumanD.setdefault(gene, []).append((pdbEntityId, taxId, taxName, unpIdMatch, seqIdentity))
                elif taxId == 10090:
                    for gene in idD[unpIdMatch]:
                        if gene in self.__hemeGeneList:
                            hemeMatchMouseD.setdefault(gene, []).append((pdbEntityId, taxId, taxName, unpIdMatch, seqIdentity))
                        if gene in self.__cdxGeneList:
                            cdxMatchMouseD.setdefault(gene, []).append((pdbEntityId, taxId, taxName, unpIdMatch, seqIdentity))
                else:
                    for gene in idD[unpIdMatch]:
                        if gene in self.__hemeGeneList:
                            hemeMatchOtherD.setdefault(gene, []).append((pdbEntityId, taxId, taxName, unpIdMatch, seqIdentity))
                        if gene in self.__cdxGeneList:
                            cdxMatchOtherD.setdefault(gene, []).append((pdbEntityId, taxId, taxName, unpIdMatch, seqIdentity))
                #
            logger.info("%s Cdx  human  %d/%d", level, len(cdxMatchHumanD), len(self.__cdxGeneList))
            logger.info("%s Heme human  %d/%d", level, len(hemeMatchHumanD), len(self.__hemeGeneList))
            logger.info("%s Cdx  mouse  %d/%d", level, len(cdxMatchMouseD), len(self.__cdxGeneList))
            logger.info("%s Heme mouse  %d/%d", level, len(hemeMatchMouseD), len(self.__hemeGeneList))
            logger.info("%s Cdx  other  %d/%d", level, len(cdxMatchOtherD), len(self.__cdxGeneList))
            logger.info("%s Heme other  %d/%d", level, len(hemeMatchOtherD), len(self.__hemeGeneList))
            #
            missGeneL = []
            for gene in self.__hemeGeneList:
                if (gene not in hemeMatchHumanD) and (gene not in hemeMatchMouseD) and (gene not in hemeMatchOtherD):
                    missGeneL.append(gene)
            logger.info(
                "%s coverage %.2f unmatched heme genes (%d) %r", level, 100.0 * (len(self.__hemeGeneList) - len(missGeneL)) / len(self.__hemeGeneList), len(missGeneL), missGeneL
            )

            missGeneL = []
            for gene in self.__cdxGeneList:
                if gene not in cdxMatchHumanD and gene not in cdxMatchMouseD and gene not in cdxMatchOtherD:
                    missGeneL.append(gene)
            logger.info(
                "%s coverage %.2f unmatched cdx genes (%d) %r", level, 100.0 * (len(self.__cdxGeneList) - len(missGeneL)) / len(self.__cdxGeneList), len(missGeneL), missGeneL
            )
            #
            tL = []
            for gene, mtupL in cdxMatchHumanD.items():
                for mtup in mtupL:
                    tL.append({"Group": "cdx", "gene": gene, "PDB_Entity_ID": mtup[0], "Taxonomy": mtup[2], "UniProt_ID": mtup[3], "Identity": mtup[4]})
            ok = mU.doExport(os.path.join(self.__workPath, "cdx-matches-human-%s.csv" % level), tL, fmt="csv")
            #
            tL = []
            for gene, mtupL in hemeMatchHumanD.items():
                for mtup in mtupL:
                    tL.append({"Group": "heme", "gene": gene, "PDB_Entity_ID": mtup[0], "Taxonomy": mtup[2], "UniProt_ID": mtup[3], "Identity": mtup[4]})
            ok = mU.doExport(os.path.join(self.__workPath, "heme-matches-human-%s.csv" % level), tL, fmt="csv")
            # ----
            tL = []
            for gene, mtupL in cdxMatchMouseD.items():
                for mtup in mtupL:
                    tL.append({"Group": "cdx", "gene": gene, "PDB_Entity_ID": mtup[0], "Taxonomy": mtup[2], "UniProt_ID": mtup[3], "Identity": mtup[4]})
            ok = mU.doExport(os.path.join(self.__workPath, "cdx-matches-mouse-%s.csv" % level), tL, fmt="csv")
            #
            tL = []
            for gene, mtupL in hemeMatchMouseD.items():
                for mtup in mtupL:
                    tL.append({"Group": "heme", "gene": gene, "PDB_Entity_ID": mtup[0], "Taxonomy": mtup[2], "UniProt_ID": mtup[3], "Identity": mtup[4]})
            ok = mU.doExport(os.path.join(self.__workPath, "heme-matches-mouse-%s.csv" % level), tL, fmt="csv")
            # -----
            tL = []
            for gene, mtupL in cdxMatchOtherD.items():
                for mtup in mtupL:
                    tL.append({"Group": "cdx", "gene": gene, "PDB_Entity_ID": mtup[0], "Taxonomy": mtup[2], "UniProt_ID": mtup[3], "Identity": mtup[4]})
            ok = mU.doExport(os.path.join(self.__workPath, "cdx-matches-other-%s.csv" % level), tL, fmt="csv")
            #
            tL = []
            for gene, mtupL in hemeMatchOtherD.items():
                for mtup in mtupL:
                    tL.append({"Group": "heme", "gene": gene, "PDB_Entity_ID": mtup[0], "Taxonomy": mtup[2], "UniProt_ID": mtup[3], "Identity": mtup[4]})
            ok = mU.doExport(os.path.join(self.__workPath, "heme-matches-other-%s.csv" % level), tL, fmt="csv")

    def testACreateDatabases(self):
        """Test case -  for PDB proteins, create sequence search database from fasta file"""
        try:
            mU = MarshalUtil()
            mU.mkdir(os.path.join(self.__workPath, "db"))
            mmS = MMseqsUtils(cachePath=self.__workPath)
            for (fastaPath, taxonPath, dbName, _) in self.__seqDataTupL:
                ok = mmS.createSearchDatabase(fastaPath, self.__seqDbTopPath, dbName, timeOut=self.__timeOut, verbose=True)
                self.assertTrue(ok)
                if taxonPath:
                    ok = mmS.createTaxonomySearchDatabase(taxonPath, self.__seqDbTopPath, dbName, timeOut=self.__timeOut)
                    self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testSearchDatabaseDatabase(self):
        """Test case -  search similar sequences (database) workflow  - special vs PDB"""
        try:
            resultDirPath = os.path.join(self.__workPath, "search-results-db")
            mU = MarshalUtil()
            mU.mkdir(resultDirPath)
            #
            for identityCutoff in [0.95, 0.90, 0.70, 0.50, 0.40]:
                mmS = MMseqsUtils(cachePath=self.__workPath)
                qTup = self.__seqDataTupL[0]
                for (_, taxonPath, dbName, minMatch) in self.__seqDataTupL[1:]:
                    resultPath = os.path.join(resultDirPath, dbName + "-results-%d.json" % int(100 * identityCutoff))
                    logger.info("Searching at level %.2f", identityCutoff)
                    mmS = MMseqsUtils(cachePath=self.__workPath)
                    ok = mmS.searchDatabase(dbName, self.__seqDbTopPath, qTup[2], resultPath, minSeqId=identityCutoff, timeOut=self.__timeOut, sensitivity=5.0)
                    self.assertTrue(ok)
                    #
                    if taxonPath:
                        mL = mmS.getMatchResults(resultPath, taxonPath, useTaxonomy=True, misMatchCutoff=-1, sequenceIdentityCutoff=identityCutoff)
                    else:
                        mL = mmS.getMatchResults(resultPath, None, useTaxonomy=False, misMatchCutoff=-1, sequenceIdentityCutoff=identityCutoff)

                    logger.info("Search result %r (%d)", dbName, len(mL))
                    self.assertGreaterEqual(len(mL), minMatch)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testProteinEntityData(self):
        """Test case - export protein entity sequence Fasta, taxonomy, and sequence details"""
        try:
            ptsW = ProteinTargetSequenceWorkflow(self.__cfgOb, self.__workPath)
            ok = ptsW.exportPDBProteinEntityFasta(minSeqLen=9)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testGeneSelectionMakeFasta(self):
        """Make a FASTA file for the uniprot accessions covering the gene list"""
        geneUniProtMappingPath = os.path.join(self.__workPath, "gene-uniprot-mapping-reviewed.json")
        mU = MarshalUtil()
        gUD = mU.doImport(geneUniProtMappingPath, fmt="json")
        logger.info("mapping list (%d)", len(gUD))
        if False:
            geneUniProtMappingPath = os.path.join(self.__workPath, "miss-gene-uniprot-mapping-reviewed.json")
            mU = MarshalUtil()
            mgUD = mU.doImport(geneUniProtMappingPath, fmt="json")
            gUD.update(mgUD)
            logger.info("mapping list plus missed (%d)", len(gUD))
            #
        unpIdList = []
        for _, tidList in gUD.items():
            unpIdList.extend(tidList)
        #
        unpIdList = sorted(set(unpIdList))
        chunkSize = 50
        numChunks = int(len(unpIdList) / chunkSize)
        logger.info("UniProt fetch input list (%d)", len(unpIdList))
        missCL = []
        sD = {}
        uUtils = UniProtUtils(saveText=False)
        for ii, chunk in enumerate(self.__chunker(unpIdList, chunkSize)):
            # if ii > 2:
            #    break
            ok, tD = uUtils.fetchSequenceList(chunk)
            logger.info("%4d/%4d (%r) fetch length (%d)", ii, numChunks, ok, len(tD))
            if not ok:
                # try again
                ok, tD = uUtils.fetchSequenceList(chunk)
            if ok:
                sD.update(tD)
            else:
                logger.info("%4d/%4d FAILED %r", ii, numChunks, chunk)
                missCL.extend(chunk)
            #
        mU = MarshalUtil()
        mU.doExport(os.path.join(self.__workPath, "data-sequences.json"), sD, fmt="json", indent=3)
        logger.info("Missed ids %r", missCL)
        #
        faPath = os.path.join(self.__workPath, "special-sequences.fa")
        faUtil = FastaUtil()
        faUtil.writeFasta(faPath, sD, makeComment=True)

    def __chunker(self, seq, size):
        return (seq[pos : pos + size] for pos in range(0, len(seq), size))

    def testFetchUniprotIdsForGenes(self):
        """Test case - find genes in UniProt"""
        try:
            # Get the unique list of genes -
            #
            gUD = {}
            missGeneL = []
            logger.info("cdx   genes %d", len(self.__cdxGeneList))
            logger.info("heme  genes %d", len(self.__hemeGeneList))
            geneList = sorted(set((self.__cdxGeneList + self.__hemeGeneList)))
            logger.info("total genes %d", len(geneList))
            uUtils = UniProtUtils(saveText=False)
            for ii, gene in enumerate(geneList):
                uniProtIdList, retCode = uUtils.doGeneLookup(gene, 9606, reviewed=True)
                if retCode != 200:
                    logger.info("return code for %r is %r", gene, retCode)
                if uniProtIdList:
                    gUD.setdefault(gene, []).extend(uniProtIdList)
                else:
                    logger.info("No UniProt result for gene %r", gene)
                    missGeneL.append(gene)
                logger.info("%4d/%4d %r (%d)", ii, len(geneList), gene, len(uniProtIdList) if uniProtIdList else 0)
            logger.info("Unmatched genes %r", missGeneL)
            #
            geneUniProtMappingPath = os.path.join(self.__workPath, "gene-uniprot-mapping-reviewed.json")
            mU = MarshalUtil()
            ok = mU.doExport(geneUniProtMappingPath, gUD, fmt="json", indent=3)
            self.assertTrue(ok)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchUniprotIdsForGenesMissed(self):
        """Test case - find genes in UniProt"""
        try:
            # Get the unique list of genes -
            #
            gUD = {}
            missGeneL = []
            geneList = ["GSK3B", "SETD2", "TMSB4X"]
            logger.info("total genes %d", len(geneList))
            uUtils = UniProtUtils(saveText=False)
            for ii, gene in enumerate(geneList):
                uniProtIdList, retCode = uUtils.doGeneLookup(gene, 9606)
                if retCode != 200:
                    logger.info("return code for %r is %r", gene, retCode)
                if uniProtIdList:
                    gUD.setdefault(gene, []).extend(uniProtIdList)
                else:
                    logger.info("No UniProt result for gene %r", gene)
                    missGeneL.append(gene)
                logger.info("%4d/%4d %r (%d)", ii, len(geneList), gene, len(uniProtIdList) if uniProtIdList else 0)
            logger.info("Unmatched genes %r", missGeneL)
            #
            geneUniProtMappingPath = os.path.join(self.__workPath, "miss-gene-uniprot-mapping.json")
            mU = MarshalUtil()
            ok = mU.doExport(geneUniProtMappingPath, gUD, fmt="json", indent=3)
            self.assertTrue(ok)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def geneCoverageSuite():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(GeneCoverageTests("testFetchUniprotIdsForGenes"))
    # suiteSelect.addTest(GeneCoverageTests("testFetchUniprotIdsForGenesMissed"))
    suiteSelect.addTest(GeneCoverageTests("testGeneSelectionMakeFasta"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = geneCoverageSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
