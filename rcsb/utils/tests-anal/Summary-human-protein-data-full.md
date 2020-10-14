
## Summary of Human Protein Data (2020-10-10)

Data tabulated (Oct 7 PDB Release) for human (taxId=9606) or human hosted (221 additional taxons) protein entities.
Humanized antibody proteins are annotated as human taxonomy (taxId=9606).

```none
    Total protein entity count: 75696  (human/human hosted)
    Unique reference sequence (UniProt) assignments: 7838
    Protein entities without reference sequence assignments: 4955
```

Protein entities with multiple distinct taxonomies have been excluded.  Proteins with multiple
human taxonomy reference sequence assignments (e.g., 2 UniProt ID within human taxonomy scope) are NOT excluded.
Some artifacts in which PDB taxonomy assignment differs from the SIFTS reassignment have been excluded.

```none
    Multi-taxonomy entities: 1173 (excluded)
    Bad SIFTS assignments (switched taxonomies): 131 (excluded)
```

Clustering performed on all protein polymer entities.  The following table gives the
essential statistics for the clustering solution.

| % Sequence Identity | Distinct Groups in the Cluster|
| :-----: | :-------: |
| 30  | 36399 |
| 40  | 41347 |
| 50  | 45782 |
| 70  | 53256 |
| 90  | 61609 |
| 95  | 65572 |
| 100 | 94209 |

For the current human protein entity cohort, the populated cluster statistics
are summarized below.

|% Sequence Identity | Cluster Groups Containing a Human Protein|
| :-----: | :-------: |
| 100 | 23295 |
|  95 | 15165 |
|  90 | 13807 |

## Protein Entity-level Data Sets

Data files containing the leading human protein representative of each cluster group are attached in
`100-first-human-entity-abbrev.csv`, `95-first-human-entity-abbrev.csv`, and `90-first-human-entity-abbrev.csv`
for sequence identities 100, 95, and 90, respectively. Some cluster groups have been extensively annotated
by UniProt and these groups may contain multiple reference sequence assignments.  A separate collection of file is
provided which contains the leading representative of each cluster as well as the leading representative for
each reference sequence assignment.  These file are attached for each level of sequence identity in
`100-first-human-entity-full.csv`, `95-first-human-entity-full.csv`, and `90-first-human-entity-full.csv`.
These files are ordered by the size of the cluster group largest to smallest. To look at the leading examples chronologically, sort the data set by PDB release year column.

|Column Name | Description |
| :-------- | :---------- |
| Cluster_ID | identifer for the cluster group |
| Cluster_Members_Total  | total number of polymer entities in the cluster group |
| Cluster_Members_Human  | number of human protein members in the cluster group |
| Cluster_UniProt_Count  | numer of distinct reference sequences assigned in the cluster group |
| PDB_Entity_ID  | PDB Entity ID for the leading human protein in the cluster group.  If multiple reference sequences are assigned to the cluster group, the leading Entity ID is for each reference sequence is given in files labeled `full`. Files labeled `abbrev` contain only the leading entity in each cluster group. |
| PDB_Release_Year  | release year for the entry containing the protein entity |
| UniProt_IDs  | reference sequence assignments in the cluster group |
| Assign_Count  | count of reference sequences assigned in the cluster group |
| PDB_Struct_title  | PDB entry structure title|
| PDB_Struct_Descr  | PDB entry structure description |
| PDB_Entity_Descr  | PDB entity description |
| UniProt_Name  | UniProt recommended protein name|
| Uniprot_Gene  | UniProt primary gene name |

## Release statistics for Human Proteins

Counts of entries containing human and human hosted proteins released by year
are summarized in the following table.  These data are include
in attached file `human-containing-entries-by-year.csv`.

| Year| Entries_Containing_Human_Proteins |
| ----- | :-------: |
| 1976 | 2 |
| 1977 | 1 |
| 1979 | 2 |
| 1981 | 2 |
| 1982 | 1 |
| 1983 | 3 |
| 1984 | 4 |
| 1985 | 1 |
| 1987 | 1 |
| 1988 | 3 |
| 1989 | 6 |
| 1990 | 30 |
| 1991 | 13 |
| 1992 | 51 |
| 1993 | 123 |
| 1994 | 262 |
| 1995 | 214 |
| 1996 | 316 |
| 1997 | 325 |
| 1998 | 484 |
| 1999 | 572 |
| 2000 | 582 |
| 2001 | 608 |
| 2002 | 657 |
| 2003 | 872 |
| 2004 | 1149 |
| 2005 | 1466 |
| 2006 | 1916 |
| 2007 | 2055 |
| 2008 | 1791 |
| 2009 | 1958 |
| 2010 | 2147 |
| 2011 | 2245 |
| 2012 | 2635 |
| 2013 | 2916 |
| 2014 | 2905 |
| 2015 | 2997 |
| 2016 | 3614 |
| 2017 | 4402 |
| 2018 | 4192 |
| 2019 | 3988 |
| 2020 | 4072 |

## Release Statistics for Leading Examples of Human Proteins

Counts of entries containing the leading example of a human protein entity by year
are summarized in the following tables cluster sequence identity 100, 95 and 90 percent.
These data are included in attached files `100-leading-human-containing-entities-by-year.csv`,
`95-leading-human-containing-entities-by-year.csv`, and `90-leading-human-containing-entities-by-year.csv`.

### Sequence Identity 100%

| Year| Entries With Leading Human Protein |
| :-----: | :-------: |
| 1976 | 2 |
| 1977 | 1 |
| 1979 | 1 |
| 1981 | 1 |
| 1982 | 1 |
| 1983 | 2 |
| 1984 | 1 |
| 1985 | 1 |
| 1987 | 1 |
| 1988 | 3 |
| 1989 | 6 |
| 1990 | 12 |
| 1991 | 7 |
| 1992 | 34 |
| 1993 | 70 |
| 1994 | 106 |
| 1995 | 105 |
| 1996 | 103 |
| 1997 | 185 |
| 1998 | 270 |
| 1999 | 310 |
| 2000 | 335 |
| 2001 | 347 |
| 2002 | 314 |
| 2003 | 406 |
| 2004 | 552 |
| 2005 | 838 |
| 2006 | 1028 |
| 2007 | 1165 |
| 2008 | 853 |
| 2009 | 829 |
| 2010 | 852 |
| 2011 | 867 |
| 2012 | 924 |
| 2013 | 1008 |
| 2014 | 1087 |
| 2015 | 1015 |
| 2016 | 1232 |
| 2017 | 1188 |
| 2018 | 1236 |
| 2019 | 1348 |
| 2020 | 1206 |

### Sequence Identity 95%

| Year| Entries With Leading Human Protein |
| :-----: | :-------: |
| 1976 | 2 |
| 1977 | 1 |
| 1979 | 1 |
| 1981 | 1 |
| 1983 | 2 |
| 1984 | 1 |
| 1985 | 1 |
| 1987 | 1 |
| 1988 | 3 |
| 1989 | 5 |
| 1990 | 9 |
| 1991 | 7 |
| 1992 | 16  |
| 1993 | 39  |
| 1994 | 69  |
| 1995 | 72  |
| 1996 | 68  |
| 1997 | 122 |
| 1998 | 162 |
| 1999 | 210 |
| 2000 | 214 |
| 2001 | 212 |
| 2002 | 210 |
| 2003 | 265 |
| 2004 | 392 |
| 2005 | 701 |
| 2006 | 771 |
| 2007 | 903 |
| 2008 | 592 |
| 2009 | 557 |
| 2010 | 570 |
| 2011 | 554 |
| 2012 | 580 |
| 2013 | 609 |
| 2014 | 621 |
| 2015 | 610 |
| 2016 | 711 |
| 2017 | 750 |
| 2018 | 737 |
| 2019 | 830 |
| 2020 | 783 |

### Sequence Identity 90%

| Year| Entries_With_Leading_Human_Protein |
| :-----: | :-------: |
| 1976 | 2 |
| 1977 | 1 |
| 1979 | 1 |
| 1981 | 1 |
| 1983 | 2 |
| 1984 | 1 |
| 1985 | 1 |
| 1987 | 1 |
| 1988 | 3 |
| 1989 | 5 |
| 1990 | 9 |
| 1991 | 7 |
| 1992 | 15 |
| 1993 | 38 |
| 1994 | 68 |
| 1995 | 68 |
| 1996 | 66 |
| 1997 | 119 |
| 1998 | 154 |
| 1999 | 206 |
| 2000 | 210 |
| 2001 | 208 |
| 2002 | 208 |
| 2003 | 260 |
| 2004 | 376 |
| 2005 | 691 |
| 2006 | 744 |
| 2007 | 868 |
| 2008 | 562 |
| 2009 | 533 |
| 2010 | 539 |
| 2011 | 518 |
| 2012 | 545 |
| 2013 | 560 |
| 2014 | 553 |
| 2015 | 551 |
| 2016 | 638 |
| 2017 | 663 |
| 2018 | 641 |
| 2019 | 751 |
| 2020 | 678 |
