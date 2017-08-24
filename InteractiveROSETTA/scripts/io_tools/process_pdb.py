#!/usr/bin/env python
# :noTabs=true:

"""
create a directory of the contents in a PDB
splits into chains, grouped chains (pairings parsed from the header),
individual HETATM PDB lines, sequence files (FASTA), etc.

verbosely:
    This method behaves slightly differently for PDB files with multiple models,
    nucleic acids, duplicate complexes, etc.
    so if you are interested in the specifics, please read the source code
    
    In short, it tries to write:
        header.txt          a text file of the header lines
        numbering_map.txt   a text file showing 1-indexed PDB numbering
        clean.pdb           only ATOM lines
        hetatm.pdb          only HETATM lines, may be split by resName
        .fa                 sequences of all peptides and nucleic acids
        subdirectories      for each protein model/subunit (similar info)
    
    does not write a text file for the "trailer" (lines after the coordinates)
    
    converts lines (ATOM or HETATM) that can be converted based on  <conversion>
    (generally) and  <na_conversion>  (specific for nucleic acids, relevant
    because RNA and DNA may require different treatment...)
    !!!WARNING!!! defaults:
        CSE     CYS     converts SelenoCysteinE to Cysteine
        HYP     PRO     converts HYdroxylProline to Proline
        CYD     CYS     does NOT convert "CYsteine Disulfides to Cysteine"
        HIP     HIS     converts "HIP" to Histidine (~double protonation)
        HID     HIS     converts "HID" to Histidine (~single delta N proton)
        HIE     HIS     converts "HIE" to Histidine (~single epsilon N proton)

    todo:
        ensure hetatm conversions step, illegal atoms!!!!
        alternate conformations (mostly supported now)
        convert DNA to Rosetta DNA
        convert ligands to params
        convert water to TP3 (or TP5)



Methods for cleaning and parsing PDB files

Most importantly, the process_pdb method does a lot to clean PDB files
from RCSB

Requires:
    
    Biopython
    
Author: Evan H. Baugh
"""

################################################################################
# IMPORT

# common modules
import optparse    # for commandline
import os
import shutil

# bigger modules
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from Bio.PDB import PPBuilder    # no longer used, much faster way to do this
#from Bio.PDB import Select    # no longer used...kinda hard to use
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# custom modules
#from helper import get_root_filename , create_directory , copy_file
#from settings import SEQFORMAT , SEQFORMAT_EXTENSION_MAP , NUCLEIC_SEQUENCE_LETTERS_MAP , NA_CODES , three2one , WATER_CONVERSION , one2three , three2three , NA_CONVERSIONS_ROSETTA
#from biopython_settings import DNAAlphabet , ProteinAlphabet

#from seq_basics import write_sequence , get_sequence

################################################################################
# SETTINGS

# unholy settings...too many...

SEQFORMAT = 'fasta'
SEQFORMAT_EXTENSION_MAP = {
    'fasta' : 'fa' ,
    'genbank' : 'gb' ,
    'clustal' : 'aln' ,
    'stockholm' : 'ann'
    }

# mapping for sequence file extensions
# update when you start using something new
SEQFORMAT_MAP = {
    'fa' : 'fasta' ,
    'fas' : 'fasta' ,
    'fasta' : 'fasta' ,
    'gbk' : 'genbank' ,
    'gb' : 'genbank' ,
    'aln' : 'clustal' ,
    'ann' : 'stockholm' ,    # Pfam uses these
    'pir' : 'pir' ,    # used by Modeller...
    'sp' : 'swiss'    # uniprot/swissprot
    }

# Biopython Alphabets
DNAAlphabet = IUPAC.unambiguous_dna    # requires Biopython
ProteinAlphabet = IUPAC.protein    # requires Biopython

# simple amino acid mapping
one2three = {
    'A':'ALA',
    'C':'CYS',
    'D':'ASP',
    'E':'GLU',
    'F':'PHE',
    'G':'GLY',
    'H':'HIS',
    'I':'ILE',
    'K':'LYS',
    'L':'LEU',
    'M':'MET',
    'N':'ASN',
    'P':'PRO',
    'Q':'GLN',
    'R':'ARG',
    'S':'SER',
    'T':'THR',
    'V':'VAL',
    'W':'TRP',
    'Y':'TYR',
    }

# the revers of above...maybe more?
three2one = {
    'ALA':'A',
    'CYS':'C',
    'ASP':'D',
    'GLU':'E',
    'PHE':'F',
    'GLY':'G',
    'HIS':'H',
    'ILE':'I',
    'LYS':'K',
    'LEU':'L',
    'MET':'M',
    'ASN':'N',
    'PRO':'P',
    'GLN':'Q',
    'ARG':'R',
    'SER':'S',
    'THR':'T',
    'VAL':'V',
    'TRP':'W',
    'TYR':'Y',
    # pseudo-standard 3 letter codes for the standard aa
    'CYD' : 'C' ,
    'CYZ' : 'C' ,
    'HID' : 'H' ,
    'HIE' : 'H' ,
    'HIP' : 'H' ,
    # just to be sure...
    'ala':'A',
    'cys':'C',
    'asp':'D',
    'glu':'E',
    'phe':'F',
    'gly':'G',
    'his':'H',
    'ile':'I',
    'lys':'K',
    'leu':'L',
    'met':'M',
    'asn':'N',
    'pro':'P',
    'gln':'Q',
    'arg':'R',
    'ser':'S',
    'thr':'T',
    'val':'V',
    'trp':'W',
    'tyr':'Y',
    'Ala':'A',
    'Cys':'C',
    'Asp':'D',
    'Glu':'E',
    'Phe':'F',
    'Gly':'G',
    'His':'H',
    'Ile':'I',
    'Lys':'K',
    'Leu':'L',
    'Met':'M',
    'Asn':'N',
    'Pro':'P',
    'Gln':'Q',
    'Arg':'R',
    'Ser':'S',
    'Thr':'T',
    'Val':'V',
    'Trp':'W',
    'Tyr':'Y',
    }

###################
# HETATM CONVERSION

# unsure about these...may include ATOM or HETATM lines...

#from http://astral.stanford.edu/scopseq-1.55/release-notes-1.55.txt
three2three = {
    'AIB' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'ALA' : 'ALA' ,    # ALA
    'ALM' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'AYA' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'BNN' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'CHG' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'CSD' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'DAL' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'DHA' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'DNP' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'FLA' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'HAC' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'PRR' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'MAA' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'TIH' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    'TPQ' : 'ALA' ,    # HETEROATOM THAT MAY BE TREATED AS ALA
    '0CS':'ALA',                    ##  0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
    '2BU':'ALA',                    ##  2BU ADE
    '2OP':'ALA',                    ##  2OP (2S  2-HYDROXYPROPANAL
    '4F3':'ALA',                    ##  4F3 ALA  CYCLIZED
    'AA4':'ALA',                    ##  AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
    'ABA':'ALA',                    ##  ABA ALA  ALPHA-AMINOBUTYRIC ACID
    'AHO':'ALA',                    ##  AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
    'AHP':'ALA',                    ##  AHP ALA  2-AMINO-HEPTANOIC ACID
    'AIB':'ALA',                    ##  AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
    'ALA':'ALA',                    ##  ALA ALA
    'ALC':'ALA',                    ##  ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
    'ALM':'ALA',                    ##  ALM ALA  1-METHYL-ALANINAL
    'ALN':'ALA',                    ##  ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
    'ALS':'ALA',                    ##  ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
    'ALT':'ALA',                    ##  ALT ALA  THIOALANINE
    'AP7':'ALA',                    ##  AP7 ADE
    'APH':'ALA',                    ##  APH ALA  P-AMIDINOPHENYL-3-ALANINE
    'AYA':'ALA',                    ##  AYA ALA  N-ACETYLALANINE
    'AYG':'ALA',                    ##  AYG ALA
    'B2A':'ALA',                    ##  B2A ALA  ALANINE BORONIC ACID
    'B3A':'ALA',                    ##  B3A ALA  (3S)-3-AMINOBUTANOIC ACID
    'BAL':'ALA',                    ##  BAL ALA  BETA-ALANINE
    'BNN':'ALA',                    ##  BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
    'C12':'ALA',                    ##  C12 ALA
    'C99':'ALA',                    ##  C99 ALA
    'CAB':'ALA',                    ##  CAB ALA  4-CARBOXY-4-AMINOBUTANAL
    'CH6':'ALA',                    ##  CH6 ALA
    'CH7':'ALA',                    ##  CH7 ALA
    'CLB':'ALA',                    ##  CLB ALA
    'CLD':'ALA',                    ##  CLD ALA
    'CLV':'ALA',                    ##  CLV ALA
    'CQR':'ALA',                    ##  CQR ALA
    'CR2':'ALA',                    ##  CR2 ALA  POST-TRANSLATIONAL MODIFICATION
    'CR5':'ALA',                    ##  CR5 ALA
    'CR7':'ALA',                    ##  CR7 ALA
    'CR8':'ALA',                    ##  CR8 ALA
    'CRK':'ALA',                    ##  CRK ALA
    'CRW':'ALA',                    ##  CRW ALA
    'CRX':'ALA',                    ##  CRX ALA
    'CSI':'ALA',                    ##  CSI ALA
    'CSY':'ALA',                    ##  CSY ALA  MODIFIED TYROSINE COMPLEX
    'CWR':'ALA',                    ##  CWR ALA
    'DAB':'ALA',                    ##  DAB ALA  2,4-DIAMINOBUTYRIC ACID
    'DAL':'ALA',                    ##  DAL ALA  D-ALANINE
    'DAM':'ALA',                    ##  DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
    'DBU':'ALA',                    ##  DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
    'DBZ':'ALA',                    ##  DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
    'DHA':'ALA',                    ##  DHA ALA  2-AMINO-ACRYLIC ACID
    'DPP':'ALA',                    ##  DPP ALA  DIAMMINOPROPANOIC ACID
'FGL':'ALA',                    ##  FGL ALA  2-AMINOPROPANEDIOIC ACID
'DYG':'ALA',                    ##  DYG ALA
'GMU':'ALA',                    ##  GMU 5MU
'HHK':'ALA',                    ##  HHK ALA  (2S)-2,8-DIAMINOOCTANOIC ACID
'HMF':'ALA',                    ##  HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
'IAM':'ALA',                    ##  IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
'IGL':'ALA',                    ##  IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
'KYN':'ALA',                    ##  KYN ALA  KYNURENINE
'LAL':'ALA',                    ##  LAL ALA  N,N-DIMETHYL-L-ALANINE
'MAA':'ALA',                    ##  MAA ALA  N-METHYLALANINE
'MDO':'ALA',                    ##  MDO ALA
'MFC':'ALA',                    ##  MFC ALA  CYCLIZED
'NAL':'ALA',                    ##  NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
'NAM':'ALA',                    ##  NAM ALA  NAM NAPTHYLAMINOALANINE
'NCB':'ALA',                    ##  NCB ALA  CHEMICAL MODIFICATION
'NRQ':'ALA',                    ##  NRQ ALA
'NYC':'ALA',                    ##  NYC ALA
'ORN':'ALA',                    ##  ORN ALA  ORNITHINE
'PIA':'ALA',                    ##  PIA ALA  FUSION OF ALA 65, TYR 66, GLY 67
'PRR':'ALA',                    ##  PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
'PYA':'ALA',                    ##  PYA ALA  3-(1,10-PHENANTHROL-2-YL)-L-ALANINE
'PYC':'ALA',                    ##  PYC ALA  PYRROLE-2-CARBOXYLATE
'PYT':'ALA',                    ##  PYT ALA  MODIFIED ALANINE
'RC7':'ALA',                    ##  RC7 ALA
'SEC':'ALA',                    ##  SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
'SIC':'ALA',                    ##  SIC ALA
'SUI':'ALA',                    ##  SUI ALA
'TIH':'ALA',                    ##  TIH ALA  BETA(2-THIENYL)ALANINE
'TPQ':'ALA',                    ##  TPQ ALA  2,4,5-TRIHYDROXYPHENYLALANINE
'UMA':'ALA',                    ##  UMA ALA
'X9Q':'ALA',                    ##  X9Q ALA
'XXY':'ALA',                    ##  XXY ALA
'XYG':'ALA',                    ##  XYG ALA

#    'ASX' : 'B' ,    # why is this here!?

    'BCS' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'BUC' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'C5C' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'C6C' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CCS' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CEA' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CME' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CSO' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CSP' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CSS' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CSX' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CSW' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CY1' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CY3' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CYG' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CYM' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'CYS' : 'CYS' ,    # CYS
    'CYQ' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'DCY' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'EFC' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'OCS' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'PEC' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'PR3' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'SCH' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'SCS' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'SCY' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'SHC' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'SMC' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
    'SOC' : 'CYS' ,    # HETEROATOM THAT MAY BE TREATED AS CYS
'5CS':'CYS',                    ##  5CS CYS
'AGT':'CYS',                    ##  AGT CYS  AGMATINE-CYSTEINE ADDUCT
'BBC':'CYS',                    ##  BBC CYS
'BCS':'CYS',                    ##  BCS CYS  BENZYLCYSTEINE
'BCX':'CYS',                    ##  BCX CYS  BETA-3-CYSTEINE
'BPE':'CYS',                    ##  BPE CYS
'BUC':'CYS',                    ##  BUC CYS  S,S-BUTYLTHIOCYSTEINE
'C3Y':'CYS',                    ##  C3Y CYS  MODIFIED CYSTEINE
'C5C':'CYS',                    ##  C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
'C6C':'CYS',                    ##  C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
'CAF':'CYS',                    ##  CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
'CAS':'CYS',                    ##  CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
'CCS':'CYS',                    ##  CCS CYS  CARBOXYMETHYLATED CYSTEINE
'CME':'CYS',                    ##  CME CYS  MODIFIED CYSTEINE
'CML':'CYS',                    ##  CML CYS
'CMT':'CYS',                    ##  CMT CYS  O-METHYLCYSTEINE
'CS1':'CYS',                    ##  CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
'CS3':'CYS',                    ##  CS3 CYS
'CS4':'CYS',                    ##  CS4 CYS
'CSA':'CYS',                    ##  CSA CYS  S-ACETONYLCYSTEIN
'CSB':'CYS',                    ##  CSB CYS  CYS BOUND TO LEAD ION
'CSD':'CYS',                    ##  CSD CYS  3-SULFINOALANINE
'CSE':'CYS',                    ##  CSE CYS  SELENOCYSTEINE
'CSO':'CYS',                    ##  CSO CYS  INE S-HYDROXYCYSTEINE
'CSR':'CYS',                    ##  CSR CYS  S-ARSONOCYSTEINE
'CSS':'CYS',                    ##  CSS CYS  1,3-THIAZOLE-4-CARBOXYLIC ACID
'CSU':'CYS',                    ##  CSU CYS  CYSTEINE-S-SULFONIC ACID
'CSW':'CYS',                    ##  CSW CYS  CYSTEINE-S-DIOXIDE
'CSX':'CYS',                    ##  CSX CYS  OXOCYSTEINE
'CSZ':'CYS',                    ##  CSZ CYS  S-SELANYL CYSTEINE
'CY0':'CYS',                    ##  CY0 CYS  MODIFIED CYSTEINE
'CY1':'CYS',                    ##  CY1 CYS  ACETAMIDOMETHYLCYSTEINE
'CY3':'CYS',                    ##  CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
'CY4':'CYS',                    ##  CY4 CYS  S-BUTYRYL-CYSTEIN
'CY7':'CYS',                    ##  CY7 CYS  MODIFIED CYSTEINE
#'CYD':'CYS',                    ##  CYD CYS
'CYF':'CYS',                    ##  CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
'CYG':'CYS',                    ##  CYG CYS
'CYQ':'CYS',                    ##  CYQ CYS
'CYR':'CYS',                    ##  CYR CYS
'CYS':'CYS',                    ##  CYS CYS
'CZ2':'CYS',                    ##  CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
'CZZ':'CYS',                    ##  CZZ CYS  THIARSAHYDROXY-CYSTEINE
'DCY':'CYS',                    ##  DCY CYS  D-CYSTEINE
'DYS':'CYS',                    ##  DYS CYS
'EFC':'CYS',                    ##  EFC CYS  S,S-(2-FLUOROETHYL)THIOCYSTEINE
'FOE':'CYS',                    ##  FOE CYS
'GT9':'CYS',                    ##  GT9 CYS  SG ALKYLATED
'GYC':'CYS',                    ##  GYC CYS
'HTI':'CYS',                    ##  HTI CYS
'KOR':'CYS',                    ##  KOR CYS  MODIFIED CYSTEINE
'M0H':'CYS',                    ##  M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
'MCS':'CYS',                    ##  MCS CYS  MALONYLCYSTEINE
'NPH':'CYS',                    ##  NPH CYS
'NYS':'CYS',                    ##  NYS CYS
'OCS':'CYS',                    ##  OCS CYS  CYSTEINE SULFONIC ACID
'OCY':'CYS',                    ##  OCY CYS  HYDROXYETHYLCYSTEINE
'P1L':'CYS',                    ##  P1L CYS  S-PALMITOYL CYSTEINE
'PBB':'CYS',                    ##  PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
'PEC':'CYS',                    ##  PEC CYS  S,S-PENTYLTHIOCYSTEINE
'PR3':'CYS',                    ##  PR3 CYS  INE DTT-CYSTEINE
'PYX':'CYS',                    ##  PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
'R1A':'CYS',                    ##  R1A CYS
'R1B':'CYS',                    ##  R1B CYS
'R1F':'CYS',                    ##  R1F CYS
'R7A':'CYS',                    ##  R7A CYS
'RCY':'CYS',                    ##  RCY CYS
'SAH':'CYS',                    ##  SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
'SC2':'CYS',                    ##  SC2 CYS  N-ACETYL-L-CYSTEINE
'SCH':'CYS',                    ##  SCH CYS  S-METHYL THIOCYSTEINE GROUP
'SCS':'CYS',                    ##  SCS CYS  MODIFIED CYSTEINE
'SCY':'CYS',                    ##  SCY CYS  CETYLATED CYSTEINE
'SHC':'CYS',                    ##  SHC CYS  S-HEXYLCYSTEINE
'SMC':'CYS',                    ##  SMC CYS  POST-TRANSLATIONAL MODIFICATION
'SNC':'CYS',                    ##  SNC CYS  S-NITROSO CYSTEINE
'SOC':'CYS',                    ##  SOC CYS  DIOXYSELENOCYSTEINE
'TEE':'CYS',                    ##  TEE CYS  POST-TRANSLATIONAL MODIFICATION
'TNB':'CYS',                    ##  TNB CYS  S-(2,3,6-TRINITROPHENYL)CYSTEINE
'TYX':'CYS',                    ##  TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
'YCM':'CYS',                    ##  YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE

    '2AS' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'ASA' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'ASB' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'ASK' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'ASL' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'ASP' : 'ASP' ,    # ASP
    'ASQ' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'BHD' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'DAS' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
    'DSP' : 'ASP' ,    # HETEROATOM THAT MAY BE TREATED AS ASP
'3MD':'ASP',                    ##  3MD ASP  2S,3S-3-METHYLASPARTIC ACID
'A0A':'ASP',                    ##  A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
'ACB':'ASP',                    ##  ACB ASP  3-METHYL-ASPARTIC ACID
'AKL':'ASP',                    ##  AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
'ASA':'ASP',                    ##  ASA ASP  ASPARTIC ALDEHYDE
'ASB':'ASP',                    ##  ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
'ASI':'ASP',                    ##  ASI ASP  L-ISO-ASPARTATE
'ASK':'ASP',                    ##  ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
'ASL':'ASP',                    ##  ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
'ASP':'ASP',                    ##  ASP ASP
'B3D':'ASP',                    ##  B3D ASP  3-AMINOPENTANEDIOIC ACID
'BFD':'ASP',                    ##  BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
'BHD':'ASP',                    ##  BHD ASP  BETA-HYDROXYASPARTIC ACID
'DAS':'ASP',                    ##  DAS ASP  D-ASPARTIC ACID
'DMK':'ASP',                    ##  DMK ASP  DIMETHYL ASPARTIC ACID
'IAS':'ASP',                    ##  IAS ASP  ASPARTYL GROUP
'OHS':'ASP',                    ##  OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
'OXX':'ASP',                    ##  OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
'PHD':'ASP',                    ##  PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
'SNN':'ASP',                    ##  SNN ASP  POST-TRANSLATIONAL MODIFICATION

    '5HP' : 'GLU' ,    # HETEROATOM THAT MAY BE TREATED AS GLU
    'CGU' : 'GLU' ,    # HETEROATOM THAT MAY BE TREATED AS GLU
    'DGL' : 'GLU' ,    # HETEROATOM THAT MAY BE TREATED AS GLU
    'GGL' : 'GLU' ,    # HETEROATOM THAT MAY BE TREATED AS GLU
    'GLU' : 'GLU' ,    # GLU
    'GMA' : 'GLU' ,    # HETEROATOM THAT MAY BE TREATED AS GLU
    'PCA' : 'GLU' ,    # HETEROATOM THAT MAY BE TREATED AS GLU
'AB7':'GLU',                    ##  AB7 GLU  ALPHA-AMINOBUTYRIC ACID
'AR4':'GLU',                    ##  AR4 GLU
'B3E':'GLU',                    ##  B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
'CGU':'GLU',                    ##  CGU GLU  CARBOXYLATION OF THE CG ATOM
'DGL':'GLU',                    ##  DGL GLU  D-GLU
'GLU':'GLU',                    ##  GLU GLU
'GMA':'GLU',                    ##  GMA GLU  1-AMIDO-GLUTAMIC ACID
'ILG':'GLU',                    ##  ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
'LME':'GLU',                    ##  LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
'MEG':'GLU',                    ##  MEG GLU  (2S,3R)-3-METHYL-GLUTAMIC ACID

    'DAH' : 'PHE' ,    # HETEROATOM THAT MAY BE TREATED AS PHE
    'DPN' : 'PHE' ,    # HETEROATOM THAT MAY BE TREATED AS PHE
    'HPQ' : 'PHE' ,    # HETEROATOM THAT MAY BE TREATED AS PHE
    'PHE' : 'PHE' ,    # PHE
    'PHI' : 'PHE' ,    # HETEROATOM THAT MAY BE TREATED AS PHE
    'PHL' : 'PHE' ,    # HETEROATOM THAT MAY BE TREATED AS PHE
'1PA':'PHE',                    ##  1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
'23F':'PHE',                    ##  23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
'4PH':'PHE',                    ##  4PH PHE  4-METHYL-L-PHENYLALANINE
'B2F':'PHE',                    ##  B2F PHE  PHENYLALANINE BORONIC ACID
'BIF':'PHE',                    ##  BIF PHE
'CHS':'PHE',                    ##  CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
'DAH':'PHE',                    ##  DAH PHE  3,4-DIHYDROXYDAHNYLALANINE
'DPH':'PHE',                    ##  DPH PHE  DEAMINO-METHYL-PHENYLALANINE
'DPN':'PHE',                    ##  DPN PHE  D-CONFIGURATION
'FCL':'PHE',                    ##  FCL PHE  3-CHLORO-L-PHENYLALANINE
'FOG':'PHE',                    ##  FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
'FRF':'PHE',                    ##  FRF PHE  PHE FOLLOWED BY REDUCED PHE
'HPE':'PHE',                    ##  HPE PHE  HOMOPHENYLALANINE
'HPH':'PHE',                    ##  HPH PHE  PHENYLALANINOL GROUP
'HPQ':'PHE',                    ##  HPQ PHE  HOMOPHENYLALANINYLMETHANE
'MEA':'PHE',                    ##  MEA PHE  N-METHYLPHENYLALANINE
'MTY':'PHE',                    ##  MTY PHE  3-HYDROXYPHENYLALANINE
'NFA':'PHE',                    ##  NFA PHE  MODIFIED PHENYLALANINE
'PBF':'PHE',                    ##  PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
'PCS':'PHE',                    ##  PCS PHE  PHENYLALANYLMETHYLCHLORIDE
'PF5':'PHE',                    ##  PF5 PHE  2,3,4,5,6-PENTAFLUORO-L-PHENYLALANINE
'PFF':'PHE',                    ##  PFF PHE  4-FLUORO-L-PHENYLALANINE
'PHA':'PHE',                    ##  PHA PHE  PHENYLALANINAL
'PHE':'PHE',                    ##  PHE PHE
'PHI':'PHE',                    ##  PHI PHE  IODO-PHENYLALANINE
'PHL':'PHE',                    ##  PHL PHE  L-PHENYLALANINOL
'PHM':'PHE',                    ##  PHM PHE  PHENYLALANYLMETHANE
'PM3':'PHE',                    ##  PM3 PHE
'PPN':'PHE',                    ##  PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
'PRQ':'PHE',                    ##  PRQ PHE  PHENYLALANINE
'PSA':'PHE',                    ##  PSA PHE
'SMF':'PHE',                    ##  SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE

    'GL3' : 'GLY' ,    # HETEROATOM THAT MAY BE TREATED AS GLY
    'GLY' : 'GLY' ,    # GLY
    'GLZ' : 'GLY' ,    # HETEROATOM THAT MAY BE TREATED AS GLY
    'GSC' : 'GLY' ,    # HETEROATOM THAT MAY BE TREATED AS GLY
    'MPQ' : 'GLY' ,    # HETEROATOM THAT MAY BE TREATED AS GLY
    'MSA' : 'GLY' ,    # HETEROATOM THAT MAY BE TREATED AS GLY
    'NMC' : 'GLY' ,    # HETEROATOM THAT MAY BE TREATED AS GLY
    'SAR' : 'GLY' ,    # HETEROATOM THAT MAY BE TREATED AS GLY
'ACY':'GLY',                    ##  ACY GLY  POST-TRANSLATIONAL MODIFICATION
'CHG':'GLY',                    ##  CHG GLY  CYCLOHEXYL GLYCINE
'CHP':'GLY',                    ##  CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
'GHP':'GLY',                    ##  GHP GLY  4-HYDROXYPHENYLGLYCINE
'GL3':'GLY',                    ##  GL3 GLY  POST-TRANSLATIONAL MODIFICATION
'GLY':'GLY',                    ##  GLY GLY
'GLZ':'GLY',                    ##  GLZ GLY  AMINO-ACETALDEHYDE
'GYS':'GLY',                    ##  GYS GLY
'IPG':'GLY',                    ##  IPG GLY  N-ISOPROPYL GLYCINE
'MEU':'GLY',                    ##  MEU GLY  O-METHYL-GLYCINE
'MPQ':'GLY',                    ##  MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
'MSA':'GLY',                    ##  MSA GLY  (2-S-METHYL) SARCOSINE
'NMC':'GLY',                    ##  NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
'PG9':'GLY',                    ##  PG9 GLY  D-PHENYLGLYCINE
'SAR':'GLY',                    ##  SAR GLY  SARCOSINE
'SHP':'GLY',                    ##  SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
'TBG':'GLY',                    ##  TBG GLY  T-BUTYL GLYCINE

    '3AH' : 'HIS' ,    # HETEROATOM THAT MAY BE TREATED AS HIS
    'DHI' : 'HIS' ,    # HETEROATOM THAT MAY BE TREATED AS HIS
    'HIC' : 'HIS' ,    # HETEROATOM THAT MAY BE TREATED AS HIS
    'HIS' : 'HIS' ,    # HIS
    'MHS' : 'HIS' ,    # HETEROATOM THAT MAY BE TREATED AS HIS
    'NEM' : 'HIS' ,    # HETEROATOM THAT MAY BE TREATED AS HIS
    'NEP' : 'HIS' ,    # HETEROATOM THAT MAY BE TREATED AS HIS
    'HID' : 'HIS' ,    # single delta N protonation
    'HIE' : 'HIS' ,    # single epsilon N protonation
'3AH':'HIS',                    ##  3AH HIS
'DDE':'HIS',                    ##  DDE HIS
'DHI':'HIS',                    ##  DHI HIS  D-HISTIDINE
'HIA':'HIS',                    ##  HIA HIS  L-HISTIDINE AMIDE
'HIC':'HIS',                    ##  HIC HIS  4-METHYL-HISTIDINE
'HIP':'HIS',                    ##  HIP HIS  ND1-PHOSPHONOHISTIDINE...or commonly used doubly protonated state
'HIQ':'HIS',                    ##  HIQ HIS  MODIFIED HISTIDINE
'HIS':'HIS',                    ##  HIS HIS
'HSO':'HIS',                    ##  HSO HIS  HISTIDINOL
'MHS':'HIS',                    ##  MHS HIS  1-N-METHYLHISTIDINE
'NEP':'HIS',                    ##  NEP HIS  N1-PHOSPHONOHISTIDINE
'NZH':'HIS',                    ##  NZH HIS
'OHI':'HIS',                    ##  OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
'PSH':'HIS',                    ##  PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE

    'DIL' : 'ILE' ,    # HETEROATOM THAT MAY BE TREATED AS ILE
    'IIL' : 'ILE' ,    # HETEROATOM THAT MAY BE TREATED AS ILE
    'ILE' : 'ILE' ,    # ILE
'B2I':'ILE',                    ##  B2I ILE  ISOLEUCINE BORONIC ACID
'DIL':'ILE',                    ##  DIL ILE  D-ISOLEUCINE
'IIL':'ILE',                    ##  IIL ILE  ISO-ISOLEUCINE
'ILE':'ILE',                    ##  ILE ILE
'ILX':'ILE',                    ##  ILX ILE  4,5-DIHYDROXYISOLEUCINE
'IML':'ILE',                    ##  IML ILE  N-METHYLATED

    'ALY' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'DLY' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'KCX' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'LLP' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'LLY' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'LYM' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'LYS' : 'LYS' ,    # LYS
    'LYZ' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'MLY' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'SHR' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
    'TRG' : 'LYS' ,    # HETEROATOM THAT MAY BE TREATED AS LYS
'6CL':'LYS',                    ##  6CL LYS  6-CARBOXYLYSINE
'ALY':'LYS',                    ##  ALY LYS  N(6)-ACETYLLYSINE
'API':'LYS',                    ##  API LYS  2,6-DIAMINOPIMELIC ACID
'APK':'LYS',                    ##  APK LYS
'AZK':'LYS',                    ##  AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
'B3K':'LYS',                    ##  B3K LYS  (3S)-3,7-DIAMINOHEPTANOIC ACID
'BLY':'LYS',                    ##  BLY LYS  LYSINE BORONIC ACID
'C1X':'LYS',                    ##  C1X LYS  MODIFIED LYSINE
'CLG':'LYS',                    ##  CLG LYS
'CLH':'LYS',                    ##  CLH LYS
'CYJ':'LYS',                    ##  CYJ LYS  MODIFIED LYSINE
'DLS':'LYS',                    ##  DLS LYS  DI-ACETYL-LYSINE
'DLY':'LYS',                    ##  DLY LYS  D-LYSINE
'DNL':'LYS',                    ##  DNL LYS  6-AMINO-HEXANAL
'FHL':'LYS',                    ##  FHL LYS  MODIFIED LYSINE
'GPL':'LYS',                    ##  GPL LYS  LYSINE GUANOSINE-5'-MONOPHOSPHATE
'IT1':'LYS',                    ##  IT1 LYS
'KCX':'LYS',                    ##  KCX LYS  CARBAMOYLATED LYSINE
'KGC':'LYS',                    ##  KGC LYS
'KST':'LYS',                    ##  KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
'LA2':'LYS',                    ##  LA2 LYS
'LCK':'LYS',                    ##  LCK LYS
'LCX':'LYS',                    ##  LCX LYS  CARBAMYLATED LYSINE
'LDH':'LYS',                    ##  LDH LYS  N~6~-ETHYL-L-LYSINE
'LET':'LYS',                    ##  LET LYS  ODIFIED LYSINE
'LLP':'LYS',                    ##  LLP LYS
'LLY':'LYS',                    ##  LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
'LSO':'LYS',                    ##  LSO LYS  MODIFIED LYSINE
'LYM':'LYS',                    ##  LYM LYS  DEOXY-METHYL-LYSINE
'LYN':'LYS',                    ##  LYN LYS  2,6-DIAMINO-HEXANOIC ACID AMIDE
'LYP':'LYS',                    ##  LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
'LYR':'LYS',                    ##  LYR LYS  MODIFIED LYSINE
'LYS':'LYS',                    ##  LYS LYS
'LYX':'LYS',                    ##  LYX LYS  N''-(2-COENZYME A)-PROPANOYL-LYSINE
'LYZ':'LYS',                    ##  LYZ LYS  5-HYDROXYLYSINE
'M2L':'LYS',                    ##  M2L LYS
'M3L':'LYS',                    ##  M3L LYS  N-TRIMETHYLLYSINE
'MCL':'LYS',                    ##  MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
'MLY':'LYS',                    ##  MLY LYS  METHYLATED LYSINE
'MLZ':'LYS',                    ##  MLZ LYS  N-METHYL-LYSINE
'OBS':'LYS',                    ##  OBS LYS  MODIFIED LYSINE
'SLZ':'LYS',                    ##  SLZ LYS  L-THIALYSINE
'XX1':'LYS',                    ##  XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE

    'BUG' : 'LEU' ,    # HETEROATOM THAT MAY BE TREATED AS LEU
    'CLE' : 'LEU' ,    # HETEROATOM THAT MAY BE TREATED AS LEU
    'DLE' : 'LEU' ,    # HETEROATOM THAT MAY BE TREATED AS LEU
    'LEU' : 'LEU' ,    # LEU
    'MLE' : 'LEU' ,    # HETEROATOM THAT MAY BE TREATED AS LEU
    'NLE' : 'LEU' ,    # HETEROATOM THAT MAY BE TREATED AS LEU
    'NLN' : 'LEU' ,    # HETEROATOM THAT MAY BE TREATED AS LEU
    'NLP' : 'LEU' ,    # HETEROATOM THAT MAY BE TREATED AS LEU
'1LU':'LEU',                    ##  1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
'2ML':'LEU',                    ##  2ML LEU  2-METHYLLEUCINE
'BLE':'LEU',                    ##  BLE LEU  LEUCINE BORONIC ACID
'BUG':'LEU',                    ##  BUG LEU  TERT-LEUCYL AMINE
'CLE':'LEU',                    ##  CLE LEU  LEUCINE AMIDE
'DCL':'LEU',                    ##  DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
'DLE':'LEU',                    ##  DLE LEU  D-LEUCINE
'DNE':'LEU',                    ##  DNE LEU  D-NORLEUCINE
'DNG':'LEU',                    ##  DNG LEU  N-FORMYL-D-NORLEUCINE
'DNM':'LEU',                    ##  DNM LEU  D-N-METHYL NORLEUCINE
'FLE':'LEU',                    ##  FLE LEU  FUROYL-LEUCINE
'HLU':'LEU',                    ##  HLU LEU  BETA-HYDROXYLEUCINE
'LED':'LEU',                    ##  LED LEU  POST-TRANSLATIONAL MODIFICATION
'LEF':'LEU',                    ##  LEF LEU  2-5-FLUOROLEUCINE
'LEU':'LEU',                    ##  LEU LEU
'LNT':'LEU',                    ##  LNT LEU
'MHL':'LEU',                    ##  MHL LEU  N-METHYLATED, HYDROXY
'MLE':'LEU',                    ##  MLE LEU  N-METHYLATED
'MLL':'LEU',                    ##  MLL LEU  METHYL L-LEUCINATE
'MNL':'LEU',                    ##  MNL LEU  4,N-DIMETHYLNORLEUCINE
'NLE':'LEU',                    ##  NLE LEU  NORLEUCINE
'NLN':'LEU',                    ##  NLN LEU  NORLEUCINE AMIDE
'NLO':'LEU',                    ##  NLO LEU  O-METHYL-L-NORLEUCINE
'PLE':'LEU',                    ##  PLE LEU  LEUCINE PHOSPHINIC ACID
'PPH':'LEU',                    ##  PPH LEU  PHENYLALANINE PHOSPHINIC ACID

    'CXM' : 'MET' ,    # HETEROATOM THAT MAY BE TREATED AS MET
    'FME' : 'MET' ,    # HETEROATOM THAT MAY BE TREATED AS MET
    'MET' : 'MET' ,    # MET
    'MSE' : 'MET' ,    # HETEROATOM THAT MAY BE TREATED AS MET
    'OMT' : 'MET' ,    # HETEROATOM THAT MAY BE TREATED AS MET
'AME':'MET',                    ##  AME MET  ACETYLATED METHIONINE
'CXM':'MET',                    ##  CXM MET  N-CARBOXYMETHIONINE
'ESC':'MET',                    ##  ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
'FME':'MET',                    ##  FME MET  FORMYL-METHIONINE
'FOR':'MET',                    ##  FOR MET
'MET':'MET',                    ##  MET MET
'MHO':'MET',                    ##  MHO MET  POST-TRANSLATIONAL MODIFICATION
'MME':'MET',                    ##  MME MET  N-METHYL METHIONINE
'MSE':'MET',                    ##  MSE MET  ELENOMETHIONINE
'MSO':'MET',                    ##  MSO MET  METHIONINE SULFOXIDE
'OMT':'MET',                    ##  OMT MET  METHIONINE SULFONE
'SME':'MET',                    ##  SME MET  METHIONINE SULFOXIDE

    'ASN' : 'ASN' ,    # ASN
    'MEN' : 'ASN' ,    # HETEROATOM THAT MAY BE TREATED AS ASN
'AFA':'ASN',                    ##  AFA ASN  N-[7-METHYL-OCT-2,4-DIENOYL]ASPARAGINE
'AHB':'ASN',                    ##  AHB ASN  BETA-HYDROXYASPARAGINE
'ASN':'ASN',                    ##  ASN ASN
'B3X':'ASN',                    ##  B3X ASN  (3S)-3,5-DIAMINO-5-OXOPENTANOIC ACID
'DMH':'ASN',                    ##  DMH ASN  N4,N4-DIMETHYL-ASPARAGINE
'DSG':'ASN',                    ##  DSG ASN  D-ASPARAGINE
'MEN':'ASN',                    ##  MEN ASN  GAMMA METHYL ASPARAGINE

    'DPR' : 'PRO' ,    # HETEROATOM THAT MAY BE TREATED AS PRO
    'PRO' : 'PRO' ,    # PRO
'1AB':'PRO',                    ##  1AB PRO  1,4-DIDEOXY-1,4-IMINO-D-ARABINITOL
'2MT':'PRO',                    ##  2MT PRO
'4FB':'PRO',                    ##  4FB PRO  (4S)-4-FLUORO-L-PROLINE
'DPL':'PRO',                    ##  DPL PRO  4-OXOPROLINE
'DPR':'PRO',                    ##  DPR PRO  D-PROLINE
'H5M':'PRO',                    ##  H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
'HY3':'PRO',                    ##  HY3 PRO  3-HYDROXYPROLINE
'HYP':'PRO',                    ##  HYP PRO  4-HYDROXYPROLINE
'LPD':'PRO',                    ##  LPD PRO  L-PROLINAMIDE
'P2Y':'PRO',                    ##  P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
'PCA':'PRO',                    ##  PCA PRO  5-OXOPROLINE
'POM':'PRO',                    ##  POM PRO  CIS-5-METHYL-4-OXOPROLINE
'PRO':'PRO',                    ##  PRO PRO
'PRS':'PRO',                    ##  PRS PRO  THIOPROLINE

    'DGN' : 'GLN' ,    # HETEROATOM THAT MAY BE TREATED AS GLN
    'GLN' : 'GLN' ,    # GLN
'DGN':'GLN',                    ##  DGN GLN  D-GLUTAMINE
'GHG':'GLN',                    ##  GHG GLN  GAMMA-HYDROXY-GLUTAMINE
'GLH':'GLN',                    ##  GLH GLN
'GLN':'GLN',                    ##  GLN GLN
'MGN':'GLN',                    ##  MGN GLN  2-METHYL-GLUTAMINE

    'ACL' : 'ARG' ,    # HETEROATOM THAT MAY BE TREATED AS ARG
    'AGM' : 'ARG' ,    # HETEROATOM THAT MAY BE TREATED AS ARG
    'ARG' : 'ARG' ,    # ARG
    'ARM' : 'ARG' ,    # HETEROATOM THAT MAY BE TREATED AS ARG
    'DAR' : 'ARG' ,    # HETEROATOM THAT MAY BE TREATED AS ARG
    'HAR' : 'ARG' ,    # HETEROATOM THAT MAY BE TREATED AS ARG
    'HMR' : 'ARG' ,    # HETEROATOM THAT MAY BE TREATED AS ARG
'2MR':'ARG',                    ##  2MR ARG  N3, N4-DIMETHYLARGININE
'AAR':'ARG',                    ##  AAR ARG  ARGININEAMIDE
'ACL':'ARG',                    ##  ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
'AGM':'ARG',                    ##  AGM ARG  4-METHYL-ARGININE
'ALG':'ARG',                    ##  ALG ARG  GUANIDINOBUTYRYL GROUP
'AR2':'ARG',                    ##  AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
'ARG':'ARG',                    ##  ARG ARG
'ARM':'ARG',                    ##  ARM ARG  DEOXY-METHYL-ARGININE
'ARO':'ARG',                    ##  ARO ARG  C-GAMMA-HYDROXY ARGININE
'BOR':'ARG',                    ##  BOR ARG
'CIR':'ARG',                    ##  CIR ARG  CITRULLINE
'DA2':'ARG',                    ##  DA2 ARG  MODIFIED ARGININE
'DAR':'ARG',                    ##  DAR ARG  D-ARGININE
'HMR':'ARG',                    ##  HMR ARG  BETA-HOMOARGININE
'HRG':'ARG',                    ##  HRG ARG  L-HOMOARGININE
'MAI':'ARG',                    ##  MAI ARG  DEOXO-METHYLARGININE
'MGG':'ARG',                    ##  MGG ARG  MODIFIED D-ARGININE
'NMM':'ARG',                    ##  NMM ARG  MODIFIED ARGININE
'OPR':'ARG',                    ##  OPR ARG  C-(3-OXOPROPYL)ARGININE
'ORQ':'ARG',                    ##  ORQ ARG  N~5~-ACETYL-L-ORNITHINE
'TYZ':'ARG',                    ##  TYZ ARG  PARA ACETAMIDO BENZOIC ACID

    'DSN' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
    'MIS' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
    'OAS' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
    'SAC' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
    'SEL' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
    'SEP' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
    'SER' : 'SER' ,    # SER
    'SET' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
    'SVA' : 'SER' ,    # HETEROATOM THAT MAY BE TREATED AS SER
'B3S':'SER',                    ##  B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
'BG1':'SER',                    ##  BG1 SER
'DHL':'SER',                    ##  DHL SER  POST-TRANSLATIONAL MODIFICATION
'DSE':'SER',                    ##  DSE SER  D-SERINE N-METHYLATED
'DSN':'SER',                    ##  DSN SER  D-SERINE
'FGP':'SER',                    ##  FGP SER
'GVL':'SER',                    ##  GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
'HSE':'SER',                    ##  HSE SER  L-HOMOSERINE
'HSL':'SER',                    ##  HSL SER  HOMOSERINE LACTONE
'MC1':'SER',                    ##  MC1 SER  METHICILLIN ACYL-SERINE
'MIS':'SER',                    ##  MIS SER  MODIFIED SERINE
'N10':'SER',                    ##  N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
'NC1':'SER',                    ##  NC1 SER  NITROCEFIN ACYL-SERINE
'OAS':'SER',                    ##  OAS SER  O-ACETYLSERINE
'OSE':'SER',                    ##  OSE SER  O-SULFO-L-SERINE
'PG1':'SER',                    ##  PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
'PYR':'SER',                    ##  PYR SER  CHEMICALLY MODIFIED
'S1H':'SER',                    ##  S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
'SAC':'SER',                    ##  SAC SER  N-ACETYL-SERINE
'SBD':'SER',                    ##  SBD SER
'SBG':'SER',                    ##  SBG SER  MODIFIED SERINE
'SBL':'SER',                    ##  SBL SER
'SDP':'SER',                    ##  SDP SER
'SEB':'SER',                    ##  SEB SER  O-BENZYLSULFONYL-SERINE
'SEL':'SER',                    ##  SEL SER  2-AMINO-1,3-PROPANEDIOL
'SEP':'SER',                    ##  SEP SER  E PHOSPHOSERINE
'SER':'SER',                    ##  SER SER
'SET':'SER',                    ##  SET SER  AMINOSERINE
'SGB':'SER',                    ##  SGB SER  MODIFIED SERINE
'SGR':'SER',                    ##  SGR SER  MODIFIED SERINE
'SOY':'SER',                    ##  SOY SER  OXACILLOYL-ACYLATED SERINE
'SUN':'SER',                    ##  SUN SER  TABUN CONJUGATED SERINE
'SVA':'SER',                    ##  SVA SER  SERINE VANADATE
'SVV':'SER',                    ##  SVV SER  MODIFIED SERINE
'SVX':'SER',                    ##  SVX SER  MODIFIED SERINE
'SVY':'SER',                    ##  SVY SER  MODIFIED SERINE
'SVZ':'SER',                    ##  SVZ SER  MODIFIED SERINE
'SXE':'SER',                    ##  SXE SER  MODIFIED SERINE

    'ALO' : 'THR' ,    # HETEROATOM THAT MAY BE TREATED AS THR
    'BMT' : 'THR' ,    # HETEROATOM THAT MAY BE TREATED AS THR
    'DTH' : 'THR' ,    # HETEROATOM THAT MAY BE TREATED AS THR
    'THR' : 'THR' ,    # THR
    'TPO' : 'THR' ,    # HETEROATOM THAT MAY BE TREATED AS THR
'AEI':'THR',                    ##  AEI THR  ACYLATED THR
'ALO':'THR',                    ##  ALO THR  ALLO-THREONINE
'BMT':'THR',                    ##  BMT THR
'CRO':'THR',                    ##  CRO THR  CYCLIZED
'CTH':'THR',                    ##  CTH THR  4-CHLOROTHREONINE
'DTH':'THR',                    ##  DTH THR  D-THREONINE
'OLT':'THR',                    ##  OLT THR  O-METHYL-L-THREONINE
'TBM':'THR',                    ##  TBM THR
'TH5':'THR',                    ##  TH5 THR  O-ACETYL-L-THREONINE
'THC':'THR',                    ##  THC THR  N-METHYLCARBONYLTHREONINE
'THR':'THR',                    ##  THR THR
'TMD':'THR',                    ##  TMD THR  N-METHYLATED, EPSILON C ALKYLATED
'TPO':'THR',                    ##  TPO THR  HOSPHOTHREONINE

    'DIV' : 'VAL' ,    # HETEROATOM THAT MAY BE TREATED AS VAL
    'DVA' : 'VAL' ,    # HETEROATOM THAT MAY BE TREATED AS VAL
    'MVA' : 'VAL' ,    # HETEROATOM THAT MAY BE TREATED AS VAL
    'VAL' : 'VAL' ,    # VAL
'B2V':'VAL',                    ##  B2V VAL  VALINE BORONIC ACID
'DIV':'VAL',                    ##  DIV VAL  D-ISOVALINE
'DVA':'VAL',                    ##  DVA VAL  D-VALINE
'MNV':'VAL',                    ##  MNV VAL  N-METHYL-C-AMINO VALINE
'MVA':'VAL',                    ##  MVA VAL  N-METHYLATED
'NVA':'VAL',                    ##  NVA VAL  NORVALINE
'VAD':'VAL',                    ##  VAD VAL  DEAMINOHYDROXYVALINE
'VAF':'VAL',                    ##  VAF VAL  METHYLVALINE
'VAL':'VAL',                    ##  VAL VAL
'VDL':'VAL',                    ##  VDL VAL  (2R,3R)-2,3-DIAMINOBUTANOIC ACID
'VLL':'VAL',                    ##  VLL VAL  (2S)-2,3-DIAMINOBUTANOIC ACID
'VME':'VAL',                    ##  VME VAL  O- METHYLVALINE

    'DTR' : 'TRP' ,    # HETEROATOM THAT MAY BE TREATED AS TRP
    'HTR' : 'TRP' ,    # HETEROATOM THAT MAY BE TREATED AS TRP
    'LTR' : 'TRP' ,    # HETEROATOM THAT MAY BE TREATED AS TRP
    'TPL' : 'TRP' ,    # HETEROATOM THAT MAY BE TREATED AS TRP
    'TRO' : 'TRP' ,    # HETEROATOM THAT MAY BE TREATED AS TRP
    'TRP' : 'TRP' ,    # TRP
'BTR':'TRP',                    ##  BTR TRP  6-BROMO-TRYPTOPHAN
'1TQ':'TRP',                    ##  1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
'23S':'TRP',                    ##  23S TRP  MODIFIED TRYPTOPHAN
'32S':'TRP',                    ##  32S TRP  MODIFIED TRYPTOPHAN
'32T':'TRP',                    ##  32T TRP  MODIFIED TRYPTOPHAN
'4DP':'TRP',                    ##  4DP TRP
'4FW':'TRP',                    ##  4FW TRP  4-FLUOROTRYPTOPHANE
'4HT':'TRP',                    ##  4HT TRP  4-HYDROXYTRYPTOPHAN
'4IN':'TRP',                    ##  4IN TRP  4-AMINO-L-TRYPTOPHAN
'6CW':'TRP',                    ##  6CW TRP  6-CHLORO-L-TRYPTOPHAN
'DTR':'TRP',                    ##  DTR TRP  D-TRYPTOPHAN
'FTR':'TRP',                    ##  FTR TRP  FLUOROTRYPTOPHANE
'HTR':'TRP',                    ##  HTR TRP  BETA-HYDROXYTRYPTOPHANE
'PAT':'TRP',                    ##  PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
'TOX':'TRP',                    ##  TOX TRP
'TPL':'TRP',                    ##  TPL TRP  TRYTOPHANOL
'TQQ':'TRP',                    ##  TQQ TRP
'TRF':'TRP',                    ##  TRF TRP  N1-FORMYL-TRYPTOPHAN
'TRN':'TRP',                    ##  TRN TRP  AZA-TRYPTOPHAN
'TRO':'TRP',                    ##  TRO TRP  2-HYDROXY-TRYPTOPHAN
'TRP':'TRP',                    ##  TRP TRP
'TRQ':'TRP',                    ##  TRQ TRP
'TRW':'TRP',                    ##  TRW TRP
'TRX':'TRP',                    ##  TRX TRP  6-HYDROXYTRYPTOPHAN
'TTQ':'TRP',                    ##  TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN

    'DTY' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'IYR' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'PAQ' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'PTR' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'STY' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'TYB' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'TYQ' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'TYR' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
    'TYS' : 'TYR' ,    # TYR
    'TYY' : 'TYR' ,    # HETEROATOM THAT MAY BE TREATED AS TYR
'1TY':'TYR',                    ##  1TY TYR
'2TY':'TYR',                    ##  2TY TYR
'3TY':'TYR',                    ##  3TY TYR  MODIFIED TYROSINE
'B3Y':'TYR',                    ##  B3Y TYR
'CRQ':'TYR',                    ##  CRQ TYR
'DBY':'TYR',                    ##  DBY TYR  3,5 DIBROMOTYROSINE
'DPQ':'TYR',                    ##  DPQ TYR  TYROSINE DERIVATIVE
'DTY':'TYR',                    ##  DTY TYR  D-TYROSINE
'ESB':'TYR',                    ##  ESB TYR
'FLT':'TYR',                    ##  FLT TYR  FLUOROMALONYL TYROSINE
'FTY':'TYR',                    ##  FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
'IYR':'TYR',                    ##  IYR TYR  3-IODO-TYROSINE
'MBQ':'TYR',                    ##  MBQ TYR
'NIY':'TYR',                    ##  NIY TYR  META-NITRO-TYROSINE
'NBQ':'TYR',                    ##  NBQ TYR
'OTY':'TYR',                    ##  OTY TYR
'PAQ':'TYR',                    ##  PAQ TYR  SEE REMARK 999
'PTH':'TYR',                    ##  PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
'PTM':'TYR',                    ##  PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
'PTR':'TYR',                    ##  PTR TYR  O-PHOSPHOTYROSINE
'TCQ':'TYR',                    ##  TCQ TYR  MODIFIED TYROSINE
'TTS':'TYR',                    ##  TTS TYR
'TY2':'TYR',                    ##  TY2 TYR  3-AMINO-L-TYROSINE
'TY3':'TYR',                    ##  TY3 TYR  3-HYDROXY-L-TYROSINE
'TYB':'TYR',                    ##  TYB TYR  TYROSINAL
'TYC':'TYR',                    ##  TYC TYR  L-TYROSINAMIDE
'TYI':'TYR',                    ##  TYI TYR  3,5-DIIODOTYROSINE
'TYN':'TYR',                    ##  TYN TYR  ADDUCT AT HYDROXY GROUP
'TYO':'TYR',                    ##  TYO TYR
'TYQ':'TYR',                    ##  TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
'TYR':'TYR',                    ##  TYR TYR
'TYS':'TYR',                    ##  TYS TYR  INE SULPHONATED TYROSINE
'TYT':'TYR',                    ##  TYT TYR
'TYY':'TYR',                    ##  TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
'YOF':'TYR',                    ##  YOF TYR  3-FLUOROTYROSINE

#    'GLX' : 'Z'    # why is this here!?
    }

####################
# NUCLEIC ACID STUFF

# for sequences...
NUCLEIC_SEQUENCE_LETTERS_MAP = {
    'A' : 'A' ,
    'G' : 'G' ,
    'C' : 'C' ,
    'T' : 'T' ,
    'U' : 'U' ,
    'a' : 'A' ,
    'g' : 'G' ,
    'c' : 'C' ,
    't' : 'T' ,
    'u' : 'U' ,
    'DA' : 'A' ,
    'DG' : 'G' ,
    'DC' : 'C' ,
    'DT' : 'T' ,
    'dA' : 'A' ,
    'dG' : 'G' ,
    'dC' : 'C' ,
    'dT' : 'T' ,
    'ADE' : 'A' ,
    'GUA' : 'G' ,
    'CYT' : 'C' ,
    'THY' : 'T' ,
    'URA' : 'U' ,
    'rA' : 'A' ,
    'rG' : 'G',
    'rC' : 'C' ,
    'rU' : 'U' ,
    # HETATM lines
    '1MA' : 'A' ,
    '1MG' : 'G' ,
    '2MG' : 'G' ,
    '7MG' : 'G' ,
    'OMG' : 'G' ,
    'YG' : 'G' ,
    '5MC' : 'C' ,
    'CB2' : 'C' ,
    'CBR' : 'C' ,
    'DC' : 'C' ,
    'OMC' : 'C' ,
    '5BU' : 'U' ,
    '5MU' : 'U' ,
    'H2U' : 'U' ,
    'PSU' : 'U' ,
    'URI' : 'U'
    }

#                line_edit = line_edit.replace( 'HO2\'', '2HO*' )
#                line_edit = line_edit.replace( 'HO5\'', '5HO*' )
#                line_edit = line_edit.replace( 'H5\'\'', '2H5*' )

#                line_edit = line_edit.replace('\'','*')
#                line_edit = line_edit.replace('OP1','O1P')
#                line_edit = line_edit.replace('OP2','O2P')

NA_CODES = {}
NA_CONVERSIONS_ROSETTA = {}

#####
# DNA

# codes whose presence indicates DNA definitively
NA_CODES['DNA'] = {
    'T' : 'T' ,
    't' : 'T' ,
    'DA' : 'A' ,
    'DG' : 'G' ,
    'DC' : 'C' ,
    'DT' : 'T' ,
    'dA' : 'A' ,
    'dG' : 'G' ,
    'dC' : 'C' ,
    'dT' : 'T' ,
    'THY' : 'T'
    }

# convert from sequence to the resName for PDB format
NA_CONVERSIONS_ROSETTA['DNA'] = {
    'A' : 'A' ,
    'G' : 'G' ,
    'C' : 'C' ,
    'T' : 'T' ,
    'ADE' : 'A' ,
    'GUA' : 'G' ,
    'CYT' : 'C' ,
    'THY' : 'T' ,
    '1MA' : 'A' ,
    '1MG' : 'G' ,
    '2MG' : 'G' ,
    '7MG' : 'G' ,
    'OMG' : 'G' ,
    'YG' : 'G' ,
    '5MC' : 'C' ,
    'CB2' : 'C' ,
    'CBR' : 'C' ,
    'DC' : 'C' ,
    'OMC' : 'C' ,
    }
    
# water! hooray!
WATER_CONVERSION = {
    'W' : 'TP3' ,
    'HOH' : 'TP3' ,
    'H2O' : 'TP3' ,
    'WAT' : 'TP3' ,
    'TP3' : 'TP3' ,
    'TP5' : 'TP3'
    }

# fun with water
#WATER_CODE = 'TP3'    # for possible use in PyRosetta
#WATER_CODES = ['W' , 'HOH' , 'H2O' , 'WAT' , 'TP3' , 'TP5']    # resNames

################################################################################
# METHODS

get_file_extension = lambda in_filename: in_filename.split( '.' )[-1]
get_file_extension.__doc__ = 'Returns the file extension of  <in_filename>\n\nin_filename.split( \'.\' )[-1]'

# hacky version
get_root_filename = lambda in_filename: in_filename[:-len( get_file_extension( in_filename ) ) - 1]
get_root_filename.__doc__ = 'Returns the \"root filename\" of  <in_filename>  (pre file extension)\n\nin_filename[:len( in_filename.split( \'.\' )[-1] ) - 1]\na little hacky...'
# better version
#get_root_filename = lambda in_filename: ''.join( [i for i in in_filename.split( '.' )[:-1]] )

# helper for creating a directory, checks and delets existing name
def create_directory( dir_name , tagline = ' to sort the data' ):
    """
    Creates the directory  <dir_name>
    
    WARNING: this will delete the directory and its contents if it already
    exists!
    
    Optionally output something special in  <tagline>
    """
    # check if it exists
    print 'Creating a new directory ' + os.path.relpath( dir_name ) + tagline
    if os.path.isdir( dir_name ):
        print 'a directory named ' + os.path.relpath( dir_name ) + ' already exists, deleting it now...'
        shutil.rmtree( dir_name )
    os.mkdir( dir_name )

# copy helper
def copy_file( filename , destination , display = False ):
    """
    Copy  <filename>  to/into  <destination>
    
    just a cp wrapper...what?
    """
    if display:    # optional
        if os.path.isdir( destination ):
            print 'placing a copy of ' + os.path.relpath( filename ) + ' into the ' + os.path.relpath( destination ) + ' directory'
        elif os.path.isfile( destination ):
            print 'copying ' + os.path.relpath( filename ) + ' to ' + os.path.relpath( destination )
    shutil.copy( filename , destination )

################################################################################
# SEQUENCE HANDLING HELPERS

# basic converters...its done a lot

# loading wrapper...basically cause "from Bio import SeqIO" is too long
def load_sequence( filename , ignore_empty = True , seqformat_map = SEQFORMAT_MAP ):
    """
    Returns the list of sequences in  <filename>  as Biopython SeqRecord
    objects
    automatically handles different file format as specified by  <seqformat_map>
    
    Optionally  <ignore_empty>  sequences (SeqID in file but no sequence)
    
    To get string, use get_sequence
    """
    # determine the file format
    seq_format = get_file_extension( filename )

    # load ALL the sequences!
    sequences = [i for i in SeqIO.parse( filename , seqformat_map[seq_format] )]
    if ignore_empty:
        sequences = [i for i in sequences if str( i.seq )]

    # or just one...
    if len( sequences ) == 1:
        sequences = sequences[0]
    
    return sequences

# general converter!
def get_sequence( sequence , seq_format = SEQFORMAT , uppercase = True , ignore_empty = True , get_ids = False ):
    """
    Returns a string or list of string depending on the input  <sequence>
    can accept:
        a filename for a  <seq_format> file
        a Biopython Seq object
        a Biopython SeqRecord object
        a string
        a list of any of the above (can be heterogenous)
        
    Optionally change the sequence to  <uppercase>  (ambiguous positions are
    sometimes lowercase)
    Optionally <ignore_empty> sequences (SeqID in file but no sequence)
    Optionally  <get_ids>  , returning a parallel list of SeqIDs and descriptions
    """
    # sort the input data type
    # for common Biopython objects
    if type( sequence ) == Seq:
        sequence = str( sequence )
    elif type( sequence ) == SeqRecord:
        seq_ids = str( sequence.id )
        seq_des = str( sequence.description )
        sequence = str( sequence.seq )

    # input file
    elif '.' in sequence:    # should never occur!
        # its a filename (?) so try to load it, it will error properly
        sequence = load_sequence( sequence , ignore_empty )
        # sort by number
        if type( sequence ) == list:    # in accordance with the above
            # optionally get the ids
            if get_ids:
                seq_ids = [str( i.id ) for i in sequence]
                seq_des = [str( i.description )*( not i.description == i.id ) for i in sequence]
            sequence = [str( i.seq ) for i in sequence]
        else:
            if get_ids:
                seq_ids = str( sequence.id )
                seq_des = str( sequence.description )*( not sequence.description == sequence.id )
            sequence = str( sequence.seq )

    # list of any of the above
    elif type( sequence ) == list:
        # then sort based on individual types...
        sequence = [get_sequence( i , seq_format , uppercase , ignore_empty , get_ids ) for i in sequence]
        if get_ids:
            seq_ids = [i[1] for i in sequence]
            seq_des = [i[2] for i in sequence]
            sequence = [i[0] for i in sequence]
    
    # should be an input single string
    else:
        seq_ids = ''
        seq_des = ''
    
    # optionally force UPPER case
    if uppercase:
        if type( sequence ) == str:
            # single sequence
            sequence = sequence.upper()
        else:
            # multiple
            sequence = [i.upper() for i in sequence]
    
    # optionally return the id and descriptions too
    if get_ids:
        return sequence , seq_ids , seq_des
    return sequence

# general writer
# return the filename
def write_sequence( sequence , out_filename = '' , seq_format = SEQFORMAT , seqid = 'unknown' , description = '' , alphabet = DNAAlphabet , seq_format_map = SEQFORMAT_EXTENSION_MAP ):
    """
    Write  <sequence>  to  <out_filename>  as  <seq_format>  using  <alphabet>
    
    Robust to sequence inputs that are:
        str (filename or sequence)
        Seq
        SeqRecord
    """
    # sort the input data type
    unknown = 1
    # for common Biopython objects
    if isinstance( sequence , str ):
        if '.' in sequence:    # should never occur, okay, I made it occur
            print 'it appears you input a path or filename...so its already a file!'
            return sequence
        sequence = SeqRecord( Seq( sequence , alphabet ) )    # already default ID of unknown
        sequence.id = seqid
        sequence.description = description
    elif isinstance( sequence , unicode ):    # hacky, unicode vs str
        sequence = str( sequence )
        if '.' in sequence:    # should never occur
            print 'it appears you input a path or filename...so its already a file!'
            return sequence
        sequence = SeqRecord( Seq( sequence , alphabet ) )    # already default ID of unknown
        sequence.id = seqid    
        sequence.description = description
    elif isinstance( sequence , Seq ):
        sequence = SeqRecord( sequence )
        sequence.id = seqid
        sequence.description = description
    elif isinstance( sequence , list ):
        # yay, do it all over again :(
        # make recursive

        # assume all members are the same type...else its an error anyway
        if isinstance( sequence[0] , str ):
            for i in xrange( len( sequence ) ):
                sequence[i] = SeqRecord( Seq( sequence[i] , alphabet ) )
                sequence[i].id = seqid + '_' + str( unknown )
                sequence[i].description = description
                unknown += 1
        elif isinstance( sequence[0] , Seq ):
            for i in xrange( len( sequence ) ):
                sequence[i] = SeqRecord( i )
                sequence[i].id = seqid + '_' + str( unknown )
                sequence[i].description = description
                unknown += 1

    # now that all are Biopython SeqRecords, write to file!
    if not out_filename:
        if type( sequence ) == list:
            out_filename = sequence[0].id + '.' + seq_format_map[seq_format]
        else:
            out_filename = sequence.id + '.' + seq_format_map[seq_format]
    
    SeqIO.write( sequence , out_filename , seq_format )
    print 'Successfully wrote the sequence(s) to ' + os.path.relpath( out_filename )
    
    return out_filename

################################################################################
# FULL RAW PROCESSING

# 1R69 - single model, single chain
# 1A17 - another random choice for testing
# 1BUW
# 1C17
# 1JYX
# 1M2V
# 1TF6
# 2C35
# 3G3O
# 1YY8 - AB and CD, single model
# 1NMR - multiple models
# 1LR1 - multiple models AND chains
# 1VTL - protein and DNA, single model
# 1UN6 - protein and RNA, single model

# the big boy...
def process_pdb( pdb_filename , seqformat = SEQFORMAT , seqformat_extension_map = SEQFORMAT_EXTENSION_MAP , conversion = three2three , na_conversion = NA_CONVERSIONS_ROSETTA , na_alphabet = DNAAlphabet , protein_alphabet = ProteinAlphabet ):
    """
    Create a directory from  <pdb_filename>  containing relevant information
    stored in the PDB file
    
    This method behaves slightly differently for PDB files with multiple models,
    nucleic acids, duplicate complexes, etc.
    so if you are interested in the specifics, please read the source code
    
    In short, it tries to write:
        header.txt          a text file of the header lines
        numbering_map.txt   a text file showing 1-indexed PDB numbering
        clean.pdb           only ATOM lines
        hetatm.pdb          only HETATM lines, may be split by resName
        .fa                 sequences of all peptides and nucleic acids
        subdirectories      for each protein model/subunit (similar info)
    
    does not write a text file for the "trailer" (lines after the coordinates)
    
    converts lines (ATOM or HETATM) that can be converted based on  <conversion>
    (generally) and  <na_conversion>  (specific for nucleic acids, relevant
    because RNA and DNA may require different treatment...)
    !!!WARNING!!! defaults:
        CSE     CYS     converts SelenoCysteinE to Cysteine
        HYP     PRO     converts HYdroxylProline to Proline
        CYD     CYS     does NOT convert "CYsteine Disulfides to Cysteine"
        HIP     HIS     converts "HIP" to Histidine (~double protonation)
        HID     HIS     converts "HID" to Histidine (~single delta N proton)
        HIE     HIS     converts "HIE" to Histidine (~single epsilon N proton)

    todo:
    ensure hetatm conversions step illegal atoms!!!!
    alternate conformations
    convert DNA to Rosetta DNA
    convert ligands to params
    convert water to TP3 (or TP5)
    """
    # process input, optionally a list
    if isinstance( pdb_filename , list ):
        print 'Multiple PDB codes detected, processing them individually...'
        # use this list comprehension, get them all!
        filenames = [process_pdb( i , seqformat , seqformat_extension_map , conversion , na_conversion , na_alphabet , protein_alphabet ) for i in pdb_filename]
        print 'Finished the whole list, enjoy!'
        return filenames
    
    ####################
    # NEW DIRECTORY ETC.
    # get root name
    pdb_filename = os.path.abspath( pdb_filename )
    root_name = get_root_filename( pdb_filename )
    best_guess = pdb_filename
            
    # make a new directory, a whole lot is gonna go here...
    create_directory( root_name , ' to sort the data' )
    # move the pdb here
    copy_file( pdb_filename , root_name )

    # oh, and go there too
    original_dir = os.getcwd()
    os.chdir( root_name )

    # "update" the target
    pdb_filename = root_name + '/' + os.path.split( pdb_filename )[-1]
    root_name = get_root_filename( pdb_filename )
    
    ##############
    # PRE CLEANING
    # does not need to know if nucleics or not
    
    # convertions!
    # ...bad...overwrite the file!...but no filename management
    convert_pdb_resnames_to_ATOM_lines( pdb_filename , pdb_filename , root_name +'_conversion_report.txt' , conversion )
    
    # produce a PDB with just the protein lines
    best_guess = clean_ATOM_lines_from_pdb( pdb_filename )
    
    # extract numbering
    # don't bother storing the map
    extract_numbering_map_from_pdb( pdb_filename , 'numbering_map.txt' )
    
    # extract HETATM lines
    clean_HETATM_lines_from_pdb( pdb_filename )
    
    # write out alternate conformations for the cleaned file
    alternate_conformations = clean_alternate_conformations_from_pdb( best_guess )
    
    ##########################
    # HEADER PARSING
    # extract info from header

    # this information is accessible from the PDBParser header...sorta...    
    # get the number of models
    models = extract_number_of_models_from_pdb_header( pdb_filename )
    
    # get the subunit complexes
    complexes = extract_duplicate_chains_from_pdb_header( pdb_filename )
        
    # write the header (?)
    # get the header
    header = extract_header_from_pdb( pdb_filename )

    ###################
    # HUNT DOWN HETATMS
    
    # use the map in the header and extracted chemical formulas to search pubchem
    
    # get map
    
    # per hetatm type
    # get formula
    # get number of residues -> needed to interpret formula...
    # search pubchem, download best sdf if exact match and at least < atoms
    # create directory for these params etc.

    ##########################
    # ASSESS NUCLEIC SITUATION

    # HERE!
    
    # choose your fate!, removes nucleic lines
    has_nucleic = clean_nucleic_acid_lines_from_pdb( pdb_filename )
    
    # get proteins if nucleics
    if has_nucleic:
        # get a PDB of protein only, use this from now on
        print 'Scanners indicate there are nucleic acid lines in ' + os.path.relpath( pdb_filename ) + '\nSadly, a lot of toys do not play well with these so a few extra steps are required...'

        # write nucleic sequences
        temp , nucleic_types = extract_nucleic_acid_sequences_from_pdb( root_name + '.nucleic.pdb' , seqformat = seqformat , alphabet = na_alphabet , seqformat_extension_map = seqformat_extension_map )
        # care not for the sequences
        
        # make a Rosetta ready nucleic PDB!!!
        # SO BAD! overwrite!
        # BAH!!!
        na_chains = split_pdb_into_chains( root_name + '.nucleic.pdb' , 0 , True )    # just 0 model...
        for i in na_chains.keys():
            # BETTER BE IN BOTH!!!
            convert_pdb_resnames_to_ATOM_lines( na_chains[i] , na_chains[i] , 'nucleic_chain_'+ i +'_conversion_report.txt' , na_conversion[nucleic_types[i]] )
        
        # check for protein :)
        has_protein = clean_protein_lines_from_pdb( pdb_filename )
        if not has_protein:
            print 'The additional features are only available for proteins\nScanner indicate that this PDB has ONLY nucleic acids (no proteins) :(\nthe remaining methods rely on the Biopython PDBParser...and things get messy with nucleic acids\nEven so, the only feature you\' missing out on is splitting into subdirectories for each chain, and since the PDB is just nucleic acid, that isn\'t as helpful'

            # premature exit
            os.chdir( original_dir )
            return best_guess

        # change the name of the best guess to .protein.pdb
        best_guess = root_name + '.protein.pdb'
        pdb_filename = root_name + '.protein.pdb'
        
        # get the nucleic chains
        nucleic_chains = extract_chains_from_pdb( root_name + '.nucleic.pdb' )
    
    ############
    # PDB PARSER
    
    # does NOT loop over ANY nucleic acid chains!
    
    # prepare to load...
    parser = PDBParser( PERMISSIVE = 1 )
    writer = PDBIO()
    
    struct = parser.get_structure( root_name , pdb_filename )

    # verify models and chains
    temp = len( struct.child_list )    # number of models
    if not temp == models:
        print 'Huh? the PDB file header claims there are ' + str( models ) + ' models but the PDB file has ' + str( temp ) + ' models...\nUsing the ACTUAL number of models (' + str( temp ) + ')'
        models = temp

    # check from reading the CHAIN
    if not complexes:
        print 'No chain/subunit information found in the header (or no header),\nassuming all individual sequences are unique i.e. if AB and copy CD, will make A, B, C, and D instead of AB and CD'
#        complexes = temp    # unecessary, automatically happens below...

    # add all new ids
    temp = struct[0].child_dict.keys()    # it better have at least 1 model...

    
    # for the nucleic case...
    if has_nucleic:
        # HERE!
        # remove nucleic lines...
        for i in xrange( len( complexes ) ):
            for j in nucleic_chains:
                if j in complexes[i]:
                    complexes[i] = complexes[i].replace( j ,'' )
        
        # sanity check...
        complexes = [i for i in complexes if i]
        
        # assume all models contain all chains...idk how this would ever NOT occur...
        # this also produces a directory for EACH chain as the default behavior!!!
        complexes += [i for i in temp if i and not i in complexes and not i in nucleic_chains]

    else:
        # normal protein stuff
        complexes += [i for i in temp if i and not i in complexes]
    
    
    # okay...this should be figured out...but isn't that big of a deal
    # found with 1JGO
#    print complexes
#    print complexes
#    complexes = [i for i in complexes if i]
    
#    input('dd')
    
    ################################
    # CREATE AND FILL SUBDIRECTORIES
    
    # again, this step is skipped for pure nucleic acid...
                
    # exit condition, only 1 model and 1 chain
    if models > 1 or len( complexes ) > 1:
        # over the models
        for model in struct.child_dict.keys():
            # over the chains
            for complx in complexes:
#                print '='*60 + complx
                # remove nucleic subunits
                # HERE!
                if has_nucleic:
                    for chain in nucleic_chains:
                        complx = complx.replace( chain , '' )    # delete the chain from the complex
            
                # check that all members are present
                chains = struct[model].child_dict.keys()
                missing = [l for l in complx if not l in chains]
                # report this!
                if missing:
                    # add models bool for str here?
                    print 'Expected model ' + str( model + 1 ) + ' to have chains ' + complx + ' but the its missing chains ' + ', '.join( missing ) + '!'
                
                # create the new directory
                # only number if more than 1 model
                dir_name = complx + str( model + 1 )*bool( models - 1 )
                new_dir = os.path.split( root_name )[0] + '/' + dir_name
                print 'Creating the subdirectory ' + os.path.relpath( new_dir )
                os.mkdir( new_dir )

                # create a copy of the complex, only the chains of interest
                # make an empty structure
                temp = Structure( 'temp' )
                temp_model = Model( model )    # and an empty model
                temp.add( temp_model )
                
                # add the complex
                for chain in complx:
                    temp[model].add( struct[model][chain] )
                    
                    # get the chain sequence
                    seqid = dir_name + ('_model_' + str( model + 1 ))*bool( models - 1 ) + '_chain_' + chain
                    seq_filename = new_dir + '/' + os.path.split( root_name )[-1] + ('_model_' + str( model + 1 ))*bool( models - 1 ) + '_chain_' + chain + '.' + seqformat_extension_map[seqformat]
                    description = '(from model ' + str( model + 1 ) + ')'
                    temp_seq = extract_protein_sequence_from_pdb( temp , True ,    # MUST insert disorder...
                        seq_filename , seqid , description , model , chain ,
                        True , seqformat , protein_alphabet , seqformat_extension_map )

                    # also, make sure at least one copy (from the first model) is in the main dir
                    seq_filename = root_name + '_chain_' + chain + '.' + seqformat_extension_map[seqformat]
                    if not os.path.exists( seq_filename ):
                        print 'Putting a copy of the sequence in the new directory'
                        # assumes all the models have the same sequence
                        write_sequence( temp_seq , seq_filename , seqformat ,
                            os.path.split( root_name )[-1] + ' chain ' + chain ,
                             description , protein_alphabet , seqformat_extension_map )

                # write out the model+chain
                writer.set_structure( temp )
                print 'Writing a copy of model ' + str( model + 1 ) + ' chain(s) ' + complx + ' to ' + new_dir + '.pdb'
                writer.save( new_dir + '/' + dir_name + '.pdb' )#, selection )
                
                # also write a cleaned PDB file, onlt ATOM lines
                clean_ATOM_lines_from_pdb( new_dir + '/' + dir_name + '.pdb' )
                
                # also write any alternate conformations
                clean_alternate_conformations_from_pdb( new_dir + '/' + dir_name + '.pdb' )
            
                # also get specific HETATMs...this is getting bulky...
                clean_HETATM_lines_from_pdb( new_dir + '/' + dir_name + '.pdb' )
                
                # no need to clean DNA
    else:
        # only 1 model AND only 1 chain
        # still write it please :)
        model = 0
        chain = complexes[0]
        
        # may seem silly, but this edge case will prevent needless re-parsing

        # get the chain sequence
        seqid = os.path.split( root_name )[-1] + '_chain_' + complexes[0]

        extract_protein_sequence_from_pdb( struct , True ,
                seqid + '.' + seqformat_extension_map[seqformat] , seqid , '' ,
                model , chain , True ,
                seqformat = seqformat , alphabet = protein_alphabet , seqformat_extension_map = seqformat_extension_map )
    
    # debug summary...
    temp = os.listdir( os.getcwd() )
    temp.sort()
    print 'New Files in the ' + root_name + ' directory :\n' + '\n'.join( ['\t'+ i for i in temp] )
    
    # return back one directoy
    os.chdir( original_dir )    # yeah...its hacky
    
    return best_guess
    
################################################################################
# HEADER STUFF

# extract header text
def extract_header_from_pdb( pdb_filename , header_filename = 'header.txt' ):
    # write the header (?)
    # get the header
    f = open( pdb_filename , 'r' )
    header = ''
    while True:    # should error from f.next() if improper input...
        # next line
        line = f.next()
        
        # exit condition
        if 'ATOM' == line[:4] or 'MODEL' == line[:5] or 'HETATM' == line[:6]:
            break
        
        header += line
    f.close()
    
    # write the header
    if header_filename:
        print 'Writing a copy of the header lines to the file ' + header_filename
        f = open( header_filename , 'w' )
        f.write( header )
        f.close()

    return header

# return any predicted shain pairs
def extract_duplicate_chains_from_pdb_header( pdb_filename ):
    # load the raw data
    f = open( pdb_filename , 'r' )
    complexes = []
    keep_going = True
    while keep_going:
        # next line
        line = f.next()

        # ...think about this...
        # check if chain info, extract the matching subunits
        if line[:6] == 'COMPND' and 'CHAIN:' in line:
            duplicate = line.split( 'CHAIN: ' )[-1].replace( ';' , '' ).strip().split( ', ' )    # ignore ";\n"
            if len( duplicate ) > 1:
                complexes.append( duplicate )
        # stop condition
        elif not ('HEADER' in line or 'TITLE' in line or 'COMPND' in line or 'CAVEAT' in line):
            keep_going = False
    f.close()

    # convert complexes
    if complexes:
        if not sum( [len( c ) - len( complexes[0] ) for c in complexes] ):
            # all are the same length
            complexes = [''.join( [c[i] for c in complexes] ) for i in xrange( len( complexes[0] ) )]
        else:
            # uh oh...
            # could be all should be unique...which puts us in exception land anyway
            # assume that last listed are aberrantly unpaired
            lowest = min( [len( c ) for c in complexes] )
            temp = [''.join( [c[i] for c in complexes] ) for i in xrange( lowest )]
            for c in complexes:
                temp += c[lowest:]
            complexes = temp
    
    return complexes

# return number of models, scanned from header
def extract_number_of_models_from_pdb_header( pdb_filename ):
    # get the number of models
    f = open( pdb_filename , 'r' )
    models = 1
    keep_going = True
    while keep_going:
        # next line
        line = f.next()
        
        # check for models
        if line[:6] == 'NUMMDL':
            models = int( line.replace( 'NUMMDL' , '' ).strip() )
            keep_going = False
        elif line[:4] == 'ATOM':
            keep_going = False
    f.close()
    
    return models

# return resolution, scanned from header
# other information? R-value? R-free?
# other places to extract the quality...?
def extract_resolution_information_from_pdb_header( pdb_filename ):
    # load it
    f = open( pdb_filename , 'r' )
    
    # ewww....should be a "for" loop that breaks...
    keep_going = True
    experimental_data = 'X-RAY DIFFRACTION'
    resolution = None
    while keep_going:
        # next line
        line = f.next()
        
        # check for models
        if line[:6] == 'EXPDTA':
#            print 'found exp data'
            experimental_data = line[6:].strip()
        elif line[:10] == 'REMARK   2':
            # check for NMR
#            print 'found remark'
#            print line
            if 'ANGSTROMS' in line:
#                print 'found resolution'
                resolution = float( line[23:].strip().split( 'ANGSTROMS' )[0].strip() )
                keep_going = False
        elif line[:4] == 'ATOM':
            keep_going = False
    f.close()
    
    return resolution , experimental_data

# return number of models, scanned from header
def extract_HETNAM_from_pdb_header( pdb_filename ):
    # get the number of models
    f = open( pdb_filename , 'r' )
    hetname_map = {}
    keep_going = True
    while keep_going:
        # next line
        line = f.next()
        
        # check for models
        if line[:6] == 'HETNAM':
            hetname = line[6:].strip().split( ' ' )
            hetkey = hetname[0]
            hetname = ''.join( [i + ' ' for i in hetname[1:]] )[:-1]
            
            hetname_map[hetkey] = hetname
        elif line[:4] == 'ATOM':
            keep_going = False
    f.close()
    
    return hetname_map

################################################################################
# DIVIDE AND JOIN

# split or join PDB files

# simple wrapper
def morph_atomName2element( atomName ):
    """
    Returns the element in  <atomName>
    
    raw PDB atomNames are supposed to have the element as the first character
    """
    element = atomName[:2].strip()
    
    # remove number characters
    for i in '0123456789':
        element = element.replace( i , '' )
    
    return element

# make sure a filename, Structure, or Model returns the Model of interest
# not tested recently...
def load_pdb( pdb , model = 0 ):
    """
    Returns the  <model>  of  <pdb>  if its a Structure object (or a filename)
    """
    # sort the input
    if isinstance( pdb , str ):
        # filename
        print 'Input filename ' + pdb + ', loading the structure now'
        parser = PDBParser( PERMISSIVE = 1 )
        pdb = parser.get_structure( 'temp' , pdb )

        # default to first one if empty...
        if not model:
            model = pdb.child_dict.keys()[0]

        print 'extracting the first model (' + str( model ) + ')'
        pdb = pdb[model]    # get the first model

    # tried doing this a prettier way...
    # check for specific methods and data types for clues...
    elif isinstance( pdb.child_dict.keys()[0] , int ):
        # its a Biopython structure
        # default to first one if empty...
        if not model:
            model = pdb.child_dict.keys()[0]

        print 'Input Biopython Structure, extracting the first model (' + str( model ) + ')'
        pdb = pdb[model]    # get the first model

    elif 'child_dict' in dir( pdb ):
        # ...could be any number of things...including what we want!
        # hooray! everything is okay
        None
    else:
        # not supported!
        raise IOError( 'That data structure is not currently supported...' )
    
    return pdb

# check the PDB for models and split into separate PDBs
def split_pdb_into_models( pdb_filename ):
    """
    Writes a single PDB file for every model in  <pdb_filename>
    uses the Biopython PDBParser and PDBIO
    """
    # make tools
    parser = PDBParser( PERMISSIVE = 1 )
    writer = PDBIO()

    pdb_filename = os.path.abspath( pdb_filename )
    root_name = get_root_filename( pdb_filename )
    struct = parser.get_structure( root_name , pdb_filename )

    # over the models
    for i in struct.child_dict.keys():
        # get just the model
        temp = Structure( 'temp' )
        temp.add( struct[i] )

        # write it
        writer.set_structure( temp )
        out_filename = root_name + '_model_' + str( i + 1 ) + '.pdb'
        print 'Model ' + str( i + 1 ) + ' written to ' + out_filename
        writer.save( out_filename )

# check the PDB for chains and split into separate PDBs
def split_pdb_into_chains( pdb_filename , model = 0 , export = False ):
    """
    Writes a single PDB file for every chain in  <pdb_filename>
    uses the Biopython PDBParser and PDBIO
    """
    # make tools
    parser = PDBParser( PERMISSIVE = 1 )
    writer = PDBIO()

    pdb_filename = os.path.abspath( pdb_filename )
    root_name = get_root_filename( pdb_filename )
    struct = parser.get_structure( root_name , pdb_filename )
    
    # assume there is only 1 model
    # over the chains
    chains = {}
    for i in struct[model].child_dict.keys():
        # get just the model
        temp = Structure( 'temp' )
        temp_mod = Model( 0 )
        temp_mod.add( struct[0][i] )
        temp.add( temp_mod )

        # write it
        writer.set_structure( temp )
        out_filename = root_name + '_chain_' + i + '.pdb'
#        chains.append( 'Chain ' + i + ' written to ' + out_filename )
        chains[i] = out_filename
        writer.save( out_filename )

    # debug output
    for i in chains.keys():
        print 'Chain ' + i + ' written to ' + chains[i]

    # optionally export
    if export:
        return chains

# add all files together in the provided order
# not tested recently...
def join_pdb_files( files , out_filename = '' ):
    """
    Combines the contents of all  <files>  and writes it out to  <out_filename>
    
    a very simple method
    """
    # default filename
    out_filename_provided = True
    if not out_filename:
        out_filename_provided = False
    
    text = ''
    for i in files:
        # open it
        f = open( i , 'r' )
        
        # add the text
        text += f.read()
        f.close()
    
        # check if the name should be added
        if not out_filename_provided:
            if '.' in i:
                out_filename += i[:i.find( '.' )]
            else:
                out_filename += i

    # write the bastard love child
    f = open( out_filename , 'w' )
    f.write( text )
    f.close()

# extract the chains from the PDB
# only considers ATOM lines, mainly for use with clean_nucleic_acid_lines_from_pdb
def extract_chains_from_pdb( pdb_filename , only = ['ATOM'] ):
    """
    Returns the chains found in  <pdb_filename>
    
    Only consider lines starting with  <only>
    """
    pdb_filename = os.path.abspath( pdb_filename )
    if os.path.exists( pdb_filename ):
        # load the data
        f = open( pdb_filename , 'r' )
        data = [i for i in f.xreadlines() if i[:6].strip() in only]
        f.close()
    
        # find unique chains
        chains = []
        for i in data:
            if not i[21] in chains:
                chains.append( i[21] )

        return chains
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False
    
# extract the name mapping
def extract_numbering_map_from_pdb( pdb_filename , out_filename = '' , only = ['ATOM'] ):
    """
    Returns a map (dict) from residues in  <pdb_filename>  that are 1-indexed
    and a reverse map (dict)
    
    Only consider lines starting with  <only>
    
    Optionally write the results to  <out_filename>
    """
    pdb_filename = os.path.abspath( pdb_filename )
    if os.path.exists( pdb_filename ):
        # load the raw data
        f = open( pdb_filename , 'r' )
        d = [i for i in f.xreadlines() if i[:6].strip() in only]
        f.close()
    
        # extract dict of pairs
        pdb_map = {}
        reverse_map = {}
        count = 0
        text = ''
        for i in d:
            # basic info
            chain = i[21]
            resseq = i[22:26].strip()
            icode = i[26]    # the icode
            key = chain + resseq + icode
        
            if not key in pdb_map.keys():
                count += 1
                pdb_map[key] = count
                reverse_map[count] = key
            
                text += key + '\t' + str( count ) + '\n'

        # optionally write to file
        # no defaulting!
        if out_filename:
            # default filename
    #        f = open( get_root_filename( pdb_filename ) + '_PDB_numbering.txt' , 'w' )
    #        f.write( ''.join( [i +'\t'+ str( pdb_map[i] ) +'\n' for i in pdb_map.keys()] )  )
            print 'Writing the PDB numbering of ' + pdb_filename + ' to ' + out_filename
            f = open( out_filename , 'w' )
            f.write( text )
            f.close()
    
        return pdb_map , reverse_map
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

# extract a protein sequence from a PDB
# make this better? specify the chain?
# for now, only works if single chain...
def extract_protein_sequence_from_pdb( pdb , include_breaks = True ,
        out_filename = '' , seqid = '' , description = '' ,
        model = 0 , chain = 'A' , export = True ,
        seqformat = SEQFORMAT , alphabet = ProteinAlphabet ,
        seqformat_extension_map = SEQFORMAT_EXTENSION_MAP ):
    """
    Returns the protein sequences found in  <pdb>  in  <model>
    
    Optionally  <export>  the sequence
    Optionally write to  <out_filename>  with  <seqid>
    
    note: does NOT scan for protein chains, it only dumps out the full
        protein sequence in the PDB file
        individual chains can be extracted using process_pdb
    """
    # ensure pdb is proper, must be a model
    pdb = load_pdb( pdb , model )    # necessary model?

    # format the chain input
    if not isinstance( chain , list ):
        chain = [chain]
   
    # over the desired chains
    # ugh...this should all be rewritten...
    sequences = []
    for c in chain:
        # get it
        if include_breaks:
            # extract the sequence as a Biopython Seq object
            # convert the model into a Structure, for getting the sequence
            for_seq = Structure( 'temp' )
            # ...oh yeah...must be model 0 and up >:[
            temp_model = Model( 0 )    # hardcoded...
            for_seq.add( temp_model )
#            for ch in pdb.child_dict.keys():
            # copy it all over directly
            for_seq[0].add( pdb[c] )

            # gap regions makred as "|"
            seq_builder = PPBuilder()
            pp = seq_builder.build_peptides( for_seq )
            seq = Seq( '|'.join( [str( frag.get_sequence() ) for frag in pp] ) , alphabet )
#            for frag in pp:
#                seq += frag.get_sequence() + '|'    # already a Biopython Seq
            seqr = SeqRecord( seq )
            seqr.description = description + ' missing residues (gap regions as \"|\")'*( '|' in seq )
        else:
            # just iterate and extract!
            seq = Seq( ''.join( [three2one[i.resname] for i in pdb.get_residues() if i.resname in three2one.keys() and i.get_parent().id == c] ) , alphabet )
            seqr = SeqRecord( seq )
            seqr.description = description
    
        # prepare to write
    #    seq = seq[:-1]
    #    seqr.description = 'missing residues (gap regions as \"|\")'*( '|' in seq )    # no need if no gaps
        seqr.id = seqid
        sequences.append( seqr )

    # optionally write the sequence
    if out_filename:
        write_sequence( sequences , out_filename , seqformat , alphabet , seqformat_extension_map )
    
    # optionally export the sequence
    if export:
        return get_sequence( sequences )
#        return str( seq )

# extract and write a file from the PDB
def extract_nucleic_acid_sequences_from_pdb( pdb_filename , out_filename = '' , NA = NUCLEIC_SEQUENCE_LETTERS_MAP , DNA = NA_CODES['DNA'] , seqformat = SEQFORMAT , alphabet = DNAAlphabet , seqformat_extension_map = SEQFORMAT_EXTENSION_MAP ):
    """
    Returns the protein sequences found in  <pdb_filename>
    
    Only consider resNames in  <NA>
    
    Optionally write to  <out_filename>
    """
    pdb_filename = os.path.abspath( pdb_filename )
    if os.path.exists( pdb_filename ):
        # load the data
        f = open( pdb_filename , 'r' )
        d = f.readlines()
        f.close()
    
        # print about fails/assumptions
        print 'Extracting nucleic sequences from ' + os.path.relpath( pdb_filename ) + '\nFor visibility, this method assumes A LOT!\n1. nucleotides are identified by a unique resSeq codes (with a proper resName)\n2. sequences are identified by unique chain IDs\n3. RNA is the default state\n4. DNA is identified by \"DG\" (etc.) OR \"T\" resi codes\n4. All sequences are continuous\n6. All sequences are recorded 5\' -> 3\' (and written to file in this order)'
    
        # check for nucleic lines - No, do while parsing
        # extract sequence
        NA_keys = NA.keys()
        DNA_keys = DNA.keys()
#        molecule = 'RNA'
        molecule_types = {}
        sequences = {}
        last = None
        for line in d:
            # must have C1 and a nucleic resi code to be considered a nucleotide
            resname = line[17:20].strip()
            resseq = line[22:27].strip()    # resseq
            if (line[:5] == 'ATOM ' or line[:4] == 'TER ') and resname in NA_keys:# and line[13:16].strip() == 'C1\'':
                # only novel lines
                if resseq == last:
                    continue
                last = resseq    # if the remainder will execute...
            
                # check for DNA
                chain = line[21]
                if [True for i in DNA_keys if i in resname]:
                    # its DNA
                    molecule_types[chain] = 'DNA'
                    # consider the whole chain DNA if ANY of the exclusive codes are present
                    # sometimes DNA is abbreviated without the "d" to designate "deoxy"
            
                # remember the letter
                if chain in sequences.keys():
                    # add the letter
                    sequences[chain] += NA[resname]    # map the code
                else:
                    # create it as well
                    sequences[chain] = NA[resname]
                    molecule_types[chain] = 'RNA'    # default

        # default out name
        root_filename = get_root_filename( pdb_filename )
        if not out_filename:
            out_filename = root_filename

        # write the sequences
        for chain in sequences.keys():
            # verify its not just a nucleotide
            seq = sequences[chain]
            if len( seq ) > 1:
                # determine the molecule type
                # record a proprt id
                seqr = SeqRecord( Seq( seq , alphabet ) )    # even if RNA (?)
                seqr.id = os.path.split( root_filename )[-1] + '_chain_' + chain
                seqr.description = molecule_types[chain]
            
                # oh yeah, write it, prints out by itself
                out_filename = seqr.id + '.' + seqformat_extension_map[seqformat]
                write_sequence( seqr , out_filename , seqformat , alphabet , seqformat_extension_map )

        return sequences , molecule_types    # empty dict will evaluate as false
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

################################################################################
# CLEANING METHODS

# HERE !!!

# a dirty input produces a cleaned output file :)
# default behavior is to produce output

# removes non ATOM lines from  <pdb_file>  and writes to  <out_file>
def clean_ATOM_lines_from_pdb( pdb_filename , out_filename = '' , HETATM_include = [] , excluded_atoms = ['CN'] , accepted_fields = ['ATOM ' , 'TER '] ):
    """
    Writes all lines in the PDB file  <pdb_filename>  beginning with "ATOM" or
    "TER" into  <out_filename>  (defaults to  <pdb_file>.clean.pdb)

    Optionally include HETATM lines with resNames in  <HETATM_include>

    Returns True if successful
    
    ...pretty much the same as:
    
        grep "ATOM" pdb_filename > out_filename

    example:
        clean_non_ATOM('1YY9.pdb')
    See also:
        Pose
        Pose.dump_pdb
        pose_from_pdb
        pose_from_rcsb
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # an optional argument for PDB files not ending in .pdb
#    if not edit:
#        edit = 255

    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()
        good = []
        for i in data:
            if [True for j in accepted_fields if i[:len( j )] == j]:
#            if i[:5] == 'ATOM ' or i[:4] == 'TER ':
                # add your preference rules for ligands, DNA, water, etc.
                # check for excluded atoms
                if i[12:16].strip() in excluded_atoms:
                    # skip it, do not add to the list
                    continue
                
                good.append( i )
            elif i[:6] == 'HETATM' and i[17:20] in HETATM_include:
                # save for later, more processsing
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No ATOM or HETATM lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.clean.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()

        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, non-ATOM lines removed\nclean data written to ' + os.path.relpath( out_filename )
        return out_filename

    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

# if you would prefer a simpler call using grep, it looks something like this
#    os.system("grep \"ATOM\" %s.pdb > %s.clean.pdb"%(pdb_file[:edit],pdb_file[:edit]))

# split the ATOM lines, only look for DNA lines
def clean_nucleic_acid_lines_from_pdb( pdb_filename , out_filename = '' , NA = NUCLEIC_SEQUENCE_LETTERS_MAP.keys() ):
    """
    Scan  <pdb_filename>  for any nucleic acid lines and writes these to
    <out_filename>
    
    defines nucleic acid resNames (three letter codes) as those with
    stripped codes in  <NA>
    
    default definition of nucleic acid resNames can be adjusted in settings.py
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()

        good = []
        for i in data:
            if (i[:5] == 'ATOM ' or i[:4] == 'TER ') and i[17:20].strip() in NA:
                # add your preference rules for ligands, DNA, water, etc.
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No nucleic acid lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.nucleic.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()
        
        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, DNA/RNA lines extracted\nclean data written to ' + os.path.relpath( out_filename )
        return out_filename

    else:
        print 'No such file or directory named '+ os.path.relpath( pdb_filename )
        return False

# split the ATOM lines, only look for not RNA/DNA lines
def clean_protein_lines_from_pdb( pdb_filename , out_filename = '' , NA = NUCLEIC_SEQUENCE_LETTERS_MAP.keys() ):
    """
    Scan  <pdb_filename>  for any nucleic acid lines and writes all "ATOM" lines
    that are NOt nucleic acids to  <out_filename>
    
    defines nucleic acid resNames (three letter codes) as those with
    stripped codes in  <NA>
    
    default definition of nucleic acid resNames can be adjusted in settings.py
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()

        good = []
        for i in data:
            if (i[:5] == 'ATOM ' or i[:4] == 'TER ') and not i[17:20].strip() in NA:
                # add your preference rules for ligands, DNA, water, etc.
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No protein lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.protein.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()
        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, protein lines extracted\nclean data written to ' + os.path.relpath( out_filename )
        return True

    else:
        print 'No such file or directory named '+ os.path.relpath( pdb_filename )
        return False

# scan for HETATMs, rewrite without all these lines, record specific ones
def clean_HETATM_lines_from_pdb( pdb_filename , out_filename = '' , only = '' , write_unique = True ):
    """
    Writes all lines in the PDB file  <pdb_filename>  beginning with "HETATM"
    into  <out_filename>  (defaults to  <pdb_filename>.hetatm.pdb)
    Optionally write PDB files for all unique residue type codes in the HETATM
    lines if  <write_unique>  is True (default True)
    
    OR
    
    Writes all lines in the PDB file  <pdb_filename>  beginning with "HETATM"
    AND with the resName  <only>

    Returns True if successful
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename

    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all HETATM
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()
        good = []
        unique = []
        for i in data:
            resn = i[17:20].strip()
            if i[:6] == 'HETATM' and (not only or resn in only):
                # save for later, more processsing
                good.append( i )
                
                # look for unique resn names
                if not only and not resn in unique:
                    unique.append( resn )
        
        # stop condition
        if not good:
            # tell the user and exit
            print 'No HETATM lines in ' + os.path.relpath( pdb_filename )
            return False
        
        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:        
            if not only:
                out_filename = root_filename + '.hetatm.pdb'
            elif only in WATER_CONVERSION.keys():    # just waters...
                out_filename = root_filename.replace( '.hetatm' , '' ) + '.waters.pdb'
            else:
                # its anything else, name based on the code
                out_filename = root_filename.replace( '.hetatm' , '' ) + '.' + only + '.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()
        
        # change this!
        if not only:
            print 'PDB ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, non-HETATM lines removed\nclean data written to ' + os.path.relpath( out_filename )
        else:
            print 'All ' + only + ' lines in PDB file ' + os.path.relpath( pdb_filename ) + ' written to ' + os.path.relpath( out_filename )

        # optionally redo for all unique members
        if not only and write_unique:
            if len( unique ) > 1:
                # do them all
#                for resn in unique:
#                    clean_HETATM_lines_from_pdb( out_filename , '' , resn )
                unique_filenames = [clean_HETATM_lines_from_pdb( out_filename , '' , resn ) for resn in unique]
                return out_filename , unique_filenames
            else:
                # only 1 HETATM type...
                unique = unique[0]
                print 'Only 1 type of HETATM found, ' + unique

                if unique in WATER_CONVERSION.keys():
                    unique = 'waters'

#                print 'Renaming ' + root_filename + '.hetatm.pdb to ' + root_filename + '.' + unique + '.pdb'
#                shutil.move( root_filename + '.hetatm.pdb' , root_filename + '.' + unique + '.pdb' )
                temp = root_filename + '.' + unique + '.pdb'
                print 'Renaming ' + os.path.relpath( out_filename ) + ' to ' + os.path.relpath( temp )
                shutil.move( out_filename , temp )
                out_filename = temp

        return out_filename
    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

# scan for alternate location fields
def clean_alternate_conformations_from_pdb( pdb_filename , remove_identifier = True ):
    """
    Writes PDB files for each of the alternate conformations found in
    <pdb_filename>
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # verify it exists
    if not os.path.exists( pdb_filename ):
        # for pipelines etc.
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False
    
    # find all alternate conformations
    f = open( pdb_filename , 'r' )
    lines = f.readlines()
    f.close()
    
    # for storage
    non_alternating = ['']
    alternate_conformations = []
    last_line_alternate = False
    index = 0
    alternate_index = -1
    conformation_names = []
    resis = set()
    for i in lines:
        # skip non ATOM lines...fix this later to support header?
        if not i[:6].strip() in ['ATOM' , 'HETATM']:
            last_line_alternate = False
            continue
        
        # sort it
        if i[16].strip():
            conformation = i[16]
            resis.add( i[21] +':'+ i[22:27].strip() )
            
            # optionally remove the alternate conformation identifier
            if remove_identifier:
                i = i[:16] + ' ' + i[17:]
            
            # did we just transition into an alt conf region?
            if last_line_alternate:
                # still in the same region
                if not conformation in alternate_conformations[alternate_index].keys():
                    alternate_conformations[alternate_index][conformation] = i
                    if not conformation in conformation_names:
                        conformation_names.append( conformation )
                else:
                    alternate_conformations[alternate_index][conformation] += i
            else:
                # in a new region
#                if alternate_conformations:
#                    conformation_names = list( set( conformation_names + alternations_conformations[-1].keys() ) )
#                    number_of_conformations = max( number_of_conformations , len( alternate_conformations[-1].keys() ) )
                alternate_index += 1
                alternate_conformations.append( {conformation : i} )
                if not conformation in conformation_names:
                    conformation_names.append( conformation )
            
            last_line_alternate = True
        else:
            # did we just transition into an alt conf region?
            if last_line_alternate:
                # entered a new region
                index += 1
                non_alternating.append( i )
            else:
                # in the same region
                non_alternating[index] += i
                
            last_line_alternate = False
        
    # exit condition
    conformation_names.sort()    # intuitive order...
    if not conformation_names:
        print 'No alternate conformations detected (17th column)'
        return False
    else:
        print 'found ' + str( len( conformation_names ) ) + ' alternate conformations: ' + ', '.join( conformation_names )
        print 'alternate locations found for residues: ' + ', '.join( list( resis ) )
    
#    print index , alternate_index , number_of_conformations
    
    # write out the alternate conformations
    conformation_filenames = []
    for i in conformation_names:
        # make a text by building from fragments
        text = ''
        for j in xrange( len( non_alternating ) - 2 ):
            text += non_alternating[j]
            if i in alternate_conformations[j].keys():
                text += alternate_conformations[j][i]
            else:
                # default to the "first" alt conf ID
                key = 0
                while not conformation_names[key] in alternate_conformations[j].keys():
                    key += 1
                key = conformation_names[key]
                text += alternate_conformations[j][key]

        # add edge case
        text += non_alternating[-1]
        
        # write the file
        out_filename = root_filename + '_conformation_' + i +'.pdb'
        print 'writing conformation ' + i + ' out to ' + os.path.relpath( out_filename ) + ' ...'
        f = open( out_filename , 'w' )
        f.write( text )
        f.close()
        
        conformation_filenames.append( out_filename )

    return conformation_filenames



    

################################################################################
# CONVERTERS

# rewrite the hetatm lines in the pdb file
def convert_pdb_resnames_to_ATOM_lines( hetatm_pdb_filename , out_filename = '' , report_filename = '' , conversion = three2three ):
    """
    Rewrites all HETATM lines in  <hetatm_pdb_filename>  found as keys in
    the dict  <conversion>  and replaces them with their values
    also rewrites the "HETATM" record as "ATOM  "
    
    used to convert HETATM lines that are proxies for amino acids    
    """
    hetatm_pdb_filename = os.path.abspath( hetatm_pdb_filename )
    
    # handle defaults
    if not out_filename:
        # override
        print 'no output filename provided, overwriting ' + hetatm_pdb_filename
        out_filename = hetatm_pdb_filename
    
    # make sure it exists
    if os.path.isfile( hetatm_pdb_filename ):
        # load in the lines
        f = open( hetatm_pdb_filename , 'r' )
        d = f.readlines()
        f.close()
    
        # change to the desired format
        converted = []
        for line in xrange( len( d ) ):
            record = d[line][:6].strip()
            resname = d[line][17:20].strip()
            # go ahead and just rewrite
            if record in ['ATOM' , 'HETATM'] and not resname in one2three.values() and resname in conversion.keys():
                new = conversion[resname]
                d[line] = d[line][:17] + new.rjust(3) + d[line][20:]

                # for records...
                temp = resname + ' lines converted to ' + new
                if not temp in converted:
                    converted.append( temp )
                
                # check the record...all to ATOM
                if record == 'HETATM':
                    d[line] = 'ATOM  '+ d[line][6:]
        
        # debug output
        if converted:
            converted = '\n'.join( converted )
            print converted
            if report_filename:
                print 'summary of converted lines written to ' + report_filename
                f = open( report_filename , 'w' )
                f.write( converted )
                f.close()
        
        # write it back
        f = open( out_filename , 'w' )
        f.writelines( d )
        f.close()
    else:
        print 'No such file named ' + os.path.relpath( hetatm_pdb_filename )
        return False


# useful?

# rewrite the water lines in the pdb file to the standard...from settings?
def convert_water_containing_pdb( hetatm_pdb_filename , conversion = WATER_CONVERSION ):
    """
    Rewrites all HETATM "water" lines in  <hetatm_pdb_filename>  to resNames
    based on  <conversion>
    
    adjust the definition of water (<look_for>) and what to switch to in
    settings.py
    
    not currently used...
    """
    hetatm_pdb_filename = os.path.abspath( hetatm_pdb_filename )
    if os.path.isfile( hetatm_pdb_filename ):
        # load in the lines
        f = open( hetatm_pdb_filename , 'r' )
        d = f.readlines()
        f.close()
    
        # change to the desired format
        for line in xrange( len( d ) ):
            resname = d[line][17:20]
            if resname.strip() in WATER_CONVERSION.keys():
                d[line] = d[line][:17] + WATER_CONVERSION[resname].rjust(3) + d[line][20:]
    
        # write it back...bad!
        f = open( hetatm_pdb_filename , 'w' )
        f.writelines( d )
        f.close()
    else:
        print 'No such file named ' + os.path.relpath( hetatm_pdb_filename )
        return False


# removes lines from  <pdb_file>  and writes to  <out_file>  ending in new
def clean_ATOM_non_new_lines_from_pdb( pdb_filename , out_filename = '' ):
    """
    Write all lines in the PDB file  <pdb_filename>  as long as the last three
    characters on the line aren't "new"
    
    used to clean Hydrogens added using Reduce
    """
    # get the file rootname
    pdb_filename = os.path.abspath( pdb_filename )
    root_filename = get_root_filename( pdb_filename )
    if not root_filename:    # in case it is improper, no "."
        root_filename = pdb_filename
    
    # an optional argument for PDB files not ending in .pdb
#    if not edit:
#        edit = 255

    # if the file exists
    if os.path.exists( pdb_filename ):
        # find all ATOM and TER lines
        f = open( pdb_filename , 'r' )
        data = f.readlines()
        f.close()
        good = []
        for i in data:
            if (i[:5] == 'ATOM ' or i[:4] == 'TER ') and not i.strip()[-3:] == 'new' and i[17:20] in one2three.values():
                good.append( i )

        # stop condition
        if not good:
            # tell the user and exit
            print 'No ATOM non-new lines in ' + os.path.relpath( pdb_filename )
            return False

        # default output file to  <pdb_filename>.clean.pdb
        if not out_filename:
            out_filename = root_filename + '.non_new.pdb'

        # write the found lines
        print 'if the file ' + os.path.relpath( out_filename ) + ' already exists, it will be overwritten!'
        f = open( out_filename , 'w' )
        f.writelines( good )
        f.close()

        print 'PDB file ' + os.path.relpath( pdb_filename ) + ' successfully cleaned, non-ATOM lines lacking \"new\" removed\nclean data written to ' + os.path.relpath( out_filename )
        return out_filename

    else:
        print 'No such file or directory named ' + os.path.relpath( pdb_filename )
        return False

################################################################################
# MAIN

if __name__ == '__main__':
    # parser object for managing input options
    parser = optparse.OptionParser()
    # essential data
    parser.add_option( '-p' , dest = 'pdb_filename' ,
        default = '' ,
        help = 'the pdb filename to process' )
    parser.add_option( '-f' , dest = 'seqformat' ,
        default = SEQFORMAT ,
        help = 'sequence file format, based on settings (!) and Biopython' )
        
    # the other options for the method...are for interactive use
    # hard to manipulate from the commandline...

    (options,args) = parser.parse_args()
    
    # check inputs
    # no edits/modifications
    # kinda silly, but I do this as "my style", easy to modify cleanly
    pdb_filename = options.pdb_filename
    seqformat = options.seqformat

    process_pdb( pdb_filename , seqformat )


################################################################################
################################################################################
# UNFINISHED!!!

# scan for repeated chains and delete them, rewrite it
def clean_redundancy_from_pdb( in_filename , out_filename = '' ):
    """
    Not currently supported
    """
    print 'This is not currently supported sorry...\nIt should look for redundant copies...though it seems the best way to do this is to read directly from the header...but...even for the same sequence, the PDB file may have slightly different coordinates..so how to choose?\nUse process_pdb instead, a separate method is not supported because of this choice problem'

# rewrite a dna or rna pdb to be rosetta friendly
def convert_nucleic_acids_for_rosetta( nucleic_pdb_filename ):
    """
    Not currently supported
    """
    print '...still researching...for whatever reason, most DNA PDB coordinates are accepted in Rosetta (and thus do not crash PyRosetta) however, I cannot get RNA loading to work no matter what (!!??!)\nthey can be ~loaded by modifying the database, but this does not seem to actually do anything, although...make_pose_from_sequence can then make RNA polymers, generating a nonstandard ResidueSet also does not work...perhaps cause the lines are ATOM?'


"""
deprecated stuff...just in case...

#    f = open( pdb_filename , 'r' )
#    complexes = []
#    keep_going = True
#    while keep_going:
        # next line
#        line = f.next()

        # ...think about this...
        # check if chain info, extract the matching subunits
#        if 'CHAIN:' in line:
#            dupl = line.split( 'CHAIN: ' )[-1].replace( ';' , '' ).strip().split( ', ' )    # ignore ";\n"
#            if len( dupl ) > 1:
#                complexes.append( dupl )
        # stop condition
#        elif not ('HEADER' in line or 'TITLE' in line or 'COMPND' in line):
#            keep_going = False        
#    f.close()

    
    # convert complexes
#    if complexes:
#        if not sum( [len( c ) - len( complexes[0] ) for c in complexes] ):
            # all are the same length
#            complexes = [''.join( [c[i] for c in complexes] ) for i in xrange( len( complexes[0] ) )]
#        else:
            # uh oh...
            # could be all should be unique...which puts us in exception land anyway
            # assume that last listed are aberrantly unpaired
#            lowest = min( [len( c ) for c in complexes] )
#            temp = [''.join( [c[i] for c in complexes] ) for i in xrange( lowest )]
#            for c in complexes:
#                temp += c[lowest:]
#            complexes = temp

#                        shutil.copy( seq_filename , seqid2 )
#                        extract_protein_sequence_from_pdb( temp , '' + seqid + '.fa' , seqid , model , chain )
                
                # so...this Biopython object just doesn't work...
                # make a new selection
    #            selection = Select()
    #            selection.accept_model( i )
    #            for l in c:
    #                selection.accept_chain( l )


    # return the filename of the "best" PDB made for Rosetta
    # also return PDB numbering map?
#    if has_nucleic:
#        pdb_filename = root_name + '/' + pdb_filename
#    else:
#        pdb_filename = root_name + '/' + root_name + '.clean.pdb'

    # extract numbering of the best
#    pdb_map , reverse_map = extract_numbering_map_from_pdb( pdb_filename , pdb_filename[:-4] + '_numbering_map.txt' )
    # uh oh, bug...dont wanna fix now
    # until this is proper, leave this to Rosetta...

#    return pdb_filename #, pdb_map



#    if not chain:
#        chain = pdb.child_dict.keys()[0]
    
    # copy the chain

    # conver the model into a Structure, for getting the sequence
#    for_seq = Structure( 'temp' )
    # ...oh yeah...must be model 0 and up >:[
#    temp_model = Model( 0 )
#    for_seq.add( temp_model )
#    for ch in pdb.child_dict.keys():
        # copy it all over directly
#        for_seq[0].add( pdb[ch] )

    # extract the sequence as a Biopython Seq object
    # gap regions makred as "|"
#    seq_builder = PPBuilder()
#    pp = seq_builder.build_peptides( for_seq )
#    seq = Seq( '' , ProteinAlphabet )
#    for frag in pp:
#        seq += frag.get_sequence() + '|'    # already a Biopython Seq

    # or...just do this...

# from making sequences for subdirectories...
#                        temp_seq = SeqRecord( Seq( temp_seq , protein_alphabet ) )
#                        temp_seq.id = os.path.split( root_name )[-1] + ' chain ' + chain
#                        temp_seq.description = '(from model ' + str( model + 1 ) + ')'

"""


