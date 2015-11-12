#!/usr/bin/env python
# :noTabs=true:

"""
these are the basic methods for downloading from UniProt and parsing the record
mainly relies on Biopython for the record parsing
methods setup to parse specific information or handle formatting etc.
UniProt updated changed their record format notably - adding internal URL link
    IDs (recently-ish), I still don't have perfect work arounds, but the basic
    features are still intact

focus on search_uniprot and download_record_from_uniprot
display_uniprot_record_summary can be helpful to find what you want
these methods are weird, designed to be used interactively for searching for
    what you want/testing how to parse it - allowing you to then write a script
    to pull down the information you want

the settings were shamelessly copied from another file, sorry for the mess, some
    of it might be useful, most of it probably won't be

there are several pre-tabulated text files that can make searching faster
    for PDB, STING, SCOP...etc. some exist, some don't, some are very outdated,
    and some exist - but are just plain incorrect
    let me know if you need a specific mapping task that might be desired-enough
    that one of these files might exist

Enjoy!

Author:
    Evan H. Baugh
"""

################################################################################
# IMPORT

# common modules
import re
import urllib2
import operator
import os

# bigger modules
from Bio import SwissProt , ExPASy , SeqIO

# custom modules
#from settings import NON_EXPERIMENTAL_QUALIFIERS , UNIPROT_FEATURE_DETAILS , VARIANT_RECORD_TYPES , UNIPROT_CROSS_REFERENCE_DATABASES_BY_TOPIC

#from helper import get_unique

################################################################################
# SETTINGS

# a bunch of settings setup to be easy to interactively modify
# most of these you won't care about

# simple settings/maps for convenience
UNIPROT_SEARCH_TERMS = [
        'citation' ,
        'clusters' ,
        'comments' ,
        'database' ,
        'domains' ,
        'domain' ,
        'ec' ,
        'id' ,
        'entry name' ,
        'existence' ,
        'families' ,
        'features' ,
        'genes' ,
        'go' ,
        'go-id' ,
        'interpro' ,
        'interactor' ,
        'keywords' ,
        'keyword-id' ,
        'last-modified' ,
        'length' ,
        'organism' ,
        'organism-id' ,
        'pathway' ,
        'protein names' ,
        'reviewed' ,
        'score' ,
        'sequence' ,
        '3d' ,
        'subcellular locations' ,
        'taxon' ,
        'tools' ,
        'version' ,
        'virus hosts' ,
        ]

# only for mapping, updated 9/25/2014
UNIPROT_MAPPING_BY_TOPIC = {
    '2D gel databases' : {
        'World-2DPAGE' : 'WORLD_2DPAGE_ID'
#        'Aarhus/Ghent-2DPAGE' : 'AARHUS_GHENT_2DPAGE_ID' ,    # no longer supported?
#        'ECO2DBASE' : 'ECO2DBASE_ID' ,    # no longer supported?
        } ,
    '3D structure databases' : {
        'DisProt' : 'DISPROT_ID' ,
        'PDB' : 'PDB_ID'
#        'HSSP' : 'HSSP_ID' ,    # no longer supported?
        } ,
    'Chemistry' : {
        'ChEMBL': 'CHEMBL_ID' ,
        'DrugBank' : 'DRUGBANK_ID' ,
        'GuidetoPHARMACOLOGY' : 'GUIDETOPHARMACOLOGY_ID'
        } ,
    'Enzyme and pathway databases' : {
        'BioCyc' : 'BIOCYC_ID' ,
        'Reactome' : 'REACTOME_ID' ,
        'UniPathWay' : 'UNIPATHWAY_ID'
        } ,
    'Gene expression databases' : {
        'CleanEx' : 'CLEANEX_ID'
#        'GermOnline' : 'GERMONLINE_ID'    # no longer supported?
        } ,
    'Genome annotation databases' : {
        'Ensembl' : 'ENSEMBL_ID' ,
        'Ensembl Transcript' : 'ENSEMBL_TRS_ID' ,
        'Ensembl Protein' : 'ENSEMBL_PRO_ID' ,
        'Ensembl Genomes': 'ENSEMBLGENOME_ID' ,
        'Ensembl Genomes Transcript' : 'ENSEMBLGENOME_TRS_ID' ,
        'Ensembl Genomes Protein' : 'ENSEMBLGENOME_PRO_ID' ,
        'GeneID' : 'P_ENTREZGENEID' ,
        'KEGG' : 'KEGG_ID' ,
        'PATRIC' : 'PATRIC_ID' ,
        'UCSC' : 'UCSC_ID' ,
        'VectorBase' : 'VECTORBASE_ID'
#        'GenomeReviews' : 'GENOMEREVIEWS_ID' ,    # no longer supported?
#        'TIGR' : 'TIGR_ID' ,    # no longer supported?
        } ,
    'Organism-specific gene databases': {
#        'AGD' : 'AGD_ID' ,
        'ArachnoServer' : 'ARACHNOSERVER_ID' ,
        'CGD' : 'CGD' ,
        'CYGD' : 'CYGD_ID' ,
        'ConoServer' : 'CONOSERVER_ID' ,
        'dictyBase' : 'DICTYBASE_ID' ,
        'EchoBASE' : 'ECHOBASE_ID' ,
        'EcoGene' : 'ECOGENE_ID' ,
        'euHCVdb' : 'EUHCVDB_ID' ,
        'EuPathDB' : 'EUPATHDB_ID' ,
        'FlyBase' : 'FLYBASE_ID' ,
        'GeneCards' : 'GENECARDS_ID' ,
        'GeneFarm' : 'GENEFARM_ID' ,
        'GenoList' : 'GENOLIST_ID' ,
        'H-InvDB' : 'H_INVDB_ID' ,
        'HGNC' : 'HGNC_ID' ,
        'HPA' : 'HPA_ID' ,
        'LegioList' : 'LEGIOLIST_ID' ,
        'Leproma' : 'LEPROMA_ID' ,
        'MGI' : 'MGI_ID' ,
        'MIM' : 'MIM_ID' ,
        'MaizeGDB' : 'MAIZEGDB_ID' ,
        'neXtProt' : 'NEXTPROT_ID' ,
        'Orphanet' : 'ORPHANET_ID' ,
        'PharmGKB' : 'PHARMGKB_ID' ,
        'PomBase' : 'POMBASE_ID' ,
        'PseudoCAP' : 'PSEUDOCAP_ID' ,
        'RGD' : 'RGD_ID' ,
        'SGD' : 'SGD_ID' ,
        'TAIR' : 'TAIR_ID' ,
        'TubercuList' : 'TUBERCULIST_ID' ,
        'WormBase' : 'WORMBASE_ID' ,
        'WormBase Transcript' : 'WORMBASE_TRS_ID' ,
        'WormBase Protein' : 'WORMBASE_PRO_ID' ,
        'Xenbase' : 'XENBASE_ID' ,
        'ZFIN' : 'ZFIN_ID'
        } ,
    'Other' : {
        'ChiTaRS' : 'CHITARS_ID' ,
        'GeneWiki' : 'GENEWIKI_ID' ,
        'GenomeRNAi' : 'GENOMERNAI_ID' ,
        'NextBio' : 'NEXTBIO_ID'
        } ,
    'Other sequence databases' : {
        'EMBL/GenBank/DDBJ' : 'EMBL_ID' ,
        'EMBL/GenBank/DDBJ CDS' : 'EMBL' ,
        'Entrez Gene (GeneID)' : 'P_ENTREZGENEID' ,
        'GI number*' : 'P_GI' ,
        'PIR' : 'PIR' ,
        'RefSeq Nucleotide' : 'REFSEQ_NT_ID' ,
        'RefSeq Protein' : 'P_REFSEQ_AC' ,
        'UniGene' : 'UNIGENE_ID'
#        'IPI' : 'P_IPI' ,    # no longer supported?
        } ,
    'PTM databases' : {
        'PhosSite' : 'PHOSSITE_ID'
        } ,
    'Phylogenomic databases' : {
        'GeneTree' : 'GENETREE_ID' ,
        'HOGENOM' : 'HOGENOM_ID' ,
        'HOVERGEN' : 'HOVERGEN_ID' ,
        'KO' : 'KO_ID' ,
        'OMA' : 'OMA_ID' ,
        'OrthoDB' : 'ORTHODB_ID' ,
        'ProtClustDB' : 'PROTCLUSTDB_ID' ,
        'TreeFam' : 'TREEFAM_ID' ,
        'eggNOG' : 'EGGNOG_ID'
        } ,
    'Polymorphism databases' : {
        'DMDM' : 'DMDM_ID'
        } ,
    'Protein family/group databases' : {
        'Allergome' : 'ALLERGOME_ID' ,
        'MEROPS' : 'MYCOCLAP_ID' ,
        'PeroxiBase' : 'PEROXIBASE_ID' ,
        'PptaseDB' : 'PPTASEDB_ID' ,
        'REBASE' : 'REBASE_ID' ,
        'TCDB' : 'TCDB_ID'
        } ,
    'Protein-protein interaction databases' : {
        'BioGrid' : 'BIOGRID_ID' ,
        'DIP' : 'DIP_ID' ,
        'MINT' : 'MINT_ID' ,
        'STRING' : 'STRING_ID'
        } ,
    'Protocols and materials databases' : {
        'DNASU' : 'DNASU_ID'
        } ,
    'UniProt': {
        'UniProtKB AC' : 'ACC' ,
        'UniProtKB ID' : 'ID' ,
        'UniProtKB AC/ID' : 'ACC+ID' ,
        'UniParc' : 'UPARC' ,
        'UniRef50' : 'NF50' ,
        'UniRef90' : 'NF90' ,
        'UniRef100' : 'NF100'
         }
    }

# make a raw copy, more list-like
UNIPROT_MAPPING = {}
for i in UNIPROT_MAPPING_BY_TOPIC.keys():
    for j in UNIPROT_MAPPING_BY_TOPIC[i].keys():
        UNIPROT_MAPPING[j] = UNIPROT_MAPPING_BY_TOPIC[i][j]
del i , j


# also some batch processing, updated 9/25/2014
# enables structuring of cross references
UNIPROT_CROSS_REFERENCE_DATABASES_BY_TOPIC = {
    '2D gel databases' : [
        'COMPLUYEAST-2DPAGE',
        'DOSAC-COBS-2DPAGE',
        'OGP',
        'REPRODUCTION-2DPAGE',
        'SWISS-2DPAGE',
        'UCD-2DPAGE',
        'World-2DPAGE'],
    '3D structure databases' : [
        'DisProt',
        'MobiDB',
        'ModBase',
        'PDB',
        'PDBj',
        'PDBsum',
        'ProteinModelPortal',
        'RCSB-PDB',
        'SMR'
        ] ,
    'Chemistry' : [
        'BindingDB' ,
        'ChEMBL' ,
        'DrugBank' ,
        'GuidetoPHARMACOLOGY'
        ] ,
    'Enzyme and pathway databases' : [
        'BioCyc',
        'BRENDA',
        'ENZYME',
        'Reactome',
        'SABIO-RK',
        'SignaLink',
        'UniPathway'
        ] ,
    'Family and domain databases' : [
        'Gene3D',
        'HAMAP',
        'InterPro',
        'PANTHER',
        'Pfam',
        'PIRSF',
        'PRINTS',
        'ProDom',
        'PROSITE',
        'ProtoNet',
        'SMART',
        'SUPFAM',
        'TIGRFAMs'
        ] ,
    'Gene expression databases' : [
        'ArrayExpress' ,
        'Bgee' ,
        'CleanEx' ,
        'Genevestigator' ,
        'ExpressionAtlas'    # not listed, but found in a record...
        ] ,
    'Genome annotation databases' : [
        'Ensembl' ,
        'EnsemblBacteria',
        'EnsemblProtists',
        'EnsemblMetazoa',
        'EnsemblFungi',
        'EnsemblPlants',
        'GeneID',
        'KEGG',
        'PATRIC',
        'UCSC',
        'VectorBase' ,
        'Proteomes'    # wtf!?
        ] ,
    'Ontologies' : [
        'GO'
        ] ,
    'Organism-specific databases' : [
        'ArachnoServer' ,
        'CGD' ,
        'ConoServer' ,
        'CTD' ,
        'CYGD' ,
        'dictyBase' ,
        'EchoBASE' ,
        'EcoGene' ,
        'euHCVdb' ,
        'EuPathDB' ,
        'FlyBase' ,
        'GenAtlas' ,
        'GeneCards' ,
        'GeneFarm' ,
        'GeneReviews' ,
        'GenoList' ,
        'Gramene' ,
        'H-InvDB' ,
        'HGNC' ,
        'HPA' ,
        'HUGE' ,
        'LegioList' ,
        'Leproma' , 
        'MaizeGDB' ,
        'MGI' ,
        'Micado' ,
        'MIM' ,
        'neXtProt' ,
        'Orphanet' ,
        'PharmGKB' ,
        'PomBase' ,
        'PseudoCAP' ,
        'RGD' ,
        'Rouge' ,
        'SGD' ,
        'TAIR' ,
        'TubercuList' ,
        'WormBase' ,
        'Xenbase' ,
        'ZFIN'
        ] ,
    'Other' : [
        'ChiTaRS' ,
        'EvolutionaryTrace' ,
        'GeneWiki' ,
        'GenomeRNAi' ,
        'NextBio' ,
        'PMAP-CutDB' ,
        'PRO' ,
        'SOURCE'
        ] ,
    'PTM databases' : [
        'PhosphoSite' ,
        'PhosSite' ,
        'UniCarbKB'
        ] ,
    'Phylogenomic databases' : [
        'eggNOG' ,
        'GeneTree' ,
        'HOGENOM' ,
        'HOVERGEN' ,
        'InParanoid' ,
        'KO' ,
        'OMA' ,
        'OrthoDB' ,
        'PhylomeDB' ,
        'TreeFam'
        ] ,
    'Polymorphism databases' : [
        'dbSNP' ,
        'DMDM'
        ] ,
    'Protein family/group databases' : [
        'Allergome' ,
        'CAZy' ,
        'GPCRDB' ,
        'IMGT' ,
        'MEROPS' ,
        'mycoCLAP' ,
        'PeroxiBase' ,
        'PptaseDB' ,
        'REBASE' ,
        'TCDB'
        ] ,
    'Protein-protein interaction databases' : [
        'BioGrid' ,
        'DIP' ,
        'IntAct' ,
        'MINT' ,
        'STRING'
        ] ,
    'Proteomic databases' : [
        'MaxQB' ,
        'PaxDb' ,
        'PeptideAtlas' ,
        'PRIDE' ,
        'ProMEX'
        ] ,
    'Protocols and materials databases' : [
        'DNASU' ,
        'SBKB'
        ] ,
    'Sequence databases': [
        'CCDS' ,
        'DDBJ' ,
        'EMBL' ,
        'GenBank' ,
        'PIR' ,
        'RefSeq' ,
        'UniGene'
        ]
    }

UNIPROT_CROSS_REFERENCE_DATABASES = sum( UNIPROT_CROSS_REFERENCE_DATABASES_BY_TOPIC.values() , [] )
UNIPROT_CROSS_REFERENCE_DATABASES = sorted( UNIPROT_CROSS_REFERENCE_DATABASES , key = lambda x : x.lower() )

# more uniprot evidence stuff...
NON_EXPERIMENTAL_QUALIFIERS = {
    'Probable' : 'Best predictions. Indicates stronger evidence than the \"Potential\" qualifier. Implies there is (at least) some experimental evidence indicating the information is expected to be found in the natural environment of a protein.' ,
    'Potential' : 'Better predictions. Indicates there is some logical or conclusive evidence that the given annotation could apply. Often used to present results from protein sequence analysis software tools, which are only annotated if the result makes sense in the biological context of the protein. A typical example is the annotation of N-glycosylation sites in secreted proteins.' ,
    'By similarity' : 'Good predictions. Biological information was experimentally obtained for a given protein (or part of it), and transferred to other protein family members within a certain taxonomic range, dependent on the biological event or characteristic.'
    }

# details on all the UniProt features
UNIPROT_FEATURE_DETAILS = {
    'DOMAIN' : ['domain' , None] ,
    'REPEAT' : ['repeat' , None] ,
    'TRANSIT' : ['localization motif' , None] ,
    'VAR_SEQ' : ['alternate sequence' , None] ,    # use this at all?
    'PEPTIDE' : ['active peptide' , None] ,
    'REGION' : ['region' , None] ,
    'NON_STD' : ['non-standard region' , None] ,
    'SITE' : ['site' , None] ,
    'COMPBIAS' : ['compositional bias' , None] ,
    'PROPEP' : ['propeptide' , None] ,
    'MOTIF' : ['linear motif' , None] ,    # always linear?
    'ACT_SITE' : ['active site' , None] ,
    'BINDING' : ['binding site' , None] ,
    'DNA_BIND' : ['DNA binding site' , None] ,
    'NP_BIND' : ['nucleotide binding site' , None] ,
    'METAL' : ['metal binding site' , None] ,
    'ZN_FING' : ['Zn finger binding' , None] ,    # others for zinc finger? domain?
    'CA_BIND' : ['Ca binding' , None] ,
    'CARBOHYD' : ['glycosylation site' , None] ,
    'LIPID' : ['lipidation site' , None] ,
    'CROSSLNK' : ['crosslink' , None] ,
    'TRANSMEM' : ['transmembrane region' , None] ,
    'INTRAMEM' : ['interior transmembrane region' , None] ,
    'TOPO_DOM' : ['peripheral membrane region' , None] ,

    'MOD_RES' : ['modification site' , None] ,

    'INIT_MET' : ['initial methionine, truncated' , None] ,    # always the case?
    'SIGNAL' : ['signal peptide' , None] ,    # others for signal?
    'DISULFID' : ['disulfide bond' , None] ,

    'CHAIN' : ['chain' , None] ,    # bah!

    # secondary structure
    'HELIX' : ['helix' , None] ,
    'STRAND' : ['strand' , None] ,
    'TURN' : ['turn' , None] ,
    'COILED' : ['coiled coil' , None] ,

    # variations
    'MUTAGEN' : ['mutagenesis variant' , None] ,
    'VARIANT' : ['natural variant' , None] ,
    
    # glitches, problems, etc.
    'CONFLICT' : ['uncertainty, CONFLICT' , None] ,
    'NON_CONS' : ['uncertainty, NON_CONS' , None] ,
    'NON_TER' : ['uncertainty, NON_TER' , None] ,
    'UNSURE' : ['uncertainty' , None]
    }

# should this live here?
# modular way to assign nonsynonymous types
VARIANT_RECORD_TYPES = {
    'missing' : lambda variant , native : 'Missing' in variant ,
    'silent' : lambda variant , native : variant == native ,
    'missense' : lambda variant , native : None ,    # the default
    'nonsense' : lambda variant , native : '*' in variant ,
    'frameshift' : lambda variant , native : len( variant ) > 1 or len( native ) > 1 ,
    'deletion' : lambda variant , native : '-' in variant
    }

################################################################################
# HELPER METHODS

# pfft, I usually need to pull in more than this puny thing...
get_unique = lambda some_list : list( set( some_list ) )

################################################################################
# METHODS

###########
# SWISSPROT
# or UNIPROT

# handle specific database search? comment on it?
# database%3A(type%3Apfam+PF00001)

# customized uniprot search interface...why is Biopython so obfuscated?
def search_uniprot( query , search_options = {} , format = 'list' , columns = [] , sort_by_score = True , include_isoforms = False , compress = False , extra_options = '' ):
    """
    REQUIRES INTERNET CONNECTION

    Returns the results of a UniProt (now contains SwissProt)  <query>
    using  <search_options>  in  <format>  with data as  <columns>
    
    Optionally  <sort_by_score>, internal text querying, HIGHYLY RECOMMENDED!!!
    Optionally  <include_isoforms>  in the results
    Optionally  <compress>  the results (currently unsupported!!!)
    Optionally include  <extra_options>  (for manually specifying details)
    
    to restrict entries, add fields to  <search_options>
    for reviewed entries, add {'reviewed':'yes'}
    for a specific  <organism>, add {'organism':<organism>}
    for a specific database cross-references, try {'database:(type' :
        <database_abbreviation> +'+'+ <entry> +')'}
    
    <search_options>  must be a dict of valid search fields as keys matching
    with values as the matching search term
    <columns>  must be a list of valid search terms

    a list of valid UniProt search terms can be loaded from settings.py
    (UNIPROT_SEARCH_TERMS)

    the preexisting Biopython method, which is included in the cookbook, does
    not appear to work, it seems that the URL it is searching for has changed...
    """
    # url
    url = 'http://www.uniprot.org/uniprot/'
    
    # format the search query
    print 'searching uniprot for \"' + query +'\"'
    if query and search_options:
        query += '+AND+'
    query = '?query=' + query.replace( ' ' , '+' )
    
    # when did this change?
    url += query + '+AND+'.join( [i +':'+ search_options[i] for i in search_options.keys()] )
    if extra_options:
        url += '+AND+' + extra_options
    
    # the format is pretty much ready...
    url += '&format=' + format
    
    # columns
    if columns:
        columns = 'columns=' + ','.join( columns )
        url += '&' + columns
    
    # optionall sort by score
    if sort_by_score:
    	url += '&sort=score'
    
    # optionally include isoforms
    if include_isoforms:
        url += '&include=on'

    # optionally compress...fix this
    if compress:
        url += '&compress=on'

    # go get the results
    print 'obtaining results from:\n\t' + url
    raw_stream = urllib2.urlopen( url )

    # read it
    result = [i for i in raw_stream.read().split( '\n' ) if i]
    print str( len( result ) ) + ' lines found'
    
    return result

# wrapper for uniprot mapping
def search_uniprot_for_mapping( ids , to_id = 'PDB_ID' , from_id = 'ACC' , format = 'tab' , biggest_load = 100 ):
    """
    REQUIRES INTERNET CONNECTION    

    Returns a mapping of  <ids>  (in  <from_id>  format) to a database in
    <to_id>  format with output structured as  <format>
    
    defaults for "tab" format and returns the "mapping" as a list of "pairs"
    [from,to]
    
    handles very large queries by breaking into  <biggest_load>  sized queries
    (prevents error due to very large urls)
    
    available formats for  <to_id>  are listed in UNIPROT_MAPPING (and IDs can
    be found using UNIPROT_MAPPING_BY_TOPIC)
    """
    # make sure list is proper
    # reverse it...split first
    if isinstance( ids , str ):
        ids = ids.split( ' ' )    # should be the delimiter...

    # too many to do at once
    if len( ids ) > biggest_load:
        print 'too many ids to search at once! (the url might get too long)\nsplitting into smaller packets...please be patient...'
        result = {}
        #result = []    # changed to dict
        for i in xrange( len( ids )/biggest_load + 1 ):
            new_ids = ids[i*biggest_load:(i + 1)*biggest_load]
            if new_ids:    # rather than considering the remainder separately, just do it this way :)
                r = search_uniprot_for_mapping( new_ids , to_id , from_id , format )    # leave out the last one...
                
                result.update( r )
                #result += r

        return result    # Im lazy...whatever
        
    # url
    url = 'http://www.uniprot.org/mapping/'
    
    # the query
    url += '?query=' + '+'.join( ids )
    
    # the format...
    url += '&format=' + format
    
    # "from" and "to"
    url += '&from=' + from_id + '&to=' + to_id

    # go get the results
    print 'obtaining results from:\n\t' + url[:50] +'...'*(len( url ) > 50)
    raw_stream = urllib2.urlopen( url )

    # read it
    result = [i for i in raw_stream.read().split( '\n' ) if i]
#    print str( len( result ) ) + ' lines found'
    
    # post-processing
    if format == 'tab':
        result = [i.split( '\t' ) for i in result[1:]]
        
        # conver to dict :)
        result_dict = {}
        for i in result:
            result_dict[i[0]] = i[1:]
        result = result_dict
    
    return result

# simple download wrapper
def download_record_from_uniprot( swissprot_id , out_filename = '' , just_sequence = '' , post_processing = False , fail_filename = 'tmp.sp' , interactive = True ):
    """
    REQUIRES INTERNET CONNECTION

    Downloads and returns the Biopython SwissProt.Record of  <swissprot_id>
    
    SwissProt records contain A LOT more information, currently no there is no
    nicer interface...just learn how to read it form what you are looking for...
    
    Optionally write the content to  <out_filename>
    Optionally return  <just_sequence>  as a Biopython SeqRecord object
        note: if also writing to file,  <just_sequence>  becomes the output ID
    Optionally restructure the record with  <post_processing>
    Optionally provide a small amount of  <interactive>  output, mainly for
        interactive use when searching for a particular protein in a list

    !!! HACK AROUND !!!
    Fall 2014 UniPro changed its API and record
    the records now have "{ECO.*}"...I'm guessing these are for
        smarter hyperlinks?
    well...now I cannot download properly...so I hacked around, filter out these
        strings etc., seems to work although I'm no IO stream expert and since
        these are containers...and not all writable, I just write out to file
        yeah...that's horrible, I need to do this hack around properly...
    """
    # is it a list of ids?
    if isinstance( swissprot_id , list ):
        result = []
        print 'downloading multiple records from UniProt...please be patient...'
        for i in swissprot_id:
            # recursion!
            # NO WRITING!!!
            result.append( download_record_from_uniprot( i , '' , just_sequence , post_processing = post_processing , interactive = False ) )

        # optional interactive sorting
#        if interactive:
#            print    # just to give space
#            for i in xrange( len( result ) ):
#                print i.entry_name.ljust( 20 ) + (','.join( i.accessions )).ljust( 20 ) + i.organism
#                display_uniprot_record_summary_brief( result[i] , front = str( i ) )

    elif isinstance( swissprot_id , str ):
        # should just be an input id
        # get it
        print 'downloading ' + swissprot_id + ' from UniProt...'
        raw_stream = ExPASy.get_sprot_raw( swissprot_id )
        # read it
        if just_sequence:
            # use SeqIO
            result = SeqIO.read( raw_stream , 'swiss' )
        elif fail_filename:    # flag for problems...
            # these "ECO" links (?), try removing these...
            print '!!! hack around !!! check to see if UniProt and Biopython are friends again...'
            
            # get the lines
            result = raw_stream.readlines()
            
            # huh...error page...
            if '<html' in ''.join( result ) and 'error' in ''.join( result ):
                # assume this is enough...
                raise IOError( swissprot_id + ' does not appear to exist!!?! try:\n' + raw_stream.url )
            
            # filter them
#            print result
            result = [''.join( re.split( '\{ECO.*\}' , i ) ) for i in result]
            # write to file...not sure how/if this will handle multiple entries?
            f = open( fail_filename , 'w' )
            f.write( ''.join( result ) )
            f.close()
#            print
#            print result
            
            # load it...requires a file handl...facepalm Biopython...
            f = open( fail_filename , 'r' )
#            result = [i for i in SwissProt.read( f )]
            result = SwissProt.read( f )
            f.close()
            
#            if len( result ) == 1:    # always one?
#                result = result[0]
                
            # cleanup the file...
            os.remove( fail_filename )    # usually unwanted...
        else:
            # ...currently inaccessible, parser cannot handle >:(
            # default to multiple scan
            result = [i for i in SwissProt.parse( raw_stream )]
            
            if len( result ) == 1:
                result = result[0]

    # optional interactive sorting
    if interactive:
        print    # just to give space
#            for i in xrange( len( result ) ):
#                print i.entry_name.ljust( 20 ) + (','.join( i.accessions )).ljust( 20 ) + i.organism
        display_uniprot_record_summary_brief( result )
    
    # optionally write the file
    if out_filename:
        if just_sequence:
            SeqIO.write( result , out_filename , just_sequence )
        else:
            print 'Sorry, the SwissProt record cannot currently be written...I lack the necessay skills to do this :('

    # optionally process the record
    if post_processing:
        # determine the isoforms
        isoforms = determine_isoforms_from_uniprot_record( result , sequences = True )
        if isoforms:
            print str( len( isoforms ) ) + ' isoforms found in the record'
            result.isoforms = isoforms

        # parse the name
        result.names = split_uniprot_name_description( result )

        # parse comments
#        comments = morph_uniprot_record_comments( result )
#        result.comments_str = result.comments    # no...too confusing...
#        result.comments_dict = comments
        result.comments_dict = morph_uniprot_record_comments( result )
        
        # restructure features as a hierarchy
        result.features_dict = morph_uniprot_record_features( result )

        # restructure cross references
        result.cross_references_dict = morph_uniprot_record_cross_references( result )

        ###########
        # summaries        
        # variant summary - compare VARIANT and MUTAGEN
        # structure summary - extract PDB, other dbs...
        # interaction summary - ...for now, this is just copying the 
        
        # pathway summary - check what/any pathways its in
        # PDB interaction - check for directly modeled interactions
    
    return result

#############################
# wrappers for multiple steps

# wrapper for search + map
def search_uniprot_and_map( query , to_id = 'PDB_ID' , search_options = {} ,
        format = 'list' , columns = [] ,
        sort_by_score = True , include_isoforms = False ,
        map_format = 'tab' , biggest_load = 100 ):
    """
    REQUIRES INTERNET CONNECTION    

    Returns a mapping of the search results for  <query>  against the UniProt
    database to a database in  <to_id>  format with output
    structured as  <format>
    
    defaults for "tab" format and returns the "mapping" as a list of "pairs"
    [from,to]
    
    available formats for  <to_id>  are listed in UNIPROT_MAPPING (and IDs can
    be found using UNIPROT_MAPPING_BY_TOPIC)
    """
    # search it
    search = search_uniprot( query , search_options , format , columns , sort_by_score , include_isoforms , False )

    print str( len( search ) ) + ' hits found, mapping them now...'
    
    # get the mao
    mapping = search_uniprot_for_mapping( search , to_id , 'ACC' , map_format , biggest_load )
    
    return mapping

#def search_uniprot( query , search_options = {} , format = 'list' , columns = [] , include_isoforms = False , compress = False ):

# wrapper for search + download
def search_uniprot_and_download( query , out_filename = '' ,
        search_options = {} , format = 'list' , columns = [] ,
        sort_by_score = True , include_isoforms = False , just_sequence = '' , post_processing = True ):
    """
    REQUIRES INTERNET CONNECTION

    Searches UniProt (contains SwissProt now) for  <query>  and downloads
    the results as a Biopython SwissProt.Record (or a list of these objects)
    
    SwissProt records contain A LOT more information, currently no there is no
    nicer interface...just learn how to read it form what you are looking for...

    uses a (limited) custom search of UniProt...details pending...
    """
    # search it
    search = search_uniprot( query , search_options , format , columns , sort_by_score , include_isoforms , False )    # hardcoded, no compression...
    
    print str( len( search ) ) + ' hits found, downloading them now...'
    
    # download
    result = download_record_from_uniprot( search , out_filename , just_sequence , post_processing )
    
    return result

# summarize uniprot
def display_uniprot_record_summary( swissprot_record ,
        show_references = True , show_organism = True ,
        show_sequence = False , show_dates = False , export = False ):
    """
    Displays a summary of the Biopython SwissProt.Record  <swissprot_record>
    
    this is a very silly method intended for interactive use only...
    
    Biopython lacks SwissProt writing capability...and I seriously doubt I could
    do better, so lets trust them and instead offer an interactive option for
    quickly reading SwissProt.Records by eye
    """
    # basic title stuffz
    text = swissprot_record.entry_name +'\n'
    text += swissprot_record.data_class +'\n'
    text += ', '.join( swissprot_record.accessions ) +' (accessions)\n'

    # names...its own section...    
    text += '\n'
    # this is soooo bad...horrible name parsing WTF "human" readable?
    names = split_uniprot_name_description( swissprot_record )
    if 'Recommended Name' in names.keys():
        text += 'Recommended Name:\n'
        if len( names['Recommended Name'] ) == 1:
            text += '\t'+ names['Recommended Name'][0][1] +'\n'
        else:
            for i in names['Recommended Name']:
                text += i[0] +'\t'+ i[1] +'\n'
    if 'Alternate Names' in names.keys():
        text += 'Alternate Names:\n'
        for i in names['Alternate Names']:
            if i[0] == 'Short':
                text += '\t'    # extra tab...
            text += '\t'+ i[1] +'\n'
    # another shameless copy
    if 'Sub Names' in names.keys():
        text += 'Sub Names:\n'
        for i in names['Sub Names']:
            if i[0] == 'Short':
                text += '\t'    # extra tab...
            text += '\t'+ i[1] +'\n'

    # details
    if 'comments_dict' in dir( swissprot_record ):    # already there
        comments = swissprot_record.comments_dict
    else:
        comments = morph_uniprot_record_comments( swissprot_record )
    if comments:
        text += '\n'
    for i in comments.keys():
    	# special cases
#    	if i == 'SIMILARITY':
    	if i == 'similarity':
    	    text += 'SIMILARITY:\n\t' + '\n\t'.join( comments[i] ) +'\n'
#    	elif i == 'INTERACTION':
    	elif i == 'interaction':
    	    text += 'INTERACTIONS:\n\t' + '\n\t'.join( ['\t'.join( [j.ljust( 15 ) , comments[i][j]['name'].ljust( 15 ) , str( comments[i][j]['experiments'] ) , comments[i][j]['IntAct']] ) for j in comments[i].keys()] ) +'\n'
        elif isinstance( comments[i] , str ):
            text += i.upper() +'\t'+ comments[i] +'\n'
        elif isinstance( comments[i] , list ):
            text += i.upper() +'\t'+ ', '.join( comments[i] ) +'\n'
        else:
            # assums its a dict
            text += i.upper() +'\t'+ ', '.join( [j +'='+ comments[i][j] for j in comments[i].keys()] ) +'\n'

    # location...
    if swissprot_record.organelle:
        text += 'organelle: ' + swissprot_record.organelle +'\n'

    # keywords
    if swissprot_record.keywords:
        text += '\nKeywords:\n\t' + '\n\t'.join( swissprot_record.keywords ) +'\n'

    # features!
    # use the morphed version?
    if swissprot_record.features:
        text += '\nFeatures:\n'
        text += '\t'.join( ['Feature' , 'start' , 'end' , 'description'] ) +'\n'
        for i in swissprot_record.features:
            text += '\t'.join( [str( j ) for j in i] ) +'\n'
    
    # references info
    if show_references:
        text += '\n'
        if swissprot_record.references:
            text += str( len( swissprot_record.references ) ) + ' references found, rather than summarizing it/them here,\njust use SwissProt.Record.references to interact with it/them\n(they are SwissProt.Reference objects)\n'
        # uh...is this correct?
        if swissprot_record.cross_references:
            text += 'Other databases:\n'
            for i in swissprot_record.cross_references:
                text += i[0].ljust( 20 ) + '\t'.join( i[1:] ) +'\n'
    
    # origin info
    if show_organism:
        text += '\n'
        if swissprot_record.gene_name:
            text += 'Gene: ' + swissprot_record.gene_name +'\n'
        text += swissprot_record.organism.strip( ' \t.' ) +'\n'
        if swissprot_record.taxonomy_id:
            text += '(NCBI taxid(s) '+ ', '.join( swissprot_record.taxonomy_id ) +')' +'\n'
        if swissprot_record.organism_classification:
            for i in xrange( len( swissprot_record.organism_classification ) ):
                text += '   '*i +'-'+ swissprot_record.organism_classification[i] +'\n'
        # host info
        if swissprot_record.host_organism or swissprot_record.host_taxonomy_id:
            # assume these are parallel...
            if isinstance( swissprot_record.host_organism , list ):
                for i in xrange( len( swissprot_record.host_organism ) ):
                    text += 'Host: ' + swissprot_record.host_organism[i] +'(NCBI taxid '+ swissprot_record.host_taxonomy_id[i] +')\n'
            else:
                # just one? does this even happen?
               text += 'Host: ' + swissprot_record.host_organism +'(NCBI taxid '+ swissprot_record.host_taxonomy_id +')\n'

    # raw sequence stuff...do we care?
    if show_sequence:
        text += '\n'
        if swissprot_record.molecule_type:
            text += swissprot_record.molecule_type +'\n'
        text += str( swissprot_record.seqinfo[0] ) +' residues\t'+ str( swissprot_record.seqinfo[1] ) +' Da\n'
        text += 'CRC32: ' + swissprot_record.seqinfo[2] +'\n'
        text += swissprot_record.sequence +'\n'    # why do ppl ever do this?
    
    # database info
    if show_dates:
        text += '\n'
        text += 'entry created: ' + swissprot_record.created[0] +'\trelease ' + str( swissprot_record.created[1] ) +'\n'
        text += 'sequence updated: ' + swissprot_record.sequence_update[0] +'\trelease ' + str( swissprot_record.sequence_update[1] ) +'\n'
        text += 'annotation(s) updated: ' + swissprot_record.annotation_update[0] +'\trelease ' + str( swissprot_record.annotation_update[1] ) +'\n'

    # oh, and print it
    print text[:-1]
    
    # optionally return it
    if export:
        return text

# simple wrapper
def display_uniprot_record_summary_brief( swissprot_record , front = '' , spread = 20 , too_long = 60 , export = True ):
    """
    Provides a brief summary of  <swissprot_record>  listing only:
        entry name
        (fist) accession
        organism
        additional names
    intended for interactive use, mainly when searching through UniProt's top
    search hits
    
    Optionally add characters in  <front>  of the output
    Optionally  <spread>  the columns to a specified number of characters
    Optionnaly specify how many characters is  <too_long>  (truncate)
    Optionally  <export>  the text generated
    
    <swissprot_record>  can be a Biopython SwissProt.Record object or a list
        of those objects (defaults  <front>  to be the list index)
    """
    if isinstance( swissprot_record , list ):
        text = ''    # make a big str by default, only return if  <export>
        for i in xrange( len( swissprot_record ) ):
            temp = front
            if not front:
                temp = str( i )
            text += display_uniprot_record_summary_brief( swissprot_record[i] , front = temp , spread = spread , export = True )
    
        if export:
            return text
    else:
        # name, accession, organism, length, gene name
#        prepare_text = lambda x : x.ljust( spread )[:spread]

        text = str( front ) + ' '*bool( front )
        text += swissprot_record.entry_name.ljust( spread )
#        text += ( ','.join( swissprot_record.accessions ) ).ljust( spread )[:spread - 1]
        text += ( swissprot_record.accessions[0] + ', ...'*bool( swissprot_record.accessions[1:] ) ).ljust( spread )
#        text += ', '.join( swissprot_record.accessions ).ljust( spread )
        organism_short_name = ' '.join( swissprot_record.organism.split( ' ' )[:2] )    # the first 2 words (?)
        if '(' in swissprot_record.organism:
            organism_short_name = swissprot_record.organism[:swissprot_record.organism.find( '(' )]
        text += organism_short_name.ljust( spread )
#        text += swissprot_record.gene_name.ljust( spread )    # wtf is this?
        text += '\n\t\t'+ ( swissprot_record.organism.strip( '.' )[:too_long] + '...'*bool( len( swissprot_record.organism ) > too_long ) ).ljust( spread )
        # assume this info will be looked up if interested
    
        # deal with alternative names...need that info...
        text += '\n\t'
        names = split_uniprot_name_description( swissprot_record )
        names = sum( [i for i in names.values()] , [] )
        names = [i[-1] for i in names if not i[0] == 'EC']
        text += '\n\t'.join( names )
        
        text += '\n'
        
        print text
        
        # optionally export
        if export:
            return text
    

################################################################################
# RECORD PARSERS

# for specific fields, creating hierarchies, etc.

# simple enough...sorta, extract VAR_SEQ features, generate the sequences
def determine_isoforms_from_uniprot_record( swissprot_record , sequences = True ):
    """
    Returns a dict of the isoforms in  <swissprot_record>  indicating the
    changes from isoform 1

    Optionally create  <sequences>  for the isoforms as well
    
    has this been tested on an insertion isoform?
    could we extract naming details from the record comments? ugh, parsing this
    """
    # extract VAR_SEQ features
    var_seq = [i for i in swissprot_record.features if i[0] == 'VAR_SEQ']
    
    # scane for "isoform" text
    # AAAA!!! handled abstractly...:(
    isoforms = {}
#    print var_seq
    for i in var_seq:
#    	print i
        if 'isoform' in i[-2].lower():
            temp = [j.strip() for j in  i[-2].lower().split( 'isoform' )]
            change = temp[0][::-1]
            change = change[change.find( ' ' ) + 1:][::-1]
            if 'missing' in change:
                change = 'deletion'
            else:
                change = change.upper()
            change = (i[1] , i[2] , change)
            temp = ['isoform ' + j[0] for j in temp[1:]]
            
            # add them
            for j in temp:
                if not j in isoforms.keys():
                    isoforms[j] = {'changes' : [change]}
                else:
                    isoforms[j]['changes'].append( change )
    # messy....extremely messy...
    
    # optionally make altered sequences
    if sequences:
        # default sequences
        for isoform in isoforms.keys():
            isoforms[isoform]['mapped sequence'] = swissprot_record.sequence
    
        # consider substitutions first
        for isoform in isoforms.keys():
            for change in isoforms[isoform]['changes']:
                if not 'deletion' in change:
                    lower = change[0]
                    upper = change[1]
                    temp = change[-1].split( '->' )[1].strip()
                    
#                    print isoforms[isoform]['mapped sequence'][161:165]
                    isoforms[isoform]['mapped sequence'] = isoforms[isoform]['mapped sequence'][:lower - 1] + temp + isoforms[isoform]['mapped sequence'][upper:]
#                    print isoforms[isoform]['mapped sequence'][161:165]
        
        # indels...um...need an example of insertion...
        for isoform in isoforms.keys():
            shifts = []
            for change in isoforms[isoform]['changes']:
                if 'deletion' in change:
                    lower = change[0]
                    upper = change[1]
                    isoforms[isoform]['mapped sequence'] = isoforms[isoform]['mapped sequence'][:lower - 1] + '-'*(upper - lower + 1) + isoforms[isoform]['mapped sequence'][upper:]

        # replace '-' character
        for isoform in isoforms.keys():
            isoforms[isoform]['sequence'] = isoforms[isoform]['mapped sequence'].replace( '-' , '' )

    return isoforms

# wrapper for the crazy naming!
def split_uniprot_name_description( swissprot_record ):
    """
    Returns a dict of lists describing the names in  <swissprot_record>
    
    <swissprot_record>  should be a Biopython Swissprot.Record object,
    specifically this function operates on Swissprot.Record.description
    
    ...probably needs some work...
    """
    # this is soooo bad...horrible name parsing WTF "human" readable?
    names = {}
    
    # hacky...
    description = swissprot_record.description.replace( 'Flags: Precursor;' , '' ).replace( '; Contains: ' , '' )
    
    if 'RecName: ' in description:
        names['Recommended Name'] = []
        temp = [i for i in description.split( 'RecName: ' ) if i][0]    # first non-empty
        if 'AltName: ' in description:
            temp = temp[:temp.find( 'AltName: ' )]    # only until here
        temp = [i.strip() for i in temp.split( ';' ) if i.strip()]
        for i in temp:
            temp2 = i.split( '=' )
            if len( temp2 ) > 1:
                names['Recommended Name'].append( [temp2[0] , temp2[1].strip()] )
            else:
                names['Recommended Name'].append( [temp2[0] , '???'] )                
    if 'AltName: ' in description:
        names['Alternate Names'] = []
        temp = description.split( 'AltName: ' )[1:]
        # over each alt name
        for i in temp:
            if 'Full' in i:
                temp2 = i.split( 'Full=' )[1]
                names['Alternate Names'].append( ['Full' , temp2[:temp2.find( ';' )].strip()] )
            temp2 = [j.strip() for j in i.split( ';' ) if j.strip() and not j[:5] == 'Full=' and not 'Flags:' in j]
            for j in temp2:
                temp3 = j.split( '=' )
#                print temp3
                if len( temp3 ) == 1:    # super hacky "fix"
                    temp3.append( '' )
                names['Alternate Names'].append( [temp3[0] , temp3[1]] )
    # shamelessly copied
    if 'SubName: ' in description:
        names['Sub Names'] = []
        temp = description.split( 'SubName: ' )[1:]
        # over each alt name
        for i in temp:
            if 'Full' in i:
                temp2 = i.split( 'Full=' )[1]
                names['Sub Names'].append( ['Full' , temp2[:temp2.find( ';' )].strip()] )
            temp2 = [j.strip() for j in i.split( ';' ) if j.strip() and not j[:5] == 'Full=' and not 'Flags:' in j]
            for j in temp2:
                temp3 = j.split( '=' )
                names['Sub Names'].append( [temp3[0] , temp3[1]] )

    return names

# silly wrapper
def morph_uniprot_record_comments( swissprot_record ):
    """
    Returns a dict summarizing the comments section of  <swissprot_record>
    
    ...just read the function if you are interested in specific details etc.

    INTERACTION comments are parsed into a dict with UniProt ID keys
    
    recently reworked it to have lower case keys...this may conflict with other
        methods
    """
    # pretty easy...
    comments = {}
    for i in swissprot_record.comments:
        temp = i.split( ':' )
        key = temp[0].lower()
        if key in comments.keys():
            # special cases...
#            if key == 'CATALYTIC ACTIVITY':
            if key == 'catalytic activity':
                # list, single entry per...
                comments[key].append( ''.join( temp[1:] ).strip( ' \t.;' ) )
#            elif key == 'PATHWAY':
            elif key == 'pathway':
                # list, potentially multiple...
                comments[key] += [j.strip() for j in ''.join( temp[1:] ).strip( ' \t.;' ).split( ';' )]
#            elif key == 'BIOPHYSICHEMICAL PROPERTIES':
            elif key == 'biophysichemical properties':
                # list, potentially multiple...
                comments[key] += [j.strip() for j in ''.join( temp[1:] ).strip( ' \t.;' ).split( ';' )]
            else:
                comments[key] += ', '+ ''.join( temp[1:] ).strip( ' \t.;' )
#            raise IOError( 'fix it! uniprot morph_uniprot_record_comments...' )
        else:
            # special cases...
#            if temp[0] == 'INTERACTION':
            if key == 'interaction':
                # split it first
                interactions = {}
                current_key = ''
                comments[key] = ' '.join( temp[1:] ).strip( ' \t.;' )    # lazy, do this first
                for j in comments[key].split( '; ' ):
                    if j[:6] == 'NbExp=':
                        # number of experiments
                        interactions[current_key]['experiments'] = int( j[6:] )
                    elif j[:7] == 'IntAct=':
                        # IntAct database ID(s)
                        interactions[current_key]['IntAct'] = j[7:]
                    else:
                        # should be a new key
                        temp2 = j.split( ' ' )
                        current_key = temp2[0]

                        # silly conditional...this is so ungly, plz rework soon
                        if not len( temp2 ) > 1:
                            temp2 = ''
                        else:
                            temp2 = temp2[1]
                        
                        interactions[current_key] = {'name' : temp2}
                comments[key] = interactions

            # have "append" issues, when do they start?
#            elif key == 'CATALYTIC ACTIVITY':
            elif key == 'catalytic activity':
                comments[key] = [''.join( temp[1:] ).strip( ' \t.;' )]
#            elif key == 'PATHWAY':
            elif key == 'pathway':
                comments[key] = [j.strip() for j in ''.join( temp[1:] ).strip( ' \t.;' ).split( ';' )]
#            elif key == 'BIOPHYSICHEMICAL PROPERTIES':
            elif key == ' biophysichemical properties':
                comments[key] = [j.strip() for j in ''.join( temp[1:] ).strip( ' \t.;' ).split( ';' )]
            else:
                comments[key] = ' '.join( temp[1:] ).strip( ' \t.;' )

    # super special cases
#    if 'SIMILARITY' in comments.keys():
#        comments['SIMILARITY'] = [j.strip() for j in comments['SIMILARITY'].split( ',' )]
    if 'similarity' in comments.keys():
        comments['similarity'] = [j.strip() for j in comments['similarity'].split( ',' )]
    if 'web resource' in comments.keys() and isinstance( comments['web resource'] , str ):
        comments['web resource'] = dict( [(j[:j.find( '=' )].strip( ' \'\"' ) , j[j.find( '=' ) + 1:].strip( ' \'\"' ).replace( 'http ' , 'http:' ) ) for j in comments['web resource'].split( ';' )] )
    if 'subcellular location' in comments.keys() and isinstance( comments['subcellular location'] , str ):
        comments['subcellular location'] = [j.strip() for j in comments['subcellular location'].split( '.' ) if j.strip()]
    # "alternate products" is too messy to go after now...
    
    return comments

# used in 2 other places, leaving it for now although the "morph" version is clearer
def extract_uniprot_record_features( swissprot_record ):
    """
    Returns a dict summarizing the features of  <swissprot_record>
    with keys as feature types and values as lists of dicts representing the
    'start' , 'end', and 'contents' of each feature
    
    ...deprecated but maintained due to dependency...
    should fix that...
    """
    feature_dict = {}
    for i in swissprot_record.features:
        feature_type = i[0]

        # try to make int if possible
        start = i[1]
        if isinstance( start , str ) and start.replace( '-' , '' ).replace( '.' , '' ).isdigit():
            start = int( start )
        end = i[2]
        if isinstance( end , str ) and end.replace( '-' , '' ).replace( '.' , '' ).isdigit():
            end = int( end )

        contents = {
            'start' : start ,
            'end' : end ,
            'contents' : [j for j in i[3:] if j]
            }
        
        # add the feature if it does not exist
        if feature_type in feature_dict.keys():
            feature_dict[feature_type].append( contents )
        else:
            feature_dict[feature_type] = [contents]
    
    # parse variants
    for i in ['VARIANT' , 'MUTAGEN']:
        if i in feature_dict.keys():
            bonus = []
            for j in xrange( len( feature_dict[i] ) ):
                #print feature_dict[i][j]
                summary = extract_uniprot_variation_details( feature_dict[i][j]['contents'][0] )
                # how to handle multiple entries? as multiple
                if isinstance( summary , list ):
                    bonus += [(k , feature_dict[i][j]['start'] , feature_dict[i][j]['end'] , feature_dict[i][j]['contents']) for k in summary[1:]]
                    summary = summary[0]
                feature_dict[i][j]['native'] = summary['native']
                feature_dict[i][j]['variant'] = summary['variant']
                feature_dict[i][j]['brief'] = summary['native'] + str( feature_dict[i][j]['start'] ) + summary['variant']
                feature_dict[i][j]['type'] = summary['type']

            # found extra
            for j in bonus:
                feature_dict[i].append( {} )
                feature_dict[i][-1]['start'] = j[1]
                feature_dict[i][-1]['end'] = j[2]
                feature_dict[i][-1]['contents'] = j[3]

                feature_dict[i][-1]['native'] = j[0]['native']
                feature_dict[i][-1]['variant'] = j[0]['variant']
                feature_dict[i][-1]['brief'] = j[0]['native'] + str( j[1] ) + j[0]['variant']
                feature_dict[i][-1]['type'] = j[0]['type']

#                summary['contents'] = feature_dict[i][j]['contents']
#                feature_dict[i][j] = summary
                
                # classify it
#                if '*' in var:
#                    feature_dict[i][j]['type'] = 'nonsense'
#                elif 'Missing' == var:
#                    feature_dict[i][j]['type'] = 'Missing'
#                elif len( wt ) > 1 or len( var ) > 1:
#                    feature_dict[i][j]['type'] = '? frameshift ?'
#                elif wt == var:
#                    feature_dict[i][j]['type'] = 'silent'
#                else:
#                    feature_dict[i][j]['type'] = 'missense'

    return feature_dict

# silly wrapper...
# wow, gone through several versions, wish this was easier
def morph_uniprot_record_features( swissprot_record , non_experimental_qualifiers = NON_EXPERIMENTAL_QUALIFIERS , uniprot_feature_details = UNIPROT_FEATURE_DETAILS , variant_features = ['natural variant' , 'mutagenesis variant'] ):
    """
    Returns a dict summarizing the features of  <swissprot_record>
    with keys as feature types and values as lists of dicts representing the
    'start' , 'end', and 'contents' of each feature
    """
    feature_dict = {}
    for i in swissprot_record.features:
        i = list( i )    # for assignment later
        feature_type = i[0]#.strip( ' .:;' )lower()
#        feature_type = uniprot_feature_details[i[0]][0]

        # determine the qualifier...if any, and clean the text
        qualifier = None
        for j in xrange( len( i ) ):
            for k in non_experimental_qualifiers.keys():
                if isinstance( i[j] , str ):
                    if k in i[j]:    # found one
                        if qualifier and not qualifier == k:    # found another one...uh oh...
                            print i
                            raise IOError( '??? found a case with multiple qualifiers' )
                        qualifier = k
#                        break

                    # remove the qualifiers...these are just unpleasant
                    i[j] = i[j].replace( k , '' ).replace( '()' , '' ).replace( '[]' , '' ).strip( ' .:;' )
        
        # try to make int if possible
        start = i[1]
        if isinstance( start , str ) and start.replace( '-' , '' ).replace( '.' , '' ).replace( '<' , '' ).replace( '>' , '' ).isdigit():
            start = int( start.replace( '<' , '' ).replace( '>' , '' ) )
        end = i[2]
#        print start , end
        if isinstance( end , str ) and end.replace( '-' , '' ).replace( '.' , '' ).replace( '<' , '' ).replace( '>' , '' ).isdigit():
            end = int( end.replace( '<' , '' ).replace( '>' , '' ) )
        
        # determine the position
        position = None
        if start and end:
            if start == end:
                position = start
            else:
                # a range
                # make a special exception for disulfides?
                position = [start , end]

        # extra details
        contents = [j.strip( ' .:;' ) for j in i[3:] if j.strip( ' .:;' )]
        subtype = None
        if contents:
            if not uniprot_feature_details[feature_type][0] in variant_features:
                subtype = contents[0] +' '+ uniprot_feature_details[feature_type][0]
#                if len( contents ) > 2:
                contents = contents[1:]    # no need for reduplication
#                elif not contents:
#                    contents = None
                    # hopefully this isn't too confusing...
                    # actually, it is, leave it as is...
#                else:
#                    contents = contents[1]
            else:
                # special handling for variants
                if isinstance( position , list ):    # no need to list it twice
#                    print i
#                    raise IOError( 'more than 1 position listed!?! this should never happen!!!' )
                    print i[3] + ', more than 1 position listed!?! this should never happen!!!'
                if qualifier:
                    print i
                    raise IOError( 'how is there a non-experimental qualifier here!!?!' )
                
                qualifier = extract_uniprot_variation_details( contents[0] )
                # include all contents, lots of hidden details poorly incorporated here as text...

#                if len( contents ) > 1:
#                    print i
#                    raise IOError( 'not currently supported :/ fix that!' )

        # store as a list
        contents = [subtype , position , qualifier , contents]
        # a dict instead? annoying to reference 'position' instead of index?
        
        # add the feature if it does not exist
        feature_type = uniprot_feature_details[feature_type][0]    # change at the last second...
        if feature_type in feature_dict.keys():
            feature_dict[feature_type].append( contents )
        else:
            feature_dict[feature_type] = [contents]

    # convert to a dict
    # must know how many of each feature, for numbering etc.
    for i in feature_dict.keys():
        # determine if subtypes are unique or need to be renamed
        subtypes = [j[0] for j in feature_dict[i]]
        if len( subtypes ) == len( get_unique( subtypes ) ):
            # all are unique, just convert into a dict :)
            feature_dict[i] = dict( [(j[0] , j[1:]) for j in feature_dict[i]] )

        elif i in variant_features:
            # special handling...
            # all should have "subtype" (0) = None
            # "qualifier" (2) should contain the parsed variant info
            variants = {}
            trouble = {}
            for j in feature_dict[i]:
                if isinstance( j[2] , list ):
                    # multiple variants
                    for k in j[2]:
                        key = k['native'] + str( j[1] ) + k['variant']

                        if key in variants.keys():
#                            print key , variants.keys()
#                            print feature_dict[i]
#                            raise IOError( 'multiple records of the this variant!!?' )
                            print 'multiple records of ' + key + ' variant!!? requires manual check...'
#                            key = 'aberrant_' + key
#                            print variants[key]
#                            print j[2] in variants[key][1]
                            if not j[-1] in variants[key][-1]:
                                variants[key][-1] += j[-1]
#                            if not j[2] in variants[key][1]:    # don't actually need this
#                                variants[key][1].append( j[2] )

                        variants[key] = j[1:]
                else:
                    # common case, derive the "brief" and use as subtype
                    key = j[2]['native'] + str( j[1] ) + j[2]['variant']

                    if key in variants.keys():
#                        print key
#                        print feature_dict[i]
#                        raise IOError( 'multiple records of the this variant!!?' )
                        print 'multiple records of ' + key + ' variant!!? requires manual check...'
#                        key = 'aberrant_' + key
#                        print variants[key]
#                        print j[2] in variants[key][1]
                        if not j[-1] in variants[key][-1]:
                            variants[key][-1] += j[-1]
#                        if not j[2] in variants[key][1]:    # don't actually need this
#                            variants[key][1].append( j[2] )

                    variants[key] = j[1:]
            feature_dict[i] = variants

        else:
            # determine names for unamed portions
            # default to renaming ALL entries in this situation...
            feature_dict[i] = dict( [( i +' '+ str( j + 1 ) , feature_dict[i][j][1:-1] + [ [feature_dict[i][j][0]] + feature_dict[i][j][-1] ]) for j in xrange( len( feature_dict[i] ) )] )
            # better way to join?
            # no, should not merge entries of the same feature name, may be
            #    incorrect with multiple domains etc. e.g. metal binding

    # third stage of processing...
    # for variants, consider possible overlap between VARIANT and MUTAGEN
#    for i in variant_features if i in feature_dict.keys()
    # NO! save this for later, this is already getting too layered with
    #    personal/interactive-favored parsing :[
    # how about some extra cleanup too?
    # NO!...at least not what was considered with that comment
    #    leave is as [position , qualifier , [contents]]
    
    # need to consider multiple multiple-mutations...very messy, but can be automated...
            
    return feature_dict
    
# simple text parsing helper for feature extraction
# ...there is a lot more going on here that it seems...
def extract_uniprot_variation_details( variant_text , position = False , variant_record_type = VARIANT_RECORD_TYPES ):
    # lots of splitting
    temp = [i.strip() for i in variant_text.replace( '->' , ' ' ).split( ' ' ) if i.strip()]
#    print temp[0] , temp[1]
    
    summary = {}
    
    if 'Missing' in temp[0]:
        summary['native'] = 'missing'
        summary['variant'] = 'missing'
        summary['type'] = 'missing'
    else:
        native = temp[0].strip()
        variant = temp[1].strip( ' .:;' )
        
        # how to handle multiple vars?
        # do it here...
        if ',' in variant:
            variant = [i.strip( ' .,:;' ) for i in variant.split( ',' )]
            
            # now multiple entries
            summary = []
            for i in variant:
                variant_type = [j for j in variant_record_type.keys() if variant_record_type[j]( i , native )]
                if len( variant_type ) > 1:
#                    print variant_text
#                    raise IOError( 'multiple types!? need to revise how variants are parsed!' )
                    print variant_text + ', multiple types!? need to revise how variants are parsed!'
                    variant_type = '?'
                elif len( variant_type ) == 1:
                    variant_type = variant_type[0]    # only one!
                else:
                    variant_type = 'missense'

                summary.append( {
                    'native' : native ,
                    'variant' : i ,
                    'type' : variant_type
                    } )
        else:
            # just one variant
            variant_type = [j for j in variant_record_type.keys() if variant_record_type[j]( variant , native )]
            if len( variant_type ) > 1:
                print variant_text
                raise IOError( 'multiple types!? need to revise how variants are parsed!' )
            elif len( variant_type ) == 1:
                variant_type = variant_type[0]    # only one!
            else:
                variant_type = 'missense'

            summary = {
                'native' : native ,
                'variant' : variant ,
                'type' : variant_type
                }
    
#    return wt , var
    return summary

# abstract hierarchy, use the mapping topics
def morph_uniprot_record_cross_references( swissprot_record , uniprot_cross_references_by_topic = UNIPROT_CROSS_REFERENCE_DATABASES_BY_TOPIC ):
    # ...faster right? create an inverse map
    database2topic = dict( sum( [[( j , i ) for j in uniprot_cross_references_by_topic[i]] for i in uniprot_cross_references_by_topic.keys()] , [] ) )

    cross_references = {}
    # look for categories
    for i in swissprot_record.cross_references:
        # debug
        if not i[0] in database2topic:
            # can happen spontaneously...
            # rather than make robust here, interactively add to the setting file
            raise IOError( '!!?! ' + i[0] + ' is not in the list of databases!!?!' )

        category = database2topic[i[0]].lower()    # I prefer lower case
        
        if category in cross_references:
            if i[0] in cross_references[category].keys():
                cross_references[category][i[0]].append( i[1:] )
            else:
                # this database has not been seen yet
                cross_references[category][i[0]] = [i[1:]]
        else:
            # not even the category has been seen
            cross_references[category] = {
                i[0] : [i[1:]]
                }

    return cross_references

###########
# summaries

# beyond simple parsing, may involve comparing between information

        # variant summary - compare VARIANT and MUTAGEN
        # structure summary - extract PDB, other dbs...
        # interaction summary - ...for now, this is just copying the 
        
        # pathway summary - check what/any pathways its in
        # PDB interaction - check for directly modeled interactions

# simple enough...
def generate_variant_summary_from_uniprot_record( swissprot_record ):
    # check for feature_dict
    if not 'features_dict' in dir( swissprot_record ):
        # should reduplicate all of the options...:/
        swissprot_record.features_dict = morph_uniprot_record_features( swissprot_record )
    
    # necessary at all?
    if 'natural variant' in swissprot_record.features_dict.keys() and 'mutagenesis variant' in swissprot_record.features_dict.keys():
        # actually need to compare them
        # start with natural variants
        variants = {}
        variants.update( swissprot_record.features_dict['natural variant'] )
        # how to call these? a single character make more sense?
        for i in variants.keys():
#            variants[i][2] = ['natural'] + variants[i][2]
            variants[i][2] = ['N'] + variants[i][2]

        for i in swissprot_record.features_dict['mutagenesis variant'].keys():
            if i in variants.keys():
                # check that they are the same otherwise?
#                variants[i][2] = ['natural and mutagenesis'] + variants[i][2][1:] + swissprot_record.features_dict['mutagenesis variant'][2]
                variants[i][2] = ['NM'] + variants[i][2][1:] + swissprot_record.features_dict['mutagenesis variant'][2]
            else:
                # add it
                variants[i] = swissprot_record.features_dict['mutagenesis variant'][i]
#                variants[i][2] = ['mutagenesis'] + variants[i][2]
                variants[i][2] = ['M'] + variants[i][2]
    elif 'natural variant' in swissprot_record.features_dict.keys():
        # just send the natural variants
        variants = swissprot_record.features_dict['natural variant']
    elif 'mutagenesis variant' in swissprot_record.features_dict.keys():
        # just send the mutagenesis variants
        variants = swissprot_record.features_dict['mutagenesis variant']
    else:
        # neither...send empty?
        variants = {}
    
    return variants

# get the PDBs as a dict
# older than the summaries idea, but performs a very similar function
def extract_pdb_codes_from_uniprot_record_cross_references( swissprot_record ,
        filter_and_sort = True , experimental_types = ['X-ray'] ,
        exclude_missing = True , include_other_structures = False ):
    """
    Returns a dict of the PDB codes found in the cross references of
    <swissprot_record>  with layers of experimental type (X-ray vs NMR)
    and the PDB code followed by the resolution (X-ray only, in Angstroms),
    relevant chain IDs, and matching lists of start and stop positions
    
    Optionally  <filter_and_sort>  counting  <experimental_types>  and
        <exclude_missing>  (put into 'sorted')
    Optionally  <include_other_structures>  (accessions from several other
        structure databases...pretty much all of the UniProt "3D structure
        databases" category of cross-reference databases)

    
    what about Gene3D?
    what about NESG? included on PMP...but annoying to hack
   
    so...SMR stands for SwissModelRepository...great! since I can remotely
    query them too...oh...it only refers to attempts...not models...useless
    ModBase searching is integrated into the UniProt website...but not
    in the record...lame...
    
    Note: some of this information may be missing or outdated for numerous
        reasons
    Note: the "matched" regions may not be accurate as modifications are made
        to crystallize the protein, residues may be missing, and the numbering
        scheme may also differ
    """
    # get the PDB cross references
    pdbs = [i for i in swissprot_record.cross_references if i[0] == 'PDB']
    # 0     "PDB"               useless
    # 1     PDB code            which it is
    # 2     method              X-ray or NMR...idk of any others...
    # 3     quality             for X-ray, ignore the " A" at the end
    # 4     chains=start-stop   ...putative...not always accurate :[
#    pmp = [i[1:] for i in swissprot_record.cross_references if i[0] == 'ProteinModelPortal']
    # should make this more modular, take the "3D structure databases" category of UniProt cross-reference databases?
#    other_structures = [i[1:] for i in swissprot_record.cross_references if i[0] in ['ProteinModelPortal' , 'SMR']]
    other_structures = [i for i in swissprot_record.cross_references if i[0] in ['DisProt' , 'MobiDB' , 'ModBase' , 'ProteinModelPortal' , 'SMR']]
    # indicates other models exist, but also that a candidate structure
    # likely exists
    
    
    # sort by experiment, code
    pdb_dict = {}
    for i in pdbs:
        code = i[1]
        experiment = i[2]
        quality_and_chain = i[3]
    
        # add entry to the dict
        if not experiment in pdb_dict.keys():
            pdb_dict[experiment] = {}
        
        # assess the quality and the chains + range info
        quality = ''
        if quality_and_chain and not quality_and_chain == '-':
            quality = float( quality_and_chain[:-2] )    # ignore the " A"

        # chains + range info...tricky...
        chains = []
        chains = {}
#        start = []
#        stop = []
        if i[4] and not i[4] == '-':
            # split by ", ": multiple spans
            for j in i[4].split( ', ' ):
                # split by "=": chains vs range
                temp = j.split( '=' )
                # split by "/": mutiple chains
                new_chains = temp[0].split( '/' )
                # split by "-": start-stop
                temp = temp[1].split( '-' )
                new_start = None
                if temp[0]:
                    new_start = int( temp[0] )
                new_stop = None
                if temp[1]:
                    new_stop = int( temp[1] )
                
                # add them
                for k in new_chains:
#                    chains.append( k )
#                    start.append( new_start )
#                    stop.append( new_stop )
                    chains[k] = [new_start , new_stop]

        # assign it
#        pdb_dict[experiment][code] = [
#            quality ,
#            chains ,
#            start , stop
#            ]
        pdb_dict[experiment][code] = {
            'chains' : chains
            }
        if quality:
            pdb_dict[experiment][code]['resolution'] = quality

    # optionally filter and sort the codes
    if filter_and_sort and pdb_dict:
        # filter by experimental types
        pdbs = {}
        for i in experimental_types:
            if i in pdb_dict.keys():
                pdbs.update( pdb_dict[i] )
    
        # filter by missing data
        if exclude_missing:
            for i in pdbs.keys():
#                if not pdbs[i][0] or not pdbs[i][1] or not pdbs[i][2] or not pdbs[i][3]:
                if not 'resolution' in pdbs[i].keys() or not pdbs[i]['resolution'] or not 'chains' in pdbs[i].keys() or not pdbs[i]['chains'] or [None for j in pdbs[i]['chains'].values() if not j[0] or not j[1]]:
                    del pdbs[i]
        
        # make and sort the list
        pdbs = sorted( pdbs.items() , key = operator.itemgetter( 1 , 0 ) )
        # descending quality, already "best" is the first one
        
        # ah, old way, instead, just add it
#        pdb_dict = pdbs
        pdb_dict['sorted'] = pdbs
        
    # optionally include PMP, SMR...no modbase :(
    if include_other_structures or not pdb_dict:
        if not pdb_dict:
            print 'no exact matches found...automatically returning results from other structure databases'
#        return pdb_dict , other_structures
        # perhaps this is too deeply used...
#        if isinstance( pdb_dict , dict ):
        if other_structures:
            pdb_dict['other structures'] = other_structures
#        elif isinstance( pdb_dict , list ):
#            pdb_dict.append( other_structures )
        
    # what about Gene3D ?
    
    return pdb_dict

# for quick title scanning...not really anything more
def display_uniprot_reference_summary( swissprot_reference , export = False ):
    """
    Displays the contents of a Bio.Swissprot.Reference object in a clean way
    should really be added as its __str__ method
    
    Optionally  <export>  the summary text
    """
    # robust to input
    if isinstance( swissprot_reference , list ):
        # assume its a list of references
        text = [display_uniprot_reference_summary( i , export ) for i in swissprot_reference]
        if export:
            return text
    elif not 'authors' in dir( swissprot_reference ) and 'references' in dir( swissprot_reference ):
        # um...it could be just the record itself...
        # recusion as list mode
        text = display_uniprot_reference_summary( swissprot_reference.references )
        if export:
            return text
    else:
        # normal treatment for a record
        text = ''
    
        # essential info
#    if 'title' in dir(    # no, these should all be here
        text += swissprot_reference.title +'\n'
        text += swissprot_reference.authors +'\n'
        text += swissprot_reference.location +'\n\n'
    
        # additional references
        for i in swissprot_reference.references:
            text += ' '.join( i ) +'\n' 
    
        if swissprot_reference.comments:
            text += '\n'
            for i in swissprot_reference.comments:
                text += ' '.join( i ) +'\n'
    
        # just display it
        print text[:-1]

        # optionally export
        if export:
            return text[:-1]

# crude method for comparing taxonomic similarity
def determine_taxonomic_distance_proxy( swissprot_record1 , swissprot_record2 , similarity = False ):
    """
    Returns the "taxonomic distance" between  <swissprot_record1>  and
    <swissprot_record2>  defined as the number of nodes* on the taxonomic tree
    that differ between the organisms these records (proteins) come from
    
    smaller distances mean more similar organisms
    
    Optionally return the  <similarity>  instead of the distance
    
    ex. a distance of 4 corresponds to at least 4 terminal nodes that differ
    between these organisms
    
    *technically, it counts the number of equivalent nodes and reports the
    first discrepancy, actually returning the difference between the smaller
    number of nodes and this similarity
    while similarity may seem more appropriate, it is not reasonably bounded
    """
    # check the records
    if isinstance( swissprot_record1 , str ):
        swissprot_record1 = download_record_from_uniprot( swissprot_record1 )
    if isinstance( swissprot_record2 , str ):
        swissprot_record2 = download_record_from_uniprot( swissprot_record2 )

    smaller = min( len( swissprot_record1.organism_classification ) , len( swissprot_record2.organism_classification ) )
    for i in xrange( smaller ):
        if not swissprot_record1.organism_classification[i] == swissprot_record2.organism_classification[i]:
            break
        
    # i is now the "similarity"
    if similarity:
        return i
    
    # just use the smaller?
    distance = smaller - i
    return distance

################################################################################
# LABELS

# user/author may disagree about annotations
# simple methods for evaluating properties

def is_membrane_protein( swissprot_record , strict = False ):
    if not strict:
        tm = bool( [i for i in swissprot_record.features if i[0] == 'TRANSMEM' or i[0] == 'LIPID'] )
    else:
        # even more strict? no "membrane" association?
        tm = bool( [i for i in swissprot_record.keywords if 'membrane' in i.lower()] )
    
    return tm

