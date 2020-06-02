## KEGG Interaction Parsing
## Anna Ritz, Feb 2013

## IMPORT STATEMENTS
import sys
import argparse
import collections
#from restful_lib import Connection
import requests ## udpate 2020
from xml.dom.minidom import *
from utilsPoirel import *
from rich_utils import *
import re
import os
from time import strftime

## STATIC VARIABLES

# Step 1: Connect to Database
KEGG_URL = "http://rest.kegg.jp"
#CONN = Connection(KEGG_URL)
NULL_ENTRY = 'None'
NAMES = {}
NAMELIST = []
IDSYMBOL = {"uniprot" : "up" , "hsa" : "hsa", "pubchem" : "pubchem"}

# get Kegg id to namespace Conversion
MAP = {}

# global dictionaries
CPD_TO_ENZY = {}    # the key is cpd, the value is set of enzymes that act on the cpd
ENZY_TO_GENES = {}  # the key is enzymes, the value is set of keggids (genes) that make the enzyme in the SPECIES (organism)

# SPECIES (e.g. 'hsa') is a global variable
SPECIES = ''

######################################################################
def create_keggId_to_namespace_dict(outdir,to_id,withcompounds): #species is the from_id, it is always column 1; to_id is always column 2
    ## http://rest.kegg.jp/conv/<target_db>/<source_db>
    ## the source_db should be the from_id, since KEGG is the source here.
    print 'creating mapping dictionary from KEGG to %s...' % (to_id)
    from_id = SPECIES  # 'hsa'
    request_string = '/conv/' + to_id + '/' + from_id
    print '  querying database...'
    #kegg_to_namespace_mapper = CONN.request_get(request_string, headers={'Accept':'text/json'})['body']
    kegg_to_namespace_mapper = requests.get(KEGG_URL+request_string, headers={'Accept':'text/json'}).text
    print '  done querying.'
    lines = kegg_to_namespace_mapper.split('\n')  # creates a list
    outfilename = outdir+SPECIES + '-KEGGtoUniprot1.txt'
    print "Writing mapper to file " + outfilename + ". # of mappings SHOULD BE: " + str(len(lines))
    out = open(outfilename,'w')
    timestamp = strftime("%Y-%m-%d %H:%M:%S")
    out.write('#' + timestamp + '\n#KEGGID\tUniprotID\n')
    skipped = 0
    numlines = 0
    for l in lines:
        entries = l.rstrip().split('\t')
        if len(entries) == 2:
            if from_id in entries[0] and IDSYMBOL[to_id] in entries[1]: # check to make sure the species and ids match
                k = entries[0]
                u = entries[1].replace(IDSYMBOL[to_id] + ':','')
                if k not in MAP:
                    MAP[k] = set()
                MAP[k].add(u)
                out.write('%s\t%s\n' % (k,u))
                numlines+=1
        else:      
            #print 'SKIPPING',entries
            skipped+=1
    out.close()
    print 'There are %d entries in the MAP dictionary' % (len(MAP.keys()))
    print 'There are %d lines in the file %s' % (numlines,outfilename)
    print 'Skipped %d entries\n' % (skipped)

    if withcompounds:
        print 'Getting Compound Mappings (for some reason off by one)...'
        from_id = 'compound' # source is compound
        to_id = 'pubchem' # hard-coded pubchem
        request_string = '/conv/' + to_id + '/' + from_id
        print '  querying database...'
        #kegg_to_namespace_mapper = CONN.request_get(request_string, headers={'Accept':'text/json'})['body']
        kegg_to_namespace_mapper = requests.get(KEGG_URL+request_string, headers={'Accept':'text/json'}).text
        print '  done querying.'
        lines = kegg_to_namespace_mapper.split('\n')  # creates a list
        outfilename =  SPECIES + '-KEGGCompoundsToPubchem.txt'
        print "Writing mapper to file " + outfilename + ". # of mappings SHOULD BE: " + str(len(lines))
        out = open(outfilename,'w')
        timestamp = strftime("%Y-%m-%d %H:%M:%S")
        out.write('#' + timestamp + '\n#KEGGCompoundID\tPubChemID\n')
        skipped = 0
        numlines = 0
        for l in lines:
            entries = l.rstrip().split('\t')
            if len(entries) == 2:
                if 'cpd' in entries[0] and IDSYMBOL['pubchem'] in entries[1]: # check to make sure the species and ids match
                    k = entries[0]
                    u = entries[1].replace(IDSYMBOL['pubchem'] + ':','')
                    if k not in MAP:
                        MAP[k] = set()
                    MAP[k].add(u)
                    out.write('%s\t%s\n' % (k,u))
                    numlines+=1
                else:
                    skipped+=1
        out.close()
        print 'There are %d entries in the MAP dictionary' % (len(MAP.keys()))
        print 'There are %d lines in the file %s' % (numlines,outfilename)
        print 'Skipped %d entries\n' % (skipped)

    
    calc_hist(MAP, 'KEGG IDs', 'Uniprot')
    return


######################################################################
def create_keggId_to_custom_namespace_dict(filename): # 'KEGGtoUniprot.txt'
    from_id = SPECIES
    lines = readColumns(filename,1,2)  #species is the from_id, it is always column 1; to_id is always column 2
    for k,u in lines:
        if from_id in k: # add check to make sure the species match
            if k not in MAP:
                MAP[k] = set()
            MAP[k].add(u)
            #print MAP[k] , u
    print len(MAP)
    return


######################################################################
def getPathwayList():
    request_string = '/list/pathway/' + SPECIES
    #species_pathways = CONN.request_get(request_string, \
    #                                    headers={'Accept':'text/json'})
    species_pathways = requests.get(KEGG_URL+request_string, headers={'Accept':'text/json'})
    species_dict = {}
    myList = species_pathways.text.split('\n')
    for l in myList:
        if '\t' in l:
            newl = l.split('\t')
            name = newl[0].split(':')[1]
            desc = newl[1]
            species_dict[name] = desc
    return species_dict

######################################################################
def getKeggDbList(keggdb):  # keggdb can take values like compound, enzyme
    print "Getting list of %s from KEGG " %(keggdb)
    request_string = '/list/' + keggdb
    #species_keggdb = CONN.request_get(request_string, headers={'Accept':'text/json'})
    species_keggdb = requests.get(KEGG_URL+request_string, headers={'Accept':'text/json'})
    keggdb_dict = {}
    myList = species_keggdb.text.split('\n')
    for l in myList:
        if '\t' in l:
            newl = l.split('\t')
            name = newl[0].split(':')[1]
            desc = newl[1]
            keggdb_dict[name] = desc
    return keggdb_dict


######################################################################
def get_keggDbEntry_write(name,outdir):
    # get TXT file of this pathway
    #txt = CONN.request_get('/get/%s' % (name), headers={'Accept':'text/json'})['body']  # cannot get kgml for compound
    txt = requests.get(KEGG_URL+'/get/%s' % (name), headers={'Accept':'text/json'}).text  # cannot get kgml for compound
    filename = '%s/enzy-comp/%s.txt' % (outdir, name)
    try:
        writeTXT_to_file(txt, filename)
    except UnicodeEncodeError:
        print 'ERROR! %s has non-unicode characters. Ignoring the unicode characters' % (name)
        txt = txt.encode("ascii", "ignore")
        writeTXT_to_file(txt, filename)

    return

######################################################################
def parseCompound(entry,outdir): # this is the kegg id of the compound
    get_keggDbEntry_write(entry)
    collect_enzy_flag = 0
    enzy_list = list()
    filename = '%s/enzy-comp/%s.txt' % (outdir, entry)
    lines = read_file(filename)
    for l in lines:
        col = l.strip().split()
        if len(col) == 0:
             continue
        if col[0] == 'ENZYME':
            collect_enzy_flag = 1
        elif col[0] == 'DBLINKS' or col[0].isupper():  # if the DBLINKS is not present/or different order (e.g. C11882) in the file, 
            #                                          #use any (uppercase) content / label at start of line to end flag
            collect_enzy_flag = 0
        if collect_enzy_flag == 1:
            enzy_list.extend(col) 
    enzy_list = [ x for x in enzy_list if (x != 'ENZYME' and x)]               #remove 'Enzyme' and any empty strings/elements
    CPD_TO_ENZY[entry] = ','.join(set(enzy_list))
    #print entry , ":" , CPD_TO_ENZY[entry]
    #print "==============\n\n"
    return


######################################################################
def parseEnzyme(entry):
    get_keggDbEntry_write(entry)
    filename = '%s/enzy-comp/%s.txt' % (outdir, entry)
    lines = read_file(filename)
    collect_name_flag =  collect_product_flag =  collect_genes_flag = 0
    enzy_name_list = list()
    enzy_product_list = list()
    enzy_genes_list = list()
    for l in lines:
        l = l.rstrip()
        if len(l) == 0:
            continue
        if ' ' != l[0]:
            col = l.split()           
            if col[0] == 'NAME':
                collect_name_flag = 1
            elif col[0] == 'CLASS' or col[0].isupper():  # if the CLASS is not present in the file, use any (uppercase) content / label at start of line to end flag
                collect_name_flag = 0
            if col[0] == 'PRODUCT':
                collect_product_flag = 1
            elif col[0] == 'COMMENT' or col[0].isupper():  # if the COMMENT is not present in the file, use any (uppercase) content / label at start of line to end flag
                collect_product_flag = 0
            if col[0] == 'GENES':
                collect_genes_flag = 1
            elif col[0] == 'DBLINKS' or col[0].isupper():  # if the DBLINKS is not present in the file, use any (uppercase) content / label at start of line to end flag
                collect_genes_flag = 0


        if collect_name_flag == 1:
            content = l.split(";")
            content = [re.sub("NAME\s+","",w) for w in content]        #remove 'NAME' and following empty strings/elements
            content = [ x for x in content if x]              # remove empty strings
            enzy_name_list.extend(content) # extend list created by splitting with ';' not space 
        # collect products of the reaction to find enzymes which created the ligand compounds
        # this information is not always present for every enzyme (e.g. 3.2.1.49)
        if collect_product_flag == 1:
            content = l.split(";")
            content = [re.sub("PRODUCT\s+","",w) for w in content]               #remove 'PRODUCT' and following empty strings/elements
            content = [ x for x in content if x]              # remove empty strings
            enzy_product_list.extend(content) # extend list created by splitting with ';' not space 
        # collect gene names for the enzyme
        if collect_genes_flag == 1:
            l = l.replace('GENES', '')  # now split with : to get the [0] organism abbrev. and [1] the genes in that organsim separated by spaces
            content = l.split(':')
            if content[0].upper().replace(' ' , '') == SPECIES.upper():
                gene_list = [ x for x in content[1].split() if x]  #remove 'GENES' and any empty strings/elements
                enzy_genes_list.extend(gene_list) # extend list created by splitting with space, the ids are kegg ids, need to be mapped to be uniprot ids 
    #print 'Enzy name(s): ' , '\t'.join(enzy_name_list)
    #print 'Enzy product(s): ' , '\t'.join(enzy_product_list)
    ## write a if condition to check if the compound is a product of the enzyme, only then add the genes
    ENZY_TO_GENES[entry] = ','.join(set(enzy_genes_list))
    #print entry , ":" , ENZY_TO_GENES[entry]
    #print "==============\n\n"
    return


######################################################################
def parsePathway(name,withcompounds,outdir):

    ## if all three files exist, print warning
    if os.path.isfile('%s/%s-entries.tsv' % (outdir, name)) and \
            os.path.isfile('%s/%s-relations.tsv' % (outdir, name)) and \
            os.path.isfile('%s/%s-reactions.tsv' % (outdir, name)):
        print '\tWARNING All files exist. Rerunning'
        

    # get KGML file of this pathway
    if os.path.isfile('%s/%s.kgml' % (outdir,name)):
        print 'WARNING: using existing KGML in directory.'
        kgml = open('%s/%s.kgml' % (outdir,name),'r').read()
    else:
        #kgml = CONN.request_get('/get/%s/kgml' % (name), headers={'Accept':'text/json'})['body']
        kgml = requests.get(KEGG_URL+'/get/%s/kgml' % (name), headers={'Accept':'text/json'}).text
        try:
            writeKGML(name,kgml,outdir)
        except UnicodeEncodeError:
            print 'ERROR! %s has non-unicode characters.  Ignoring the unicode characters' % (name)
            kgml = kgml.encode("ascii", "ignore")
            writeKGML(name,kgml,outdir)

    # Parse it
    parsed_kgml = parseString(kgml)
    
    # Construct text version of the pathway
    entries,relations,reactions,groups = handlePathwayMap(parsed_kgml,withcompounds)

    # write pathway to file
    writeEntries(name,entries,withcompounds,outdir)
    writeRelations(name,relations,withcompounds,outdir)
    writeReactions(name,reactions,withcompounds,outdir)
    if len(groups) > 0:
        writeGroups(name,groups,entries,withcompounds,outdir)
    return

######################################################################
# From Allison
def handlePathwayMap(pathway,withcompounds):
    entrytag = pathway.getElementsByTagName("entry")

    print '\t\t--> Handling Gene Entries'
    pathway_entries,grouplist,groupSets = handlePathwayEntries(entrytag,withcompounds) ## Parse Entries from KGML
    print '\t\tGroupSets:',groupSets
    # documenting the groups
    print '\t\t--> Handling Groups'
    pathway_groups = groupSets ##handlePathwayGroups(pathway.getElementsByTagName("reaction"),entrylist)

    # need first element of pathway_entries for filtering relations and reactions.
    entrylist = [x[0] for x in pathway_entries]
    
    # parse relations from KGML
    print '\t\t--> Handling Relations'
    pathway_relations = handlePathwayRelations(pathway.getElementsByTagName("relation"),grouplist,entrylist)
    # parse reaction from KGML
    print '\t\t--> Handling Reactions'
    pathway_reactions = handlePathwayReactions(pathway.getElementsByTagName("reaction"),entrylist)

    return pathway_entries,pathway_relations,pathway_reactions,pathway_groups

######################################################################
# From Allison
def handlePathwayEntries(entries,withcompounds):
    elements = []
    typedict = {'ortholog':0,'enzyme':0,'reaction':0,'reaction':0,'gene':0,\
                    'group':0,'compound':0,'map':0,'other':0}
    missingnamespace = 0
    tot = 0
    groups = []
    groupSets = {}
    multimaps = 0
    for entry in entries:
        entryid = entry.attributes.getNamedItem("id").value
        entryname = entry.attributes.getNamedItem("name").value
        entrytype = entry.attributes.getNamedItem("type").value
        if entry.attributes.getNamedItem("link") == None:
            link = NULL_ENTRY
        else:
            link = entry.attributes.getNamedItem("link").value

        ## Currently skips any type that is NOT gene or group.
        typedict[entrytype]+=1
        if entrytype == 'gene':  # parse Uniprot
        
            # query database for Uniprot names.  
            ## Check if compounds have Uniprot ID
            entrynamelist = clean(entryname).split(',')
            entrynamelist = [e for e in entrynamelist if e != 'undefined']
            tot+=len(entrynamelist)
            uniprotlist = []
            entrylist = []
            for entry in entrynamelist:
                if entry in MAP.keys(): # in pathway hsa03320, entry 62 (hsa:7316) has no uniprot id. Hence this implementation creates a line in entry file with
                    # kegg entryid but no gene id (hsa:7316) nor uniprot id  http://www.kegg.jp/dbget-bin/www_bget?hsa:7316
                    if len(MAP[entry]) > 1:
                        multimaps+=1
                    uniprotlist.append('|'.join([e for e in MAP[entry]]))
                    entrylist.append(entry)
                else:
                    missingnamespace+=1
            if len(entrylist) == 0: # this means none of the entry ids have any namespace (uniprot) info.  ## added by Rich
                entrylist.extend(entrynamelist)  # to avoid the case of line in entry file with kegg entryid but no gene id nor namespace info.
            mylist = [entryid,','.join(entrylist),entrytype,','.join(uniprotlist),link]
            elements.append(clean(mylist))

        elif entrytype == 'group': # add to list to parse later
            groups.append(entry)
            groupSets[entryid] = [g.attributes.getNamedItem("id").value for g in entry.getElementsByTagName("component")] ## keep track of entries assigned to each group


        elif withcompounds and entrytype == 'compound': # parse PubChem

            # query database for PubChemID 
            ## Check if compounds have PubChemID
            entrynamelist = clean(entryname).split(',')
            entrynamelist = [e for e in entrynamelist if e != 'undefined']
            tot+=len(entrynamelist)
            pubchemlist = []
            entrylist = []
            for entry in entrynamelist: 
                if entry in MAP.keys():
                    # hsa04512 pathway, entry 27 (or 106) has no pubchem info for compound, hence a line in entry file with kegg entryid but no
                    # compound id (gl:G02170) nor pubchem id       http://www.kegg.jp/dbget-bin/www_bget?G02170
                    if len(MAP[entry]) > 1:
                        multimaps+=1
                    pubchemlist.append('|'.join([e for e in MAP[entry]]))
                    entrylist.append(entry)
                else:
                    missingnamespace+=1
            if len(entrylist) == 0: # this means none of the entry ids have any namespace (pubchem) info.  ## added by Rich
                entrylist.extend(entrynamelist)  # to avoid the case of line in entry file with kegg entryid but no compound id nor namespace info.
            mylist = [entryid,','.join(entrylist),entrytype,','.join(pubchemlist),link]
            elements.append(clean(mylist))

    for key in typedict.keys():
        print '\t\t%s: %d entries' % (key,typedict[key])
    print '\t\t%d of %d gene entries missing Namespace ID' % (missingnamespace,tot)
    print '\t\t%d of %d gene entries have at least one KeggID mapped to many NameSpace IDs.' % (multimaps,tot)
    print '\t\t%d Groups total with largets group having %s entries' % (len(groupSets.keys()),max([0] + [b for a,b in groupSets.items()]))
    return elements,groups, groupSets


######################################################################
# From Allison
def handlePathwayRelations(relations,groups,entrylist):
    interactions = []
    tot = 0
    kept = 0
    group_to_geneid_dict = dict()
    
    # handle group entries
   
    interactiontype='group-entry'
    entrytype = 'PPrel'
    print '\t\tProcessing %d groups...' % (len(groups))
    
    interactiontype_value = 'NA' 
    
    for group in groups:
        groupid = group.attributes.getNamedItem("id").value
        
        # get component ids that are in the group (internal IDs)
        ids = [g.attributes.getNamedItem("id").value for g in group.getElementsByTagName("component")]
        #print group.attributes.getNamedItem("id").value,ids
        
        allowed_gene_ids_for_groups = list()
        # generate all-vs-all interactions between ids
        for ida in ids:
            if ida not in entrylist:
                continue
            allowed_gene_ids_for_groups.append(ida)  #only use valid ids for groups
            for idb in ids:
                if idb not in entrylist:
                    continue
                if ida != idb:
                    mylist = [ida,idb,entrytype,interactiontype, interactiontype_value]
                    interactions.append(clean(mylist))
                    kept+=1
        group_to_geneid_dict[groupid] = ','.join(allowed_gene_ids_for_groups)
    print '\t\t%d lines from groups in relations file' % (kept)
    
    kept = 0          # reset to count the number of relations kept
    lines_group_counter = 0 # to count number of lines added from processing group gene or group group entries
    lines_nogroup_counter = 0 # to count number of lines added from processing gene gene entries
    
    for relation in relations:
        entry1 = relation.attributes.getNamedItem("entry1").value
        entry2 = relation.attributes.getNamedItem("entry2").value
        tot+=1
        if (entry1 not in entrylist and entry1 not in group_to_geneid_dict.keys()) or (entry2 not in entrylist and entry2 not in group_to_geneid_dict.keys()):
            continue
        kept+=1

        entrytype = relation.attributes.getNamedItem("type").value
        subtypelist = relation.getElementsByTagName("subtype")
        if len(subtypelist) == 0:
            interactiontype = NULL_ENTRY
        else:
            # sometimes there's a space: replace with dash
            ## add code to even get the value attribute for subtype tag
            subtypelistname = [el.attributes.getNamedItem("name").value for el in subtypelist]
            subtypelistvalue = [el.attributes.getNamedItem("value").value for el in subtypelist]
            
            for i in range(len(subtypelistname)): 
                subtypelistname[i] = re.sub(r"\s+","-",subtypelistname[i])
                subtypelistvalue[i] = re.sub(r"\s+","-",subtypelistvalue[i])

            interactiontype = ",".join(subtypelistname)
            interactiontype_value = ",".join(subtypelistvalue)
        
        if entry1 in group_to_geneid_dict.keys() and entry2 not in group_to_geneid_dict.keys():
            for e in group_to_geneid_dict[entry1].split(','):
                mylist = [e, entry2, entrytype, interactiontype, interactiontype_value]
                interactions.append(clean(mylist))
                lines_group_counter+=1 
        elif entry1 not in group_to_geneid_dict.keys() and entry2 in group_to_geneid_dict.keys():
            for e in group_to_geneid_dict[entry2].split(','):
                mylist = [entry1, e, entrytype, interactiontype, interactiontype_value]
                interactions.append(clean(mylist))
                lines_group_counter+=1 
        elif entry1 in group_to_geneid_dict.keys() and entry2 in group_to_geneid_dict.keys():
            for ida in group_to_geneid_dict[entry1].split(','):
                for idb in group_to_geneid_dict[entry2].split(','):
                    mylist = [ida,idb,entrytype,interactiontype, interactiontype_value]
                    interactions.append(clean(mylist))
                    lines_group_counter+=1 
        else:
            mylist = [entry1, entry2, entrytype, interactiontype, interactiontype_value]    
            interactions.append(clean(mylist))
            lines_nogroup_counter+=1 
    
    print '\t\t%d lines from processing group gene or group group entries in relations file' % (lines_group_counter)
    print '\t\t%d lines from processing gene gene entries in relations file' % (lines_nogroup_counter)
    print '\t\t%d of %d relations between entries were kept' % (kept,tot)
 
    return interactions

######################################################################
# From Allison
def handlePathwayReactions(reactions,entrylist):
    interactions = []
    tot = 0
    kept = 0
    for reaction in reactions:
        entryid = reaction.attributes.getNamedItem("id").value
        entryname = reaction.attributes.getNamedItem("name").value
        entrytype = reaction.attributes.getNamedItem("type").value
        substrateid = reaction.getElementsByTagName("substrate")[0].attributes.getNamedItem("id").value
        productid = reaction.getElementsByTagName("product")[0].attributes.getNamedItem("id").value
        tot+=1
        if substrateid not in entrylist or productid not in entrylist:
            continue
        kept+=1
        
        mylist = [entryid,entryname,entrytype,substrateid,productid]
        interactions.append(clean(mylist))

    print '\t\t%d of %d reactions between entries we kept' % (kept,tot)
    return interactions

######################################################################
def writeEntries(name,entries,withcompounds,outdir):
    if withcompounds:
        filename = '%s/%s-withcompounds-entries.tsv' % (outdir, name)
    else:
        filename = '%s/%s-entries.tsv' % (outdir, name)
    outfile = open(filename,'w')
    outfile.write('#KEGG Markup Language manual: http://www.kegg.jp/kegg/xml/docs/\n')
    outfile.write('#EntryID\tEntryName(s)\tType\tNamespaceMapping\tURL\n')
    for line in entries:
        outfile.write("%s\n" % "\t".join(line))
    outfile.close()
    print '\tEntry File written to %s' % (filename)
    return

######################################################################
def writeRelations(name,relations,withcompounds,outdir):
    if withcompounds:
        filename = '%s/%s-withcompounds-relations.tsv' % (outdir, name)
    else:
        filename = '%s/%s-relations.tsv' % (outdir, name)
    outfile = open(filename,'w')
    outfile.write('#KEGG Markup Language manual: http://www.kegg.jp/kegg/xml/docs/\n')
    outfile.write('#EntryID1\tEntryID2\tEntryType\tInteractionType\tInteractionType_Value\n')
    for line in relations:
        outfile.write("%s\n" % "\t".join(line))
    outfile.close()
    print '\tRelations File written to %s' % (filename)
    return

######################################################################
def writeReactions(name,reactions,withcompounds,outdir):
    if withcompounds:
        filename = '%s/%s-withcompounds-reactions.tsv' % (outdir, name)
    else:
        filename = '%s/%s-reactions.tsv' % (outdir, name)
    outfile = open(filename,'w')
    outfile.write('#KEGG Markup Language manual: http://www.kegg.jp/kegg/xml/docs/\n')
    outfile.write('#EntryID\tEntryName\tEntryType\tSubstrateID\tProductID\n')
    for line in reactions:
        outfile.write("%s\n" % "\t".join(line))
    outfile.close()
    print '\tReactions File written to %s' % (filename)
    return

#####################################################################
def writeGroups(name,groups,entries,withcompounds,outdir):
    if len(groups.keys()) < 1:
        return
    entryDict = {e[0]:e for e in entries}
    if withcompounds:
        filename = '%s/%s-withcompounds-groups.tsv' % (outdir, name)
    else:
        filename = '%s/%s-groups.tsv' % (outdir, name)
    outfile = open(filename,'w')
    outfile.write('#KEGG Markup Language manual: http://www.kegg.jp/kegg/xml/docs/\n')
    outfile.write('#GroupID\tNumEntries\tEntryIDs\tKeggIDs\tNamespaceMapping\n')
    for g,tok in groups.items():
        outfile.write("%s\n" %('\t'.join([g, str(len(tok)),";".join(tok), ';'.join([entryDict[t][1] for t in tok if t in entryDict]),';'.join([entryDict[t][3] for t in tok if t in entryDict]) ])))
    outfile.close()
    print '\tGroup File written to %s' % (filename)
    return

######################################################################
# removes extra spaces and double commas
def clean(myvar):
    if type(myvar)==type([]):
        for i in range(len(myvar)):
            myvar[i] = re.sub(r"\s+",",",myvar[i])
            myvar[i] = re.sub(",,",",",myvar[i])
    else: # string
        myvar = re.sub(r"\s+",",",myvar)
        myvar = re.sub(",,",",",myvar)
    return myvar

######################################################################
def writeKGML(name,kgml,outdir):
    filename = '%s/%s.kgml' % (outdir, name)
    outfile = open(filename,'w')
    outfile.write("%s\n" % kgml)
    outfile.close()
    print '\tKGML written to %s' % (filename)
    return



######################################################################
def main(args):
  
    # When parse_args() is called
    # optional arguments will be identified by the - prefix
    # the remaining arguments will be assumed to be positional (required)
    parser = argparse.ArgumentParser(description='Parse the XML file.')
    parser.add_argument('-s', '--species', help="Required. case insensitive symbol of the species/organism. e.g. hsa for human")
    parser.add_argument('-a', '--allpathways', action='store_true', help="get the interactions for all KEGG pathways. Either this or -p is required.")
    parser.add_argument('-p', '--pathway', action='store', help="get the interactions for only this pathway. Either this or -a is required")
    parser.add_argument('-o', '--outdir', action='store', metavar='STR', help="Required. Output directory.")
    parser.add_argument('-c', '--convertKeggIdTo', action='store', help="convert kegg id to this case insensitive id/namespace(uniprot, ncbi-gi, ncbi-geneid) on the fly. This mapper will be written to a file with timestamp.")
    parser.add_argument('-u', '--convertKeggIdUsingCustomList', action='store', help="convert kegg id to the id (namespace) in the provided two column tab delmited mapping file")
    parser.add_argument('-m', '--compounds', action='store_true', default=False, help="parse compounds (pubchem id) in addition to genes")
    parser.add_argument('-l', '--logoutput', action='store_true', help="instead of printing to output screen, write it to log file")
    #parser.add_argument('-g', '--genesymbol', action='store_true', help="write the interactions using gene symbols")
    #parser.add_argument('-e', '--entrezid', action='store_true', help="write the interactions using entrez ids")  # this is default
    #parser.add_argument('-o', '--outstring', action='store', default="today-netpath-interactions.txt", help="output string to be prefixed to filenames")

    args = parser.parse_args()
    if args.species == None :
        parser.print_help()
        sys.exit('\n--species required\n')
    if args.allpathways == False and args.pathway == None :
        parser.print_help()
        sys.exit('\nEither --allpathways (-a) or --pathway (-p) <orgNumb> required\n')

    if args.outdir == None:
        parser.print_help()
        sys.exit('\nERROR: -o (--outdir) option must be specified.\n')

    # SPECIES is a global variable. Set it.
    global SPECIES
    SPECIES = str(args.species).lower() # convert to lower case 

    old_stdout = log_file = ''
    if args.logoutput:
        old_stdout = sys.stdout
        log_filename = strftime("%Y-%m-%d") + '-' + SPECIES   + '-log_file.txt'
        log_file = open(log_filename ,"w")
        sys.stdout = log_file
        
    # Get List of species Pathways
    species_dict = getPathwayList()
    if args.pathway != None and args.pathway not in species_dict:
        args.pathway = SPECIES + args.pathway  # if number fails, try hsaNumber; this condition would be reached only if -p was a Number
        if args.pathway != None and args.pathway not in species_dict:
            sys.exit('ERROR! %s is not in species dictionary.\nValid values are: %s\n\n' % (args.pathway,str(species_dict.keys())))
   
    if args.convertKeggIdTo != None:
        create_keggId_to_namespace_dict(args.outdir,str(args.convertKeggIdTo).lower(),args.compounds)  # on the fly
    elif args.convertKeggIdUsingCustomList != None:
        create_keggId_to_custom_namespace_dict(args.convertKeggIdUsingCustomList)  # from mapperfilename
    
    if not os.path.isdir(args.outdir): # make species dir is absent
        print 'Warning: %s does not exist. Creating...' % (args.outdir)
        os.makedirs(args.outdir)
        enzy_comp_dir = args.outdir + '/enzy-comp'
        if args.compounds and not os.path.isdir(enzy_comp_dir):
            os.makedirs(enzy_comp_dir)
    
    if args.compounds:
        enzyme_dict = getKeggDbList('enzyme')
        print "Writing %d enzymes info. to files and creating dictionary[enzyme] = gene"  % (len(enzyme_dict))
        for entry in enzyme_dict:
            parseEnzyme(entry)
        calc_hist(ENZY_TO_GENES, 'KEGG IDs', 'Gene')
        print "\nwriting ENZY_TO_GENES dict to file\n"
        filename = '%s-ENZY_TO_GENES.txt' % (SPECIES)
        writeDICT(ENZY_TO_GENES, filename)
        
        compound_dict = getKeggDbList('compound')
        print "Writing %d compounds info. to files and creating dictionary[compound] = enzyme"  % (len(compound_dict))
        for entry in compound_dict:    
            parseCompound(entry)
        calc_hist(CPD_TO_ENZY, 'KEGG IDs', 'Enzyme')
        print "\nwriting CPD_TO_ENZY dict to file\n"
        filename = '%s-CPD_TO_ENZY.txt' % (SPECIES)
        writeDICT(CPD_TO_ENZY, filename)

    
    if args.allpathways:  # this evaluates to True if provided, else False
        # Parse each Pathway
        pathnum = 0
        for name in species_dict.keys():
            print '#%d of %d: PATHWAY  %s: %s' % (pathnum+1,len(species_dict.keys()),name,species_dict[name])
            parsePathway(name,args.compounds,args.outdir)
            pathnum += 1
    elif args.pathway != None:
        name = args.pathway
        print 'PATHWAY  %s: %s' % (name,species_dict[name])
        parsePathway(name,args.compounds,args.outdir)

    print 'DONE'
    
    if args.logoutput:
        sys.stdout = old_stdout
        log_file.close()
    
    print 'DONE'


if __name__=='__main__':
    main(sys.argv)
