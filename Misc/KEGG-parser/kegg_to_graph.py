import sys
from utilsPoirel import *
import argparse
import collections
import re
import os
from time import strftime
from rich_utils import *
import itertools
from csbdb import Interface

## Mapping file used by Chris/Anna in Signaling Pathway Prediction
##mapfile = '/data/poirel/research/signaling-pathways/data/human-ppi/human-gene-map.txt'

####### ANNA NO LONGER USES THESE #################
# KEGG Interactions of these types will be IGNORED from the network construction
IGNORED_TYPES = ['None','missing-interaction']
# KEGG Interactions of these types (edge subtypes) will be BIDIRECTED
# 'group-entry' are complexes that have been expanded to cliques (our term)
BIDIRECTED_TYPES = ['state-change','binding/association','dissociation', \
                       'compound','hidden-compound','group-entry', \
                        'indirect-effect']
# KEGG Interactions of these types (edge subtypes) will be DIRECTED
# Interactions are only directed if ALL edge subtypes are directed.
DIRECTED_TYPES = ['activation','inhibition','expression','repression', \
                      'ubiquitination','methylation','phosphorylation', \
                      'dephosphorylation']
########## END ANNA NO LONGER USES THESE ################

# SPECIES (e.g. 'hsa') is a global variable
# This is set when the user specifies the species with '-s'
SPECIES = ''

CPD_to_NAMESPACE_DICT = dict()

##################################################################################
def main(args):
  
    # When parse_args() is called
    # optional arguments will be identified by the - prefix
    # the remaining arguments will be assumed to be positional (required)
    desc = '''
Parse the XML file retrieved from the KEGG database.  -s is required, and one of -a/-p is required.

When just these three arguments are specified, the graph for signaling pathway prediction (chris/anna) is generated.  Constructing this graph requires the *-withcompounds-* files from getinteractions.py, but only takes the genes as nodes (it ignores compounds).  Edges that include non-genes are ignored, as are GErel edges (TF->gene target).

To create graph for intercellular signaling project (Rich) (i.e. map compound to uniprot ids of proteins acting on it), use all these args. (-r , -m , -e and -g). Generate the files for -e and -g from the getinteractions.py . Does not ignore GErel.

'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-s', '--species', help="Required. case insensitive symbol of the species/organism. e.g. hsa for human. Must be specified.")
    parser.add_argument('-a', '--allpathways', action='store', help="get the interactions for all KEGG pathways in the file. the file can be genereated by listpathways.py. Either this or -p is required.")
    parser.add_argument('-p', '--pathway', action='store', help="get the interactions for only this pathway. Either this or -a is required")
    parser.add_argument('-i', '--indir', action='store', metavar='STR', help="Input directory that contains '*.kgml,*-entries.tsv,*-interactions.tsv, and *-relations.tsv files.  Must be specified.")
    parser.add_argument('-o', '--outdir', action='store', metavar='STR', help="Output directory.  Must be specified.")
    parser.add_argument('-c', '--convertKeggIdTo', action='store', help="convert kegg id to this case insensitive id/namespace(uniprot, ncbi-gi, ncbi-geneid) on the fly. This mapper will be written to a file with timestamp.")
    parser.add_argument('-u', '--convertKeggIdUsingCustomList', action='store', help="convert kegg id to the id (namespace) in the provided two column tab delmited mapping file")
    parser.add_argument('-e', '--compoundToEnzyme', action='store', help="compound to enzyme acting on it mapping, two column tab delmited mapping file")
    parser.add_argument('-g', '--enzymeToKeggGene', action='store', help="enzyme to kegg entry id of gene forming it mapping, two column tab delmited mapping file")
    parser.add_argument('-m', '--compounds', action='store_true', default=False, help="parse compounds (pubchem id) in addition to genes")
    parser.add_argument('-r', '--signalingAndMetabolic', action='store_true', default=False, help="parsing specific for Rich's intercellular signaling project ")
    parser.add_argument('-l', '--logoutput', action='store_true', help="instead of printing to output screen, write it to log file")

    args = parser.parse_args()
    if args.species == None :
        parser.print_help()
        sys.exit('\n--species required\n')
    if args.allpathways == None and args.pathway == None :
        parser.print_help()
        sys.exit('\nEither --allpathways (-a) <file-containing-pathways-list> or --pathway (-p) <orgNumb> required\n')

    # SPECIES is a global variable. Set it.
    global SPECIES
    SPECIES = str(args.species).lower() # convert to lower case 

    # write log to file, not std out
    old_stdout = log_file = ''
    if args.logoutput:
        old_stdout = sys.stdout
        log_filename = strftime("%Y-%m-%d") + '-' + SPECIES  + '-log_file_edges.txt'
        log_file = open(log_filename ,"w")
        sys.stdout = log_file
    
    if args.indir == None:
        parser.print_help()
        sys.exit('\nERROR: -i (--indir) option must be specified.\n')
        
    if args.outdir == None:
        parser.print_help()
        sys.exit('\nERROR: -o (--outdir) option must be specified.\n')

    if not os.path.isdir(args.indir): # make species dir if absent
       sys.exit('\nInput directory is absent. Check dir or try running getinteractions.py\n') 
    
    if not os.path.isdir(args.outdir): # make species dir if absent
        print 'Warning: output directory is absent. Generating.'
        os.makedirs(args.outdir)

    if args.compounds:
        if args.convertKeggIdUsingCustomList == None or args.compoundToEnzyme == None or args.enzymeToKeggGene == None:
            parser.print_help()
            sys.exit('\n if using -m (--compounds), then -u (--convertKeggIdUsingCustomList) , -e (--compoundToEnzyme) and -g (enzymeToKeggGene) required. \n (i) The compounds are mapped to the enzymes acting on it, \n (ii) the enzymes are mapped to the kegg entry of genes creating it, \n (iii) these kegg ids are mapped to requested namespace \n')
        elif args.convertKeggIdUsingCustomList != None and args.compoundToEnzyme != None and args.enzymeToKeggGene != None:
            create_compound_to_namespace_dict(args.convertKeggIdUsingCustomList , args.compoundToEnzyme , args.enzymeToKeggGene)

    GL_TOT = GL_counter_PPrel = GL_counter_ECrel = GL_counter_GErel = GL_counter_PCrel = 0
    GL_numnotinentrylist = 0
    GL_numpairskeggids = 0
    GL_numpairsuniprotids = 0

    # initialize CSBDB db:
    db = Interface()
                    
    pathnum = 0
    if args.allpathways != None: 
        # read file
        lines = read_file(args.allpathways)
        # Parse each Pathway
        for l in lines:
            if '#' in l:
                continue
            pathid , pathname = l.strip().split(" --> ")
            pathid = pathid.replace('path:' , '')
            print '#%d of %d: PATHWAY  %s: %s' % (pathnum+1, len(lines), pathid , pathname)

            ##################################
            #### NETWORK CONSTRUCTION ########
            if args.signalingAndMetabolic:
                # Rich's Network (for Intercellular Communication in Rat)
                (tot, interaction_type_counter_PPrel, interaction_type_counter_ECrel, interaction_type_counter_GErel, interaction_type_counter_PCrel, numnotinentrylist, numpairskeggids, numpairsuniprotids) = create_rich_network(pathid, args.compounds,args.indir,args.outdir)
                GL_TOT = GL_TOT + tot
                GL_counter_PPrel = GL_counter_PPrel + interaction_type_counter_PPrel
                GL_counter_ECrel = GL_counter_ECrel + interaction_type_counter_ECrel
                GL_counter_GErel = GL_counter_GErel + interaction_type_counter_GErel
                GL_counter_PCrel = GL_counter_PCrel + interaction_type_counter_PCrel
                GL_numnotinentrylist = GL_numnotinentrylist  + numnotinentrylist
                GL_numpairskeggids = GL_numpairskeggids + numpairskeggids
                GL_numpairsuniprotids = GL_numpairsuniprotids + numpairsuniprotids
            else:
                # Anna's Network (for Signaling Pathway Prediction)
                create_network(pathid,db,args.indir,args.outdir)
            ##################################
            pathnum += 1
    elif args.pathway != None:
        pathid = args.pathway
        print 'PATHWAY  %s ' % (pathid)
        ##################################
        #### NETWORK CONSTRUCTION ########
        if args.signalingAndMetabolic:
            # Rich's Network (for Intercellular Communication in Rat)
            create_rich_network(pathid, args.compounds,args.indir,args.outdir) 
        else:
            # Anna's Network (for Signaling Pathway Prediction)
            create_network(pathid,db,args.indir,args.outdir)
        ##################################

    if args.signalingAndMetabolic:
        print '\nRelation Statistics over all %d PATHWAYS:' %(pathnum)
        print '%d relations were skipped b/c they were not in namespace list i.e. do not have mapping to namespace' % (GL_numnotinentrylist)
        print '%d relations parsed' % (GL_TOT)
        print '  %d pairs of internal KEGG ids (e.g. hsa:4040,hsa:4041)' % (GL_numpairskeggids)
        print '  %d PCrel , %d PPrel , %d EC rel, %d GErel' %(GL_counter_PCrel , GL_counter_PPrel , GL_counter_ECrel , GL_counter_GErel)
        print '  %d pairs of uniprot ids (e.g. Q9BQB4,O75581)' % (GL_numpairsuniprotids) 

    print 'DONE'
   
    # reset std out
    if args.logoutput:
        sys.stdout = old_stdout
        log_file.close()
    return

#############################################################
#### to do:   # while creating the kegg_entryid_to_namespace_dict, check if the entry is a gene, 
              #     if compound write code to handle that, replace compound with the uniprot id of the gene making the enzyme acting on compound
def create_compound_to_namespace_dict(keggToNamespaceFile, compToEnzyFile, enzyToKeggGeneFile ):
    # read file and convert list to dictionary  # list_to_dict(List, delim='\t', joiner=',', string="current")
    keggToNamespaceDict = list_to_dict( read_file(keggToNamespaceFile), '\t' , '|' , keggToNamespaceFile)
    compToEnzyDict = list_to_dict( read_file(compToEnzyFile), '\t' , ',' , compToEnzyFile )
    enzyToKeggGeneDict = list_to_dict( read_file(enzyToKeggGeneFile), '\t' , ',' , enzyToKeggGeneFile)
       
    print "Created dictionaries"

    # method in rich_utils.py
    calc_hist(keggToNamespaceDict, 'KEGG IDS' , 'NAMESPACE IDS')
    calc_hist(compToEnzyDict, 'COMPOUNDS' , 'ENZYMES')
    calc_hist(enzyToKeggGeneDict, 'ENZYMES' , 'KEGG IDS')

    for cpd , v in compToEnzyDict.items():
        #print cpd , v
        gene_namespaceids = list()
        for enz in v.split(','):  # when cpd has no enzy info, see case 1
            genes_keggids = 'NA'    # when enz does not hav kegg geneid  info, see case 2
            if enz in enzyToKeggGeneDict:
                genes_keggids = enzyToKeggGeneDict[enz]
            # remove the "(gene_name)" from the value e.g. 6683(SPAST),84056(KATNAL1),11104(KATNA1)
            filtered_genes_keggids = re.sub("\(.*?\)" , "" , genes_keggids )
            if '(' in filtered_genes_keggids:
                sys.exit('%s has (gene name)' %filtered_genes_keggids )
            #print "compound: %s\tkeggid_gene_making_enzy: %s" %(cpd , filtered_genes_keggids)
            for kegg_g_id in filtered_genes_keggids.split(','):
                kegg_g_id_key = SPECIES.lower() + ':' + kegg_g_id   # e.g. hsa:448835
                if kegg_g_id_key in keggToNamespaceDict:
                    gene_namespaceids.append(keggToNamespaceDict[kegg_g_id_key])
                else:                 # when no uniprot id, do nothing
                    pass
                    #print kegg_g_id_key
        #print "compound: %s\tnamespaceid_gene_making_enzy: %s" %(cpd , ",".join(gene_namespaceids))
        # case 1: cpd has no enzy info.  the dict contains cpd with no uniprot info
        # case 2: cpd has 1 enzy which does not have a gene keggid for that org. the dict contains cpd with no uniprot info
        # case 3: cpd has 1 enzy which has a gene keggid for that org. but no uniprot id,  the dict contains cpd with no uniprot info
        # case 4: cpd has 1 or more enzy, which maps to two gene keggids 1,2; where 1 maps to A and B ; 2 maps to C; then 1,2 would be printed as A|B,C
        #                                                    # gene keggids 85465,10390 maps to Q9C0D9,B3KN25|Q9Y6K0; where 10390 maps to B3KN25|Q9Y6K0
        CPD_to_NAMESPACE_DICT[cpd] = ",".join(gene_namespaceids)
    filen = SPECIES + '-CPD_TO_UNIPROT.txt'
    writeDICT_to_file(CPD_to_NAMESPACE_DICT , filen)
    return

#############################################################
def create_rich_network(pathway, withcompounds,indir,outdir):
    print 'PARSING %s' % (pathway)
    # get KEGG nodes
    entryid_to_namespace = {}
    entryid_to_entrynames = {}
    genes = {}
    entrytypes = {}
    lines = ''
    outfile = ''
    if withcompounds:
        lines = readColumns('%s/%s-withcompounds-entries.tsv' % (indir, pathway),1,2,3,4)
    else:
        lines = readColumns('%s/%s-entries.tsv' % (indir, pathway),1,2,3,4)
    for l in lines:
        entryid = l[0]
        entrytypes[entryid] = l[2]
        if entrytypes[entryid] == 'gene':
            entryid_to_entrynames[entryid] = condense_list( l[1].split(',') , '|')       
            entryid_to_namespace[entryid] = condense_list( l[3].split(',') , '|')
        elif entrytypes[entryid] == 'compound':
            # fill using compound to gene dictionary
            # the kegg entry id can map to multiple cpd ids
            all_cpds =  l[1].split(',')
            all_genes = list()
            for cmpd in all_cpds:
                cmpd = cmpd.replace("cpd:" , "")
                if cmpd in CPD_to_NAMESPACE_DICT:
                    all_genes.extend(CPD_to_NAMESPACE_DICT[cmpd].split(','))
                elif 'gl:' in cmpd:
                    print 'this compound is a glycan, so it does not have info of gene creating the enzyme acting on it'
                    continue
                elif 'dr:' in cmpd:
                    print 'this compound is a drug, not interested in info of gene interacting with it'
                    continue
                else:
                    sys.exit('Cannot find the genes for cmpd %s' %cmpd)
            entryid_to_entrynames[entryid] =  condense_list( all_cpds , '|')
            entryid_to_namespace[entryid] =  condense_list( all_genes , '|')
        entryid_to_entrynames[entryid] =  [e for e in entryid_to_entrynames[entryid] if e and e != 'NA'] # to not print empty or NA elem
        entryid_to_namespace[entryid] =  [e for e in entryid_to_namespace[entryid] if e and e != 'NA'] # to not print empty or NA elem

    # get KEGG edges
    if withcompounds:
        lines = readColumns('%s/%s-withcompounds-relations.tsv' % (indir, pathway),1,2,3,4,5)
        outfile = open('%s/%s-withcompounds-edges.txt' % (outdir, pathway),'w')
    else:
        lines = readColumns('%s/%s-relations.tsv' % (indir, pathway),1,2,3,4,5)
        outfile = open('%s/%s-edges.txt' % (outdir, pathway),'w')
    outfile.write('#Tail\tHead\tEdgeType\tEdgeSubtype\tActual_interactor_types\n')
    
    tot = interaction_type_counter_PPrel = interaction_type_counter_ECrel = interaction_type_counter_GErel = interaction_type_counter_PCrel = 0
    numnotinentrylist = 0
    numpairskeggids = 0
    numpairsuniprotids = 0
    
    for l in lines:
        edgetype = l[3]
        
        # skip if both the interactors are compounds
        if entrytypes[l[0]] == entrytypes[l[1]] and entrytypes[l[0]] == 'compound':
            #print ('%s: %s\t%s\t%s\t%s\t%s') %(l[2] , ','.join(entryid_to_entrynames[l[0]]) , ','.join(entryid_to_entrynames[l[1]]) , ','.join(entryid_to_namespace[l[0]]) , ','.join(entryid_to_namespace[l[1]]) , 'SKIPPED')
            continue
                
        ## get whether it's bidirected or not. Also check the edge types.
        ignored = True
        bidir = False
        etlist = edgetype.split(',')
        for e in etlist:
            print e
            if edgetype not in IGNORED_TYPES:
                ignored = False
            if e in BIDIRECTED_TYPES:
                bidir = True
            elif e not in DIRECTED_TYPES and e not in IGNORED_TYPES:
                sys.exit('ERROR: %s not in list' % (e))            
        if ignored:
            continue
     
        tot = tot + 1
        # mappings of kegg entry ids to namepspaces
        entry1 = entryid_to_namespace[l[0]]
        entry2 = entryid_to_namespace[l[1]]
        
        # e.g. (i) glycans may not have enzyme -> kegg gene id -> uniprot ids, so the list would be empty
        # (ii) some compounds do not have uniprot ids C00076
        # (iii) some kegg genes do not have uniprot ids rno:24207, rno:24552
        if len(entry1) == 0 or len(entry2) == 0:  
            #print ('%s: %s\t%s\t%s\t%s\t%s') %(l[2] , ','.join(entryid_to_entrynames[l[0]]) , ','.join(entryid_to_entrynames[l[1]]) , ','.join(entry1) , ','.join(entry2) , 'SKIPPED')
            numnotinentrylist = numnotinentrylist + 1
            continue

        numpairskeggids+=1
         
        if "PCrel" in l[2]:
            interaction_type_counter_PCrel = interaction_type_counter_PCrel + 1
            #print ('PCrel_Interaction: %s\t%s\t%s\t%s') %(','.join(entryid_to_entrynames[l[0]]) , ','.join(entryid_to_entrynames[l[1]]) , ','.join(entry1) , ','.join(entry2))
        elif "PPrel" in l[2]:
            interaction_type_counter_PPrel = interaction_type_counter_PPrel + 1
        elif "ECrel" in l[2]:
            interaction_type_counter_ECrel = interaction_type_counter_ECrel + 1
        elif "GErel" in l[2]:
            interaction_type_counter_GErel = interaction_type_counter_GErel + 1

        ## write all-vs-all 
        for n1 in entry1:
            for n2 in entry2:
                # kegg includes interactions where though the interaction type is PCrel, both interactors are proteins
                # below I assign the true types of the interactors
                interactor_types_true = entrytypes[l[0]] + '_' + entrytypes[l[1]]
                outfile.write('%s\t%s\t%s\t%s\t%s\n' % (n1,n2,l[2],l[3], interactor_types_true))
                numpairsuniprotids = numpairsuniprotids + 1

                if bidir:
                    interactor_types_true = entrytypes[l[1]] + '_' + entrytypes[l[0]]
                    outfile.write('%s\t%s\t%s\t%s\t%s\n' %  (n2,n1,l[2],l[3], interactor_types_true))
                    numpairsuniprotids = numpairsuniprotids + 1
        
    outfile.close()

    print '\nRelation Statistics:'
    print '%d relations were skipped b/c they were not in namespace list i.e. do not have mapping to namespace' % (numnotinentrylist)
    print '%d relations parsed' % (tot)
    print '  %d pairs of internal KEGG ids (e.g. hsa:4040,hsa:4041)' % (numpairskeggids)
    print '  %d PCrel , %d PPrel , %d EC rel, %d GErel' %(interaction_type_counter_PCrel , interaction_type_counter_PPrel , interaction_type_counter_ECrel , interaction_type_counter_GErel)
    print '  %d pairs of uniprot ids (e.g. Q9BQB4,O75581)' % (numpairsuniprotids)
    
    print '\n'
    
    return (tot, interaction_type_counter_PPrel , interaction_type_counter_ECrel , interaction_type_counter_GErel , interaction_type_counter_PCrel, numnotinentrylist, numpairskeggids, numpairsuniprotids)
        
#############################################################
## Anna's Network
def create_network(origpathway,db,indir,outdir):
    # strip pathway ID of prefix (map04310 --> 04310, etc.)
    matchObj = re.search(r'(\d+)',origpathway)
    if matchObj:
        pathway = matchObj.group(1)
    else:
        pathway = origpathway

    print '\n'+'-'*50
    print 'PARSING %s WITH ID %s\n' % (origpathway,pathway)

    # populate mapping as we query the database
    mapping = {}
    
    # make sure files exist 
    # first look for files where getinteractions.py is run WITHOUT -m (compound) option.
    # then look for files where getinteractions.py is run WITH -m (compound) option.
    fileprefix = '%s/%s%s' % (indir,SPECIES,pathway)
    entryfile = '%s-entries.tsv' % (fileprefix)
    relationsfile = '%s-relations.tsv' % (fileprefix)
    if not os.path.isfile(entryfile) or not os.path.isfile(relationsfile):
        fileprefix2 = '%s/%s-withcompounds-%s' % (indir,SPECIES,pathway)
        entryfile2 = '%s-entries.tsv' % (fileprefix2)
        relationsfile2 = '%s-relations.tsv' % (fileprefix2)
        if not os.path.isfile(entryfile2) or not os.path.isfile(relationsfile2):
            print 'Error! getinteractions.py must be run. (the -m option produces files labeled as "withcompounds").'
            if not os.path.isfile(entryfile):
                print '%s does not exist' % (entryfile)
            if not os.path.isfile(relationsfile):
                print '%s does not exist' % (relationsfile)
            if not os.path.isfile(entryfile2):
                print '%s does not exist' % (entryfile2) 
            if not os.path.isfile(relationsfile2):
                print '%s does not exist' % (relationsfile2) 
            return
        else:
            entryfile = entryfile2
            relationsfile = relationsfile2

    # get KEGG nodes
    entry = {} # entrid -> list of KEGG ids
    genes = {} # KEGG id -> list of Uniprot ids
    numnongenes = 0
    tot = 0
    multikeggids = 0
    nomapcounter = 0
    lines = readColumns(entryfile,1,2,3,4)
    for (entryid,keggids,entrytype,uniprotids) in lines:
        # Only retain entries that are of type 'gene'
        if entrytype != 'gene': 
            numnongenes+=1
            continue
        tot+=1

        # make sure the number of IDs and the number of names are the same.
        ids = keggids.split(',')
        names = uniprotids.split(',')
        if len(ids) != len(names):
            print 'WARNING! Length of idlist and namelist are different.\n%d != %d: %s and %s' % \
                (len(ids),len(names),','.join(ids),','.join(names))
            nomapcounter+=1
            continue

        # if there is more than one name for the internal KEGG id, count.
        if len(names)>1:
            multikeggids+=1

        # populated entryID -> KEGG ID list and KEGGID -> uniprot ID list.
        entry[entryid] = ids
        for i in range(len(ids)):
            genes[ids[i]] = names[i].split('|')

    print 'Entry Statistics:'
    print '%d entries skipped (non-genes)' % (numnongenes)
    print '%d entries parsed' % (tot)
    print '  %d entries contain multiple KEGG ids (usually complexes)' % (multikeggids)
    print '%d KEGG ids total. Counting how many Uniprot IDs map to each KEGG id' % (len(genes.keys()))
    print '%d ENTRIES that have at least one missing Uniprot mapping' % (nomapcounter)
    calc_hist(genes,'KEGG IDs','Uniprot ID')

    # get KEGG edges
    lines = readColumns(relationsfile,1,2,3,4)
    outfile = open('%s/%s%s-edges.txt' % (outdir, SPECIES, pathway),'w')
    outfile.write('#Tail\tHead\tTailName\tHeadName\tEdgeType\tEdgeSubtype\n')

    ignoredcounts = {}
    ignoredcounts['tot'] = 0
    keptcounts = {}
    keptcounts['tot'] = 0
    numnotinentrylist = 0
    numpairskeggids = 0
    numpairsuniprotids = 0
    numpairscommonnames = 0
    numgroupedges = 0 # num of edges from expanding a complex
    groupnodes = set() # nodes from expanding a complex -> 1
    relationinteractionse = {} # relation -> interactionlist (keggid)
    relationinteractionsu = {} # relation -> interactionlist (uniprot)
    relationedges = {} # relation -> edge list 
    for (entryid1,entryid2,entrytype,interactiontype) in lines:
  
        ## if entryids 1 and 2 aren't in entrylist, skip (they are compounds)
        if entryid1 not in entry or entryid2 not in entry:
            numnotinentrylist+=1
            continue
    
        ## if this edge is ignored, increment counters and continue.
        if isIgnoredEdge(entryid1,entryid2,entrytype,interactiontype):
            ignoredcounts['tot'] += 1
            if entrytype not in ignoredcounts:
                ignoredcounts[entrytype] = 0
            ignoredcounts[entrytype]+=1
            if interactiontype not in ignoredcounts:
                ignoredcounts[interactiontype] = 0
            ignoredcounts[interactiontype]+=1
            continue
        
        # determine edge direction
        edgeDirection = determineEdgeDirection(entryid1,entryid2,entrytype,interactiontype)

        # increment counters for edges that are kept
        keptcounts['tot'] += 1
        if entrytype not in keptcounts:
            keptcounts[entrytype] = 0
        keptcounts[entrytype]+=1
        if interactiontype not in keptcounts:
            keptcounts[interactiontype] = 0
        keptcounts[interactiontype]+=1

        ## if this is of type 'group-type', add to groupnodes
        if 'group-entry' in interactiontype:
            groupnodes.add(entryid1)
            groupnodes.add(entryid2)

        entry1 = entry[entryid1]
        entry2 = entry[entryid2]
        relationinteractionse['(%s %s)' % (entryid1,entryid2)] = [] # edges of (entry1,entry2) from interaction(entryid1,entryid2)
        relationinteractionsu['(%s %s)' % (entryid1,entryid2)] = [] # edges of (uniprot1,uniprot2) from interaction(entryid1,entryid2)
        relationedges['(%s %s)' % (entryid1,entryid2)] = [] # bi-directed edges by (uniprot1,uniprot2) from interaction(entryid1,entryid2)
        ## write all-vs-all 
        for e1,e2 in itertools.product(entry1,entry2):  
            relationinteractionse['(%s %s)' % (entryid1,entryid2)].append((e1,e2))
            numpairskeggids+=1

            uniprot1 = genes[e1]
            uniprot2 = genes[e2]

            ## NEW (Aug 2013)
            ## For each gene, find ONLY reviewed Uniprot IDs.
            ## If no Reviewed IDs exist, then use unreviewed.
            # make sure list of uniprot ids go to reviewed IDs
            uniprot1 = set()
            for g in genes[e1]:
                if db.is_reviewed(g) == True:
                    uniprot1.add(g)
            if len(uniprot1)==0: # no reviewed IDs
                uniprot1 = genes[e1] # add all unreviewed IDs

            uniprot2 = set()
            for g in genes[e2]:
                if db.is_reviewed(g) == True:
                    uniprot2.add(g)
            if len(uniprot2)==0: # no reviewed IDs
                uniprot2 = genes[e2] # add all unreviewed IDs
                
            ## Go through all PAIRS of reviewed uniprot IDs
            for u1,u2 in itertools.product(uniprot1,uniprot2):
                numpairsuniprotids+=1

                # Get GeneName from Uniprot ID
                if u1 not in mapping:    
                    # can't assume human
                    mappedset = db.map_id(u1,'UniProtKB','GeneName')
                    if len(mappedset)==0:
                        continue
                    else:
                        mapping[u1] = mappedset.pop()
                if u2 not in mapping:
                    # can't assume human
                    mappedset = db.map_id(u2,'UniProtKB','GeneName')
                    if len(mappedset)==0:
                        continue
                    else:
                        mapping[u2] = mappedset.pop()

                        
                relationinteractionsu['(%s %s)' % (entryid1,entryid2)].append((u1,u2))
                relationedges['(%s %s)' % (entryid1,entryid2)].append((u1,u2))
                numpairscommonnames+=1
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % \
                                  (u1,u2,mapping.get(u1,u1), \
                                       mapping.get(u2,u2), \
                                       entrytype,interactiontype))
                
                
                if edgeDirection=='Undirected': 
                    relationedges['(%s %s)' % (entryid1,entryid2)].append((u2,u1))
                    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % \
                                      (u2,u1,mapping.get(u2,u2),\
                                           mapping.get(u1,u1),\
                                           entrytype,interactiontype))
    outfile.close()

    print 'NOW Counting how many COMMON NAMES map to each KEGG ID'
    mappedgenes = {}
    for g in genes:
        mappedgenes[g] = [u for u in genes[g] if u in mapping]
    calc_hist(mappedgenes,'KEGG IDs','Common Name')

    print '\nRelation Statistics:'
    print '%d relations were skipped b/c they were not gene-gene interactions' % (numnotinentrylist)
    print '%d relations have interaction types we ignore' % (ignoredcounts['tot'])
    print '%d relations parsed' % (keptcounts['tot'])
    print '  %d pairs of internal KEGG ids (e.g. hsa:4040,hsa:4041)' % (numpairskeggids)
    print '  %d pairs of uniprot ids (e.g. Q9BQB4,O75581)' % (numpairsuniprotids)
    print '  %d pairs of common names (e.g. SOST,LRP6)' % (numpairscommonnames)

    edgetypes = ['ECrel','PPrel','PCrel','GErel']
    print '\nEdge Types and SubTypes of IGNORED relations:'
    print ' Types:'
    for c in sorted(ignoredcounts):
        if c in edgetypes:
            print '   %s: %d' % (c,ignoredcounts[c])
    print ' SubTypes:'
    for c in sorted(ignoredcounts):
        if c not in edgetypes and c != 'tot':
            print '   %s: %d' % (c,ignoredcounts[c])

    print '\nEdge Types and SubTypes of relations we KEPT:'
    print ' Types:'
    for c in sorted(keptcounts):
        if c in edgetypes:
            print '   %s: %d' % (c,keptcounts[c])
    print ' SubTypes:'
    for c in sorted(keptcounts):
        if c not in edgetypes and c != 'tot':
            print '   %s: %d' % (c,keptcounts[c])

    print '\nEdge Statistics:'
    print '%d edges created by complex expansion (%d nodes in these complexes)' % (numgroupedges,len(groupnodes))
    print 'Histogram of # of Interactions by Relation'
    calc_hist(relationinteractionse,'Relations','Interactions')
    print 'Histogram of # of Uniprot-Mapped Interactions by Relation (ensuring that both UniprotIDs are in Mapping file)'
    calc_hist(relationinteractionsu,'Relations','Uniprot-Mapped Interactions')
    print 'Histogram of # of Edges by Relation (bidirected edges count as 2)'
    calc_hist(relationedges,'Relations','Edges')

    return

######################################################################
## Determine whether edge is ignored from KEGG Edge Type and KEGG Edge SubType
## Edge information and quotes are from http://www.kegg.jp/kegg/xml/docs/
## Input: 
##    u,v (ints) internal KEGG Entry IDs.
##    entrytype (string) One of PPrel, ECrel, GErel, PCrel
##    interactiontypes (set of strings) a set of KEGG Subtypes (compound,activation,inhibition,etc)
## OUTPUT:
##    True if ignored, False otherwise
def isIgnoredEdge(u,v,entrytype,interactiontypes):
    
    ## DETERMINE WHETHER EDGE SHOULD BE IGNORED
    # GErel is gene expression interaction (TF -> Target Gene)
    # Ignore these edges
    if entrytype == 'GErel':
        return True

    # KEEP INDIRECT-EFFECT EDGES NOW
    # If the edge has an 'indirect-effect' subtype, ignore.
    # 'indirect-effect' = "indirect effect without molecular details"
    #elif 'indirect-effect' in interactiontypes:
    #    return True

    # If the edge has a 'state-change' subtype, ignore.
    elif 'state-change' in interactiontypes:
        if u != v: # Verify that u=v
            print 'WARNING! Edge with entry IDs %s - %s has "state-change" subtype' % (u,v)
        return True

    # If the edge has a 'missing-interaction' subtype, ignore.
    # 'missing-interaction' = missing interaction due to mutation, etc.
    elif 'missing-interaction' in interactiontypes:
        return True

    # If the edge has an 'expression' subtype, ignore.
    # FOR ALL HUMAN: 'expression' only appears alone or with 'compound'.
    # I hand-checked a few, and they all seem to be mislabeled as PPrel (should be GErel)
    elif 'expression' in interactiontypes:
        return True

    # A small handful (~10 for all human pathway entries) have NO subtype.  
    # Ignore these.
    elif interactiontypes=='None':
        return True

    # If we get this far, the edge is not ignored
    return False

######################################################################
## Determine Edge Direction from KEGG Edge Type and KEGG Edge SubType
## Edge information and quotes are from http://www.kegg.jp/kegg/xml/docs/
## Input: 
##    u,v (ints) internal KEGG Entry IDs.
##    entrytype (string) One of PPrel, ECrel, GErel, PCrel
##    interactiontypes (set of strings) a set of KEGG Subtypes (compound,activation,inhibition,etc)
## OUTPUT:
##    edgeDir (string) 'directed' or 'undirected'
def determineEdgeDirection(u,v,entrytype,interactiontypes):
    
    # If edge has 'activation' or 'inhibition' edge type, let this set the direction.
    # activation/inhibition: 'positive and negative effects which may be associated 
    # with molecular information below'
    if 'activation' in interactiontypes or 'inhibition' in interactiontypes:
        edgeDir = 'Directed'
        
    # If the edge is a molecular event (phos,dephos,glyco,ubiq,or methyl),
    # it is a directed edge.  Doesn't matter if it also has an undirected subtype.
    elif 'phosphorylation' in interactiontypes or 'dephosphorylation' in interactiontypes or \
            'glycosylation' in interactiontypes or 'ubiquitination' in interactiontypes or \
            'methylation' in interactiontypes:
        edgeDir = 'Directed'
        
    # NEW, now that we are not ignoring indirect-effect
    # If edge is 'indirect-effect', then it is directed.
    elif 'indirect-effect' in interactiontypes:
        edgeDir = 'Directed'

    # Compounds are tricky.  For now, add directed edges for compound and 'hidden compound'.  
    # compound: "shared with two successive reactions (ECrel) or intermediate of two 
    # interacting proteins (PPrel)"
    elif 'compound' in interactiontypes:
        edgeDir = 'Directed'

    # Otherwise, the edge is an undirected edge.
    # binding/association, dissociation
    # group-entry: KEGG entry that is labeled as 'group', which is interpreted as a complex.
    # group-entry is the ONLY label that is NOT from the KEGG manual (we add it when making edges
    # to represent the complex)
    elif 'binding/association' in interactiontypes or 'dissociation' in interactiontypes \
            or 'group-entry' in interactiontypes:
        edgeDir = 'Undirected'
        
    else: 
        sys.exit('ERROR! Edge direction cannot be established for edge (%s,%s) with entrytype %s and subtype %s' % \
                     (u,v,entrytype,interactiontypes))
    
    return edgeDir


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

if __name__ == '__main__':
    main(sys.argv)
