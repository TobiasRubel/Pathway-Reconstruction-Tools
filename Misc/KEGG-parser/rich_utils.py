######################################################################
# READ FILE, ARG: FILENAME
# RETURN LIST
#####################################################################
def read_file(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    print 'Read from file %s ' %(filename)
    return lines


#####################################################################
# WRITE TXT TO FILE, ARG: TXT , FILENAME 
####################################################################
def writeTXT_to_file(txt, filename):
    outfile = open(filename,'w')
    outfile.write("%s\n" % txt)
    outfile.close()
    print 'TXT written to file %s' % (filename)
    return

#####################################################################
# WRITE DICTIONARY TO FILE, ARG: DICTIONARY , FILENAME
######################################################################
def writeDICT_to_file(DICT, filename):
    outfile = open(filename,'w')
    for k, v in DICT.items():
        outfile.write("%s\t%s\n" % (k , v))
    outfile.close()
    print 'DICT written to file %s' % (filename)
    return

######################################################################
# CALC HISTOGRAM OF KEYS WITH # OF VALUES, 
# ARG: DICT , OPTIONAL ARG: strings for key and value
######################################################################
def calc_hist(DICT, k='KEYS' , v='VALUES'):
    hist = {}  # counts the number of kegg ids with particular number of values
    for key in DICT:
        l = len(DICT[key])
        if l not in hist:
            hist[l] = 0
        hist[l]+=1
    for key in sorted(hist):
        print '%d %s have %d %s entries' % (hist[key], k, key, v)
    print '\n'
    return

######################################################################
# CONVERT LIST TO DICTIONARY, 
# ARG: LIST
# OPTIONAL ARG: DELIMITER FOR SPLITTING INTO KEY & VALUE, JOINER for multiple mappings , string for filename & 0-indexed column to use as value
######################################################################
def list_to_dict(List, delim='\t', joiner=',', string="current", val_column=1):
    local_dict = dict()
    multiple_mappings = 'no'
    for l in List:
        if '#' in l:
            continue
        key = value = ''
        array = l.strip().split(delim)
        key = array[0]
        if len(array) == 2:
            value = array[1]
        else:
            # to be decided
            #value = 'NA'
            value = array[val_column]
        if key in local_dict:
            #print '%s key already exists in %s file. Adding %s to its value in dict ' %(key, string, value)
            multiple_mappings = ''
            current_values = list()
            current_values.append(local_dict[key])
            current_values.append(value)
            current_values = condense_list(current_values , joiner)
            local_dict[key] = joiner.join(current_values)
        else:
            local_dict[key] = value
    print "There is %s multiple mappings, multiple values to the key are separated by %s" % (multiple_mappings, joiner)
    return local_dict

######################################################################
# splits the elements of the list with requested delimiter, 
# returns list with unique elements
######################################################################
def condense_list(List , delim=','):
    condensed_list = list()
    for x in List:
        condensed_list.extend(x.split(delim))
    return sorted(list(set(condensed_list)))
