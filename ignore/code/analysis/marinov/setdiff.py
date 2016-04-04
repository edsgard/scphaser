

from argparse import ArgumentParser

def slurp_list(list_file, list_keyfield):

    f = set()
    for line in open(list_file):
        fields = line.strip().split()
        f.add(fields[list_keyfield])
    
    return f

if '__main__' == __name__:

    #Parse arguments
    parser = ArgumentParser(description='filter tofilter_file on list_file matches')
    parser.add_argument('list_file')
    parser.add_argument('tofilter_file')
    parser.add_argument('--list_keyfield', type = int)
    parser.add_argument('--file_keyfield', type = int)
        
    args = parser.parse_args()
    list_keyfield = args.list_keyfield - 1 #0-based indexing
    file_keyfield = args.file_keyfield - 1 
    
    #Slurp the filter list as a set
    f = slurp_list(args.list_file, list_keyfield)
    
    #Filter tofilter_file entries if in f
    for line in open(args.tofilter_file):

        line = line.strip()
                
        #get fields of interest
        fields = line.split()
        key = fields[file_keyfield]
        if key not in f:
            print(line)
