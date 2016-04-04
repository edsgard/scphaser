

from argparse import ArgumentParser

def slurp_genome(genome_file):

    g = {}
    n_chroms = 0
    seq = ''
    for line in open(genome_file):
        line = line.strip()
        if line.startswith('>'):
            if n_chroms != 0:
                g[chrom] = seq
                seq = ''
            n_chroms = n_chroms + 1
            chrom = line[1:]            
        else:
            seq = seq + line
    g[chrom] = seq
    
    return g

if '__main__' == __name__:

    #Parse arguments
    parser = ArgumentParser(description='print haplotypes given variant positions and haplo-genomes')
    parser.add_argument('pos_file')
    parser.add_argument('g1_file')
    parser.add_argument('g2_file')
    parser.add_argument('--chr_fieldpos', type = int)
    parser.add_argument('--g1_fieldpos', type = int)
    parser.add_argument('--g2_fieldpos', type = int)
        
    args = parser.parse_args()
    chr_fieldpos = args.chr_fieldpos - 1
    g1_fieldpos = args.g1_fieldpos - 1 #0-based indexing
    g2_fieldpos = args.g2_fieldpos - 1
    
    #Slurp the genomes as a dict with chr as key
    g1 = slurp_genome(args.g1_file)
    g2 = slurp_genome(args.g2_file)
    
    #Map pos_file positions line-by-line
    for line in open(args.pos_file):
        line = line.strip()
        
        #get fields of interest
        fields = line.split('\t')
        chrom = fields[chr_fieldpos]
        g1_pos = int(fields[g1_fieldpos])
        g2_pos = int(fields[g2_fieldpos])

        #map pos to allele
        g1_allele = g1[chrom][g1_pos]
        g2_allele = g2[chrom][g2_pos]

        print("\t".join([line, g1_allele, g2_allele]))
