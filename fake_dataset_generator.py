import yaml, sys, os, random


def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def generate_barcodes(barcode_length, n, min_barcode_distance):
    barcodes = []
    while(len(barcodes) < n):
        barcode = ''.join(random.choice('ATCG') for _ in range(barcode_length))
        if all(hamming_distance(barcode, b) > min_barcode_distance for b in barcodes):
            barcodes.append(barcode)
    return barcodes

def print_fasta_line(f, id, seq):
    f.write('>{}\n{}\n'.format(id, seq))

def print_fastq_line(f, id, seq, qual):
    f.write('@{}\n{}\n+\n{}\n'.format(id, seq, qual))  
    
def reverse_complement(s):
    return s.translate(str.maketrans('ATCG', 'TAGC'))[::-1]  
    
    
def generate_fastq_files(path, sample_profile, construct=None, max_barcode_muts=None):
    """Write a fastq file with the given parameters
    
    Arguments:
        path {str} -- where to write the fastq files
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'reference' -- sequence to use for each construct [str]
                'number_of_reads' -- number of reads to generate for each construct [list]
                'mutations' -- number of mutations to introduce in each read [list]
                'deletions' -- number of deletions to introduce in each read [list]
                'insertions' -- number of insertions to introduce in each read [list]
        construct {str} -- if specified, only generate fastq files for this construct
    """
    if construct is not None:
        f_prefix = construct
    else:
        f_prefix = path.split('/')[-1]
    fastq1_name, fastq2_name = os.path.join(path, f_prefix + '_R1.fastq'), os.path.join(path, f_prefix + '_R2.fastq')
    with open(fastq1_name, 'w') as f1, open(fastq2_name, 'w') as f2:
        # write placeholder reads
        for c, v in sample_profile.items():
            for i in range(v['number_of_reads']):
                if max_barcode_muts is not None:
                    bs, be = v['barcode_start'], v['barcode_start'] + len(v['barcodes'])
                    if len([m for m in v['mutations'][i] if m >= bs and m < be]) > max_barcode_muts:
                        continue
                sequence = sample_profile[c]['reads'][i]
                
                print_fastq_line(f1, '{}:{}'.format(c, i), sequence, 'F'*len(v['reference']))
                print_fastq_line(f2, '{}:{}'.format(c, i), reverse_complement(sequence), 'F'*len(v['reference']))
    
def generate_fasta_file(filename, sample_profile):
    """Write a fasta file with the given parameters
    
    Arguments:
        filename {str} -- where to write the fasta file
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'reads' -- sequence to use for each construct [list]
    """
    with open(filename, 'w') as f:
        for c, v in sample_profile.items():
            print_fasta_line(f, c, v['reference'])
            
def generate_dataset(path, construct_UMI_partition, L=150, UMI_start = 50, barcode_start = 100):
    """Generates fastq pair and fasta for the given construct_UMI_partition.

    Args:
        path (_type_): where to write the files
        construct_UMI_partition (dict): 
            'construct':
                'barcode': barcode sequence [str]
                'partition':
                    {UMI: # of reads} 
        L = length of reads to generate [int]
    """
    fasta = open(os.path.join(path, 'reference.fasta'), 'w')
    fastq1 = open(os.path.join(path, 'reads_R1.fastq'), 'w')
    fastq2 = open(os.path.join(path, 'reads_R2.fastq'), 'w')
    
    for construct in construct_UMI_partition:
        sequence = ''.join([random.choice('ATCG') for _ in range(L)])
        sequence = sequence[:barcode_start] + construct_UMI_partition[construct]['barcode'] + sequence[barcode_start + len(construct_UMI_partition[construct]['barcode']):]
        for UMI in construct_UMI_partition[construct]['partition']:
            temp_sequence = sequence[:UMI_start] + UMI + sequence[UMI_start + len(UMI):]
            print_fasta_line(fasta, '{}_{}'.format(construct, UMI), temp_sequence)
            for i in range(construct_UMI_partition[construct]['partition'][UMI]):
                print_fastq_line(fastq1, '{}_{}_{}'.format(construct, UMI, i), temp_sequence, 'F'*L)
                print_fastq_line(fastq2, '{}_{}_{}'.format(construct, UMI, i), reverse_complement(temp_sequence), 'F'*L)
        
    fasta.close(), fastq1.close(), fastq2.close()
    
if __name__ == '__main__':
    
    construct_UMI_partition = {
        'construct1': {
            'barcode': 'ATCGATCG',
            'partition': {
                'ATCG': 100,
                'TAGC': 100
            }
        },
        'construct2': {
            'barcode': 'ATCGATCG',
            'partition': {
                'ATCG': 100,
                'TAGC': 100
            }
        }
    }
    
    path = 'data'
    os.makedirs(path, exist_ok=True)
    
    generate_dataset(path, construct_UMI_partition)