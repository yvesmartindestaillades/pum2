import os 
from fake_dataset_generator import generate_dataset, reverse_complement

# Read in the fasta file
def read_fasta_line(fasta):

    output_construct_UMI_partition = {}

    while line := fasta.readline():
        if line[0] == '>':
            construct_UMI = line[1:].strip()
            construct, UMI = construct_UMI.split('_')
            if construct not in output_construct_UMI_partition:
                output_construct_UMI_partition[construct] = {
                    'barcode': None,
                    'partition': {}
                }
            output_construct_UMI_partition[construct]['partition'][UMI] = 0
        else:
            sequence = line.strip()
            output_construct_UMI_partition[construct]['barcode'] = sequence[barcode_start:barcode_start+8]
    
    return output_construct_UMI_partition
    
    
def read_next_sequence(fastq):
    while line := fastq.readline():
        if line[0] == '@':
            sequence = fastq.readline().strip()
            return sequence

def is_barcode_match(barcode, sequence, barcode_start, reverse=False):
    if reverse:
        return barcode == reverse_complement(sequence)[barcode_start:barcode_start+len(barcode)] #TODO make it more efficient
    return barcode == sequence[barcode_start:barcode_start+len(barcode)]

if __name__ == '__main__':
    
    # Generate a dataset with 2 constructs, each with 2 UMIs, and 100 reads for each UMI
    input_construct_UMI_partition = {
        'construct1': {
            'barcode': 'ATCGATCG',
            'partition': {
                'ATCG': 100,
                'TAGC': 100
            }
        },
        'construct2': {
            'barcode': 'CCCACTGG',
            'partition': {
                'ATCG': 100,
                'TAGC': 100
            }
        }
    }
    
    L=150
    UMI_start = 50
    barcode_start = 100
    UMI_len = 4
    path = 'data'
    os.makedirs(path, exist_ok=True)
    
    generate_dataset(path, input_construct_UMI_partition)
    
    
    
    # Now, we want to read the data back in and compare it to the input
    fasta = open('data/reference.fasta', 'r')
    fastq1 = open('data/reads_R1.fastq', 'r')
    fastq2 = open('data/reads_R2.fastq', 'r')
    
    output_construct_UMI_partition = read_fasta_line(fasta)
    
    # read fastq1 and fastq2 4 lines at a time and update the counts in output_construct_UMI_partition
    while (f1_seq := read_next_sequence(fastq1)) and (f2_seq := read_next_sequence(fastq2)):
    
        # Try to find a barcode match
        for construct in output_construct_UMI_partition:
            if is_barcode_match(output_construct_UMI_partition[construct]['barcode'], f1_seq, barcode_start) or \
                is_barcode_match(output_construct_UMI_partition[construct]['barcode'], f2_seq, barcode_start, reverse=True):
                UMI = f1_seq[UMI_start:UMI_start+UMI_len]
                if UMI in output_construct_UMI_partition[construct]['partition']:
                    output_construct_UMI_partition[construct]['partition'][UMI] += 1
                    break
                
    fasta.close(), fastq1.close(), fastq2.close()

    # Compare the input and output
    assert input_construct_UMI_partition == output_construct_UMI_partition, 'The input and output do not match'
    print('Success!')