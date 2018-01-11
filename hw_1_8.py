import urllib.request
import gzip
import itertools

# Our nucleotide alphabet
alphabet = ["A", "C", "G", "T"]


def get_compressed_text_from_URL(URL):

    data_compressed = urllib.request.urlopen(URL).read()
    data_uncompressed = gzip.decompress(data_compressed)
    text_decoded = data_uncompressed.decode("UTF-8")
    text_lines = text_decoded.splitlines()

    return text_lines


# URLS for the genome and annotation information for Deinococcus radiodurans
# (a highly radiation resistant bacteria) Strain R1
FASTA_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/565/"\
    "GCA_000008565.1_ASM856v1/GCA_000008565.1_ASM856v1_genomic.fna.gz"
GFF_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/565/"\
    "GCA_000008565.1_ASM856v1/GCA_000008565.1_ASM856v1_genomic.gff.gz"

# A map containing nucleic acid sequences by name - we will load the genome here
genome_segments = {}

FASTA_lines = get_compressed_text_from_URL(FASTA_URL)

current_genome_segment = ""

for FASTA_line in FASTA_lines:
    
    # Check if there's an angle bracket character; this indicates a new
    # genome segment
    if len(FASTA_line) > 0 and FASTA_line[0] == ">":

        # The first space-separated part of the FASTA line is the name
        # of the genome segment
        genome_segment_name = FASTA_line[1:].split()[0]
        genome_segments[genome_segment_name] = ""
        current_genome_segment = genome_segment_name
    else:
        # If we haven't come across a genome segment yet, keep looking
        if current_genome_segment == "":
            continue

        # Otherwise append the current line to the current genome segment
        genome_segments[current_genome_segment] += FASTA_line

GFF_lines = get_compressed_text_from_URL(GFF_URL)

# Sets of coding sequences, organized by feature type
coding_sequences = {}

for line in GFF_lines:

    # Skip blank or comment lines
    if len(line) < 1 or line[0] == "#":
        continue

    annotation_elements = line.split()

    feature_type = annotation_elements[2]

    if feature_type not in coding_sequences:
        coding_sequences[feature_type] = []

    # Now we extract the sequence from the genome segment
    genome_segment_name = annotation_elements[0]
    start = int(annotation_elements[3])
    end = int(annotation_elements[4])

    sequence = genome_segments[genome_segment_name][start-1:end-1]
    coding_sequences[feature_type].append(sequence)

# Pre-populate our list of codons
codons = itertools.product(alphabet, repeat=3)
codons = ["".join(x) for x in codons]

# Reserve NNN for all invalid codons
codons.append("NNN")

# Let's tally up codon usage for each feature
for feature_type in coding_sequences:

    codon_count = {}

    for codon in codons:
        codon_count[codon] = 0

    # Loop through each sequence for this feature type
    for sequence in coding_sequences[feature_type]:

        # Count off in codons (intervals of length 3)
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]

            # Check to see if this codon has any invalid nucleotides. If so, we
            # flag it as 'NNN'
            has_invalid_nucleotide = False
            for nucleotide in codon:
                if nucleotide not in alphabet:
                    has_invalid_nucleotide = True
            if has_invalid_nucleotide:
                codon = "NNN"

            codon_count[codon] += 1

    total_codon_count = sum(codon_count.values())

    print("Codon usage for feature '%s' (%i sequence(s))" %
          (feature_type, total_codon_count))

    for codon in codons:
        codon_percent = codon_count[codon] / total_codon_count * 100
        print("%s: %i (%.2f%%)" %
              (codon, codon_count[codon], codon_percent))

