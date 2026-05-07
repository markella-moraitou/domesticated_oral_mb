# This script is adapted from https://github.com/MGXlab/gtdb2cat
# Originally written by Nikos Pappas

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import os, glob

#################################
####  CREATE TAXONOMY FILES  ####
#### based on taxonomy.ipynb ####
#################################

ncbi_field_sep = '\t|\t'
ncbi_line_sep = '\t|\n'

# Define the taxonomy ranks
official_ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# Get unique taxonomies for bacteria and archaea
archaea_tax = Path("ar53_taxonomy_r220_reps.tsv")

def get_rank_level(rank_string):
    if rank_string.startswith('d__'):
        return 'domain'
    elif rank_string.startswith('p__'):
        return 'phylum'
    elif rank_string.startswith('c__'):
        return 'class'
    elif rank_string.startswith('o__'):
        return 'order'
    elif rank_string.startswith('f__'):
        return 'family'
    elif rank_string.startswith('g__'):
        return 'genus'
    elif rank_string.startswith('s__'):
        return 'species'
    else:
        print("Unknown rank for: {}".format(rank_string))
        
def create_entries_dic(taxonomy_fp):
    '''
    Parse a taxonomy file into a dictionary
    
    This is for getting each unique entry to use as
    a name and assign a unique id to.
    
    Official rank names are used keys, their members 
    are prefixed with 'p__', 'f__' etc....
    
    
    Return:
      unique_ids: dict: Keys are string, values are sets
                      {
                         'domain' :                            { 
                           'd__Archaea', 
                           'd__Bacteria'
                            },
                          'phylum' : 
                            {
                            'p__phylum1',
                            'p__phylum2'
                            },
                          ...
                         }
    '''
    
    entries_dic = dict(zip(official_ranks, [[]]*len(official_ranks)))
    
    total_entries = 0
    
    with open(taxonomy_fp, 'r') as fin:
        for line in fin:
            total_entries +=1            
            fields = [f.strip() for f in line.split('\t')]
            lineage = fields[1]
            for rank in lineage.split(';'):
                level = get_rank_level(rank)                
                if len(entries_dic[level]) == 0:
                    entries_dic[level] = [rank]
                else:
                    entries_dic[level].append(rank)
                    
    unique_ids = {k: set(v) for k,v in entries_dic.items()}
    
    print("Total entries: {}".format(total_entries))
    
    return unique_ids
    
archaeal_entries = create_entries_dic(archaea_tax)

for rank, entry in archaeal_entries.items():
    print(rank, len(entry))
    
bacteria_tax = Path("bac120_taxonomy_r220_reps.tsv")

bacterial_entries = create_entries_dic(bacteria_tax)

for rank, entry in bacterial_entries.items():
    print(rank ,len(entry))

all_entries = {}
for rank in archaeal_entries:
    all_entries[rank] = archaeal_entries[rank].union(bacterial_entries[rank])

sum([len(v) for v in all_entries.values()])

# Double check with calculations from above
for rank, entry in all_entries.items():
    print(rank ,len(entry))

# Create names.dmp file
offsets = dict(zip(official_ranks, [2, 6, 501, 1501, 5001, 20001, 80001]))

names_dic = {}
for rank in all_entries:
    for i, entry in enumerate(all_entries[rank]):
        # add the offset of the rank to the rolling counter
        uid = i + offsets[rank]
        names_dic[uid] = entry

len(names_dic)

# Create a taxonomy dir if it is not in there
taxonomy_dir = Path('taxonomy')
if not taxonomy_dir.exists():
    taxonomy_dir.mkdir()

# Specify the path to the output file
names_dmp = taxonomy_dir / Path('names.dmp')

# Dump entries in the dump
with open(names_dmp, 'w') as fout:
    # Write the root first
    root_string = ncbi_field_sep.join(['1', 'root', '', 'scientific name']) + ncbi_line_sep
    fout.write(root_string)
    for uid in names_dic:
        uid_string = ncbi_field_sep.join([str(uid), names_dic[uid], '', 'scientific name']) + ncbi_line_sep
        fout.write(uid_string)

# Create nodes.dmp file
inv_names_dic = {v:k for k,v in names_dic.items()}

len(inv_names_dic)

assert len(names_dic) == len(inv_names_dic)

nodes_dmp = 'nodes.dmp'

seen_ranks = []

root_node_string = ncbi_field_sep.join(['1', '1', 'no rank']) + ncbi_line_sep
with open(nodes_dmp, 'w') as fout:
    ## Include root information so CAT parsing works
    fout.write(root_node_string)
    
    for taxonomy_fp in [archaea_tax, bacteria_tax]:
        with open(taxonomy_fp, 'r') as fin:
            for line in fin:
                fields = [f.strip() for f in line.split('\t')]
                lineage = fields[1]
                lineage_list = lineage.split(';')
                for i, rank in enumerate(lineage_list):
                    if i == 0 and (rank not in seen_ranks):
                        child_id = inv_names_dic[rank]
                        parent_id = 1
                        # Use superkingdom instead of domain archaea and bacteria to work
                        # better with CAT add names.
                        rank_line = ncbi_field_sep.join(map(str, [child_id, parent_id, 'superkingdom'])) + ncbi_line_sep
                        fout.write(rank_line)                        
                        seen_ranks.append(rank)
                    elif i >=1 and (rank not in seen_ranks):
                        child_id = inv_names_dic[rank]
                        parent_name = lineage_list[i-1]
                        parent_id = inv_names_dic[parent_name]
                        rank_line = ncbi_field_sep.join(map(str, [child_id, parent_id, official_ranks[i]])) + ncbi_line_sep
                        fout.write(rank_line)
                        seen_ranks.append(rank)
                    else:
                        pass
        print("Parsed taxonomy: {}".format(taxonomy_fp))
#                 print("Not sure what to do here: {}".format(line))
                
len(set(seen_ranks))

##################################
####  CREATE SEQUENCE FILES   ####
#### based on sequences.ipynb ####
##################################

protdir = Path('proteins')
if not protdir.exists():
    protdir.mkdir()

# Create list of accessions to download data for
# Archaea RefSeq
os.system("cut -f1 ar53_taxonomy_r220_reps.tsv | grep ^RS | sed -e 's/^RS_//g' > archaea_refseq.txt")
os.system("ncbi-genome-download -s refseq -F protein-fasta  -A archaea_refseq.txt -p 8 -v -o proteins archaea 2> archaea_refseq.log")

# Archaea GenBank
os.system("cut -f1 ar53_taxonomy_r220_reps.tsv | grep ^GB | sed -e 's/^GB_//g' > archaea_genbank.txt")
os.system("ncbi-genome-download -s genbank -F protein-fasta -A archaea_genbank.txt -p 8 -v -o proteins archaea 2> archaea_genbank.log")

# Bacteria RefSeq
os.system("cut -f1 bac120_taxonomy_r220_reps.tsv | grep ^RS | sed -e 's/^RS_//g' > bacteria_refseq.txt")
os.system("ncbi-genome-download -s refseq -F protein-fasta -A bacteria_refseq.txt -p 8 -vv -o proteins bacteria 2> bacteria_refseq.log")

# Bacteria GenBank
os.system("cut -f1 bac120_taxonomy_r220_reps.tsv | grep ^GB | sed -e 's/^GB_//g' > bacteria_genbank.txt")
os.system("ncbi-genome-download -s genbank -F protein-fasta -A bacteria_genbank.txt -p 8 -v -o proteins bacteria 2> bacteria_genbank.log")

# Data directories
bac_prot_refseq_dir = Path("proteins/refseq/bacteria")
bac_prot_genbank_dir = Path("proteins/genbank/bacteria")
ar_prot_refseq_dir = Path("proteins/refseq/archaea")
ar_prot_genbank_dir = Path("proteins/genbank/archaea/")

# Translated taxonomies from taxonomy.ipynb

names_dmp = Path("taxonomy/names.dmp")

# Raw taxonomies
# They contain genome to taxonomy mapping

archaea_tax = Path("bac120_taxonomy_r220_reps.tsv")
bacteria_tax = Path("ar53_taxonomy_r220_reps.tsv")

# Load id to taxonomy mappings in one dir
ar_tax_dic = {}
with open(archaea_tax, 'r') as fin:
    for line in fin:
        fields = [f.strip() for f in line.split('\t')]
        ar_tax_dic[fields[0]] = fields[1]

bac_tax_dic = {}
with open(bacteria_tax, 'r') as fin:
    for line in fin:
        fields = [f.strip() for f in line.split('\t')]
        bac_tax_dic[fields[0]] = fields[1]

all_tax = {**ar_tax_dic, **bac_tax_dic}

# Create mapping of available data dirs
ar_prot_refseq = [d for d in ar_prot_refseq_dir.iterdir() if d.is_dir()]
ar_prot_genbank = [d for d in ar_prot_genbank_dir.iterdir() if d.is_dir()]
bac_prot_refseq = [d for d in bac_prot_refseq_dir.iterdir() if d.is_dir()]
bac_prot_genbank = [d for d in bac_prot_genbank_dir.iterdir() if d.is_dir()]

all_data_dirs = ar_prot_refseq + ar_prot_genbank + bac_prot_refseq + bac_prot_genbank

data_dirs = {d.name: d for d in all_data_dirs}

names_dic = {}
with open(names_dmp, 'r') as fin:
    for line in fin:
        fields = [f.strip() for f in line.split('\t')]
        names_dic[int(fields[0])] = fields[2]

inv_names_dic = {v:k for k,v in names_dic.items()}

assert len(names_dic) == len(inv_names_dic)

def get_protein_fasta_from_dir(input_dir):
    '''
    Helper function to get proteins, if available, for a genome
    '''
    try:
        faa = list(input_dir.glob('*.faa.gz'))[0]
        return faa
    except IndexError:
#         print("No proteins found for {} (path: {})".format(input_dir.name, input_dir))
        pass

def transform_fasta_headers(input_fasta_gz, uid):
    seq_counter = 1
    try:
        with gzip.open(input_fasta_gz, 'rt') as fin:
            for record in SeqIO.parse(fin, "fasta"):
                new_id = "{}_{}".format(uid, seq_counter, record.description)
                transformed_record = SeqRecord(record.seq,
                                                id=new_id,
                                               description=record.description
                                              )
                seq_counter += 1
                yield transformed_record
    except EOFError:
        # It appears some downloaded files are malformatted
        print(input_fasta_gz)
        raise

db_dir = Path("db")

if not db_dir.exists():
    db_dir.mkdir(exist_ok=True)

gtdb_nr_gz = db_dir / Path("gtdb.nr.gz")
gtdb_nr_fa = db_dir / Path("gtdb.nr.fa")
prot2acc_gz = db_dir / Path("prot.accession2taxid.gz")
prot2acc = db_dir / Path("prot.accession2taxid.txt")

# Produce accession2taxid file
missing_proteins = []
valid_proteins = []
processed_genome_ids = 0

with open(gtdb_nr_fa, 'w') as faa, open(prot2acc, 'wt') as p2acc:
    
    p2acc.write('{}\t{}\n'.format('accession.version', 'taxid'))
    
    for genome, gtdb_lineage in all_tax.items():
        genome_acc = genome[3:]
        gtdb_taxid = gtdb_lineage.split(';')[-1]
        uid = inv_names_dic[gtdb_taxid]
        
        try:
            genome_proteins = get_protein_fasta_from_dir(data_dirs[genome_acc])
            if not genome_proteins:
                missing_proteins.append(genome_acc)
            else:
                for rec in transform_fasta_headers(genome_proteins, genome):
                    SeqIO.write(rec, faa, "fasta")
                    p2acc.write("{}\t{}\n".format(rec.id, uid))
        except KeyError:
    #         print("No data dir found for {}".format(genome_acc))
            missing_proteins.append(genome_acc)
        finally:
            processed_genome_ids += 1
            if processed_genome_ids % 10 == 0:
                print("Processed {} / {} genome ids\r".format(processed_genome_ids, len(all_tax)),
                      end='')

print("Processed genome IDs:")
processed_genome_ids

print("missing proteins:" + str(len(missing_proteins)))

def create_prot2taxid(fasta_fp, inv_names_dic, all_taxonomies, prot2taxid_gz):
    with gzip.open(fasta_fp, 'rt') as fin, gzip.open(prot2taxid_gz, 'wt') as fout:
        fout.write('accession.version\ttaxid\n')
        for line in fin:
            if line.startswith('>'):
                # Protein id
                prot_id = line.split()[0].replace('>', '')
                # Genome accession
                genome_acc = '_'.join(prot_id.split('_')[:-1])
                
                gtdb_lineage = all_taxonomies[genome_acc]
                gtdb_name = gtdb_lineage.split(';')[-1]
                uid = inv_names_dic[gtdb_name]
                fout.write("{}\t{}\n".format(prot_id, uid))

create_prot2taxid(gtdb_nr_fa, inv_names_dic, all_tax, prot2acc_gz)