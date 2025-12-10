from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass

input_file = input("Enter the path of the fasta file:")
output_file = input("Enter output fasta file name or leave it blank for auto-generation") 
if not output_file:
    output_file = f"SixFrameProteins_{input_file.rsplit('/', 1)[-1].split('.')[0]}.fasta"

records = SeqIO.parse(input_file, "fasta")

@dataclass
class Protein:
    sequence: str
    length: int

def print_wrapped_sequence(seq, wrap_len=120):
    for i in range(0, len(seq), wrap_len):
        print(seq[i:i+wrap_len])

def extract_proteins(protein_seq):
    proteins = []
    current = ""
    recording = False 
    for aa in protein_seq:
        if aa == "M" and not recording:
            current = "M"
            recording = True     
        elif aa == "*" and recording:
            proteins.append(Protein(sequence=current, length=len(current)))
            recording = False   
        elif aa and recording:
            current += aa       
    return proteins

all_proteins = []
for record in records:
    dna_id = record.id
    dna = record.seq.upper()
    print(f"\n{dna_id}")
    print_wrapped_sequence(dna)

    dna_rc = dna.reverse_complement()
    proteins_record =[]
    
    for strand, seq in [("forward", dna), ("reverse", dna_rc)]:
        for frame in range(3):
             print(f"\n==={strand.upper()} STRAND, FRAME {frame+1}===")
             framed_dna = seq[frame:]
             trim_len = len(framed_dna) - (len(framed_dna) % 3)
             coding_dna = framed_dna[:trim_len]

             translated = coding_dna.translate(to_stop=False)
             print_wrapped_sequence(translated)
             
             proteins = extract_proteins(translated)
             for i, p in enumerate(proteins, 1):
                 label = f"{dna_id}_{strand}_frame{frame+1}_Protein{i} ({p.length} aa)"
                 print(label)
                 print(p.sequence)

                 protein_record = SeqRecord(Seq(p.sequence), id=f"{label}", description="")
                 proteins_record.append(protein_record)

    all_proteins.extend(proteins_record)

with open(output_file, "w") as output_handle:
    SeqIO.write(all_proteins, output_handle, "fasta")
    print(f"All SixFrameProteins for MultipleDNASequences saved to: {output_file}")
        
                 

        

        
        
