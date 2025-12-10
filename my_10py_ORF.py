from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass

@dataclass
class ORF:
    frame: str
    start: int
    end: int
    sequence: str
    protein: str
    
input_file = "/Users/gtan/Downloads/sequence.gb"
record = SeqIO.read(input_file, "genbank")
dna_seq = record.seq.upper()
reverse_seq = dna_seq.reverse_complement()
sequence_len = len(dna_seq)

print("Sequence ID:", record.id)
print("Description:", record.description)
print("Sequence length:", len(dna_seq))
print("First 120 bases:")
print(dna_seq[:120])

def find_orfs(protein_seq: str, frame_label: str, offset: int, is_reverse: bool) -> list[ORF]:
    orfs = []
    aa_start = None
    for i, aa in enumerate(protein_seq):
        if aa == "M" and aa_start is None:
            aa_start = i
        elif aa == "*" and aa_start is not None:
            prot = protein_seq[aa_start:i]
            dna_start = offset + aa_start*3
            dna_end = offset + i*3 +3
            if is_reverse:
                dna_start, dna_end = sequence_len - dna_end, sequence_len - dna_start
            nt_seq = (reverse_seq if is_reverse else dna_seq)[dna_start:dna_end]
            orfs.append(ORF(frame=frame_label, start=dna_start, end=dna_end, sequence=str(nt_seq), protein=prot))
            aa_start = None
    return orfs

all_orfs = []

for frame in range(3):
    framed = dna_seq[frame:]
    trim_len = len(framed) - (len(framed) % 3)
    coding = framed[:trim_len]
   
    protein = coding.translate(to_stop=False)
    orfs = find_orfs(str(protein), f"+{frame+1}", offset=frame, is_reverse=False)
    all_orfs.extend(orfs)

for frame in range(3):
    framed = reverse_seq[frame:]
    trim_len = len(framed) - (len(framed) % 3)
    coding = framed[:trim_len]

    protein = coding.translate(to_stop=False)
    orfs = find_orfs(str(protein), f"-{frame+1}", offset=frame, is_reverse=True)
    all_orfs.extend(orfs)

print(f"\nFound {len(all_orfs)} ORFs:\n")
for i, orf in enumerate(all_orfs, 1):
    print(f"ORF{i}:Frame {orf.frame}, Start{orf.start}, End{orf.end}, Length{orf.end -orf.start}")
    print(f"Protein:{orf.protein}")
    print()
      

