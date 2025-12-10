from Bio import SeqIO

records = list(SeqIO.parse("/blast/db/236606_A07.fasta", "fasta"))
assert len(records) == 2, "Please make sure the file contains exactly two sequences."

def print_wrapped_seq(seq, wrap_len = 120):
    for i in range(0, len(seq), wrap_len):
        print(seq[i:i+wrap_len])

wrap_len = 120
seq1 = str(records[0].seq)
print("\n> Seq1")
print_wrapped_seq(seq1, wrap_len)

seq2_rc = str(records[1].seq.reverse_complement())
print("\n>Seq2_reverse_complement")
print_wrapped_seq(seq2_rc, wrap_len)
print()

def bases_match(b1, b2):
    return b1 == 'N' or b2 == 'N' or b1 == b2

def find_overlap(s1, s2, min_length=10, max_mismatch_ratio=0.5):
    max_len = min(len(s1), len(s2))
    for i in range(max_len, min_length - 1, -1):
        mismatch = 0
        for a, b in zip(s1[-i:], s2[:i]):
            if not bases_match(a, b):
                mismatch += 1
        if mismatch / i <= max_mismatch_ratio:
            return i
    return 0

def visualize_overlap(s1, s2, overlap_len):
    overlap_seq1 = s1[-overlap_len:]
    overlap_seq2 = s2[:overlap_len]
    
    # First line: last part of seq1
    line1 = overlap_seq1
    # Second line: match markers
    line2 = ''.join(['|' if bases_match(a, b) else ' ' for a, b in zip(overlap_seq1, overlap_seq2)])
    # Third line: first part of seq2_rc
    line3 = overlap_seq2 
    
    return f"\nOverlap visualization (length: {overlap_len}):\n{line1}\n{line2}\n{line3}\n"

overlap_len = find_overlap(seq1, seq2_rc)

if overlap_len:
    merged_seq = seq1 + seq2_rc[overlap_len:]
    print(f"✅ Overlap of {overlap_len} bases found (allowing Ns). Sequences merged.")
    print(visualize_overlap(seq1, seq2_rc, overlap_len))
    print(">merged_seq")
    print_wrapped_seq(merged_seq, wrap_len)
    print()
else:
    merged_seq = seq1 + "N" * 10 + seq2_rc
    print("⚠️ No sufficient overlap found. Sequences joined with 'N' padding.")

with open("/blast/db/236606_A07_contig2.fasta", "w") as f:
    f.write(">merged_contig\n")
    for i in range(0, len(merged_seq), 60):
        f.write(merged_seq[i:i+60] + "\n")

print("📁 Merged sequence saved as '236606_A07_contig2.fasta' in Documents folder. \n")
