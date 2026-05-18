import re

def read_fasta(path):
    """Read a FASTA-like file and return a list of uppercase DNA sequences (ACGT only).
    Lines beginning with '>' are treated as headers and ignored; only A/C/G/T letters are kept.
    """
    seqs = []
    try:
        with open(path, 'r') as f:
            current = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current:
                        s = ''.join(current).upper()
                        s = re.sub('[^ACGT]', '', s)
                        if s:
                            seqs.append(s)
                        current = []
                else:
                    current.append(line)
            if current:
                s = ''.join(current).upper()
                s = re.sub('[^ACGT]', '', s)
                if s:
                    seqs.append(s)
    except FileNotFoundError:
        raise
    return seqs
