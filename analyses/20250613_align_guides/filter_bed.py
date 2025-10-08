import sys


enh = set()
with open(sys.argv[1]) as fh:
    for line in fh:
        enh.add(line.strip())


with open(sys.argv[2]) as fh:
    for line in fh:
        if line.startswith('#') or line.startswith("chrom"):
            print(line.strip())
            continue
        fields = line.strip().split('\t')
        chrom, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
        if name in enh:
                print(line.strip())

