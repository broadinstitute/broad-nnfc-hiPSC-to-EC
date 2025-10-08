from collections import Counter, defaultdict

# Initialize counters for each position: -3, -2, -1
pos_counts = [Counter() for _ in range(3)]

with open("sequences.txt") as f:
    next(f)  # skip header
    for line in f:
        seq = line.strip()
        if len(seq) >= 3:
            last3 = seq[-3:]
            for i, base in enumerate(last3):
                pos_counts[i][base] += 1

# Output results
print(f"{'Base':<5} {'Pos1':>5} {'Pos2':>5} {'Pos3':>5}")
for base in "ACGT":
    counts = [pos_counts[i][base] for i in range(3)]
    print(f"{base:<5} {counts[0]:>5} {counts[1]:>5} {counts[2]:>5}")

