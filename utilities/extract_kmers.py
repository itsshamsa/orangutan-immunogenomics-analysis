sequence = input("Paste sequence: ").strip().upper().replace(" ", "")
positions = input("Mutation positions (comma-separated): ")

positions = [int(x.strip()) for x in positions.split(",")]

k = 31
half = k // 2

for position in positions:

    index = position - 1

    if index < half or index + half >= len(sequence):
        print(f"Position {position}: not enough flanking sequence.")
        continue

    original_kmer = sequence[index-half:index+half+1]

    print(f"\n### Position {position} (reference base = {sequence[index]}) ###")

    for nt in "ACGT":
        kmer = original_kmer[:half] + nt + original_kmer[half+1:]
        print(f">pos{position}_{nt}")
        print(kmer)
