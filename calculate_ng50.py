def main():
    genome_length = 3_200_000_000
    file1 = open('phasstphase_modified.vcf', 'r')
    Lines = file1.readlines()
    phase_blocks = []
    count = 0
    current_phase_block = 0
    current_phase_block_length = 0
    current_phase_block_start = 0
    # Strips the newline character
    for line in Lines:
        split_array = line.strip().split("\t")
        count += 1
        if count % 5000 == 0:
            print(count)
        if split_array[0][0] == "c":
            location = int(split_array[1])
            haplotype = split_array[9].split(":")[0]
            assert(len(haplotype) == 3)
            ps_available = False
            PS_pos = 0
            for entry in split_array[8].split(":"):
                if entry == "PS":
                    ps_available = True
                    break
                PS_pos += 1
            if ps_available:
                if split_array[9].split(":")[PS_pos] == ".":
                    continue
                this_phase_block = int(split_array[9].split(":")[PS_pos])
                if this_phase_block == current_phase_block:
                    current_phase_block_length = location - current_phase_block_start
                    if current_phase_block_length > 100_000:
                        print(location, this_phase_block, "dsd")
                        print(line)
                else:
                    if current_phase_block_length > 0:
                        phase_blocks.append(current_phase_block_length)
                    current_phase_block = this_phase_block
                    current_phase_block_length = 0
                    current_phase_block_start = location
    print(phase_blocks)
    print(calculate_n50(phase_blocks))
    print(calculate_ng50(phase_blocks, genome_length))
    return



def calculate_n50(block_lengths):
    """
    Calculate N50 from a list of block lengths.
    
    Parameters:
    block_lengths (list): List of integers representing the lengths of blocks (contigs or phased blocks).

    Returns:
    int: N50 value.
    """
    # Sort block lengths in descending order
    sorted_blocks = sorted(block_lengths, reverse=True)
    
    # Calculate the total length of all blocks
    total_length = sum(sorted_blocks)
    
    # Calculate N50, the block length where 50% of the total length is covered
    cumulative_sum = 0
    half_total_length = total_length / 2
    for block_length in sorted_blocks:
        cumulative_sum += block_length
        if cumulative_sum >= half_total_length:
            return block_length
    
    return None  # If N50 is not found


def calculate_ng50(phased_blocks, genome_length):
    """
    Calculate NG50 from a list of phased block lengths.
    
    Parameters:
    phased_blocks (list): List of integers representing the lengths of the phased blocks.
    genome_length (int): Total length of the genome.

    Returns:
    int: NG50 value.
    """
    # Sort the phased blocks by length in descending order
    sorted_blocks = sorted(phased_blocks, reverse=True)
    
    # Calculate the cumulative sum of the block lengths
    cumulative_sum = 0
    
    # Calculate NG50, which is the block length where 50% of the genome length is covered
    half_genome = genome_length / 2
    for block_length in sorted_blocks:
        cumulative_sum += block_length
        if cumulative_sum >= half_genome:
            return block_length
    
    return None  # If the NG50 is not found

if __name__ == "__main__":
    main()