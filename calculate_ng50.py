from collections import defaultdict

def main():
    # Genome length to be used in NG50
    genome_length = 2_820_000_000
    file1 = open('/data1/phasstphase_test/hg38/vcf_test_wg_new_hic/phasstphase.vcf', 'r')
    Lines = file1.readlines()
    phase_blocks_all = []
    count = 0
    current_phase_block = 0
    current_phase_block_length = 0
    current_phase_block_start = 0
    # Read the vcf and get phase block data
    for line in Lines:
        split_array = line.strip().split("\t")
        count += 1
        if count % 5000 == 0:
            print(count)
        if split_array[0][0] == "c":
            location = int(split_array[1])
            chromosone = split_array[0]
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
                    #if current_phase_block_length > 1000_000:
                    #    print(location, this_phase_block, "dsd")
                else:
                    if current_phase_block_length > 0:
                        phase_blocks_all.append((chromosone, current_phase_block, current_phase_block_length))
                    current_phase_block = this_phase_block
                    current_phase_block_length = 0
                    current_phase_block_start = location
    # Merge the similar phase blocks
    merged_dict = defaultdict(int)
    for entry in phase_blocks_all:
        key = (entry[0], entry[1])
        merged_dict[key] += entry[2]
    phase_blocks_all_merged = [(key[0], key[1], value) for key, value in merged_dict.items()]
    phase_blocks_all_merged = sorted(phase_blocks_all_merged, reverse=True, key=lambda x: x[2])
    # Use the merged list to make the phase blocks (just lengths)
    phase_blocks = []
    for entry in phase_blocks_all_merged:
        print(entry)
        phase_blocks.append(entry[2])
    print(calculate_n50(phase_blocks))
    print(calculate_ng50(phase_blocks, genome_length))
    return

def calculate_n50(block_lengths):
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