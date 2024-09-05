# Using readlines()
file1 = open('phasstphase_modified.vcf', 'r')
Lines = file1.readlines()
count = 0
lines_to_write = []
# Strips the newline character
for line in Lines:
    count += 1
    if count % 1000 == 0:
        print(count)
    split_array = line.strip().split("\t")
    if split_array[0][0] == "c":
        #print("yes")
        #print(split_array[9])
        # extract the three points
        haplotype = split_array[9].split(":")[0]
        assert(len(haplotype) == 3)
        if haplotype == "0|1":
            haplotype = "1|0"
        haplotype_reverse = haplotype[::-1]
        new_split9 = haplotype_reverse + split_array[9][3:]
        split_array[9] = new_split9
        new_line = "\t".join(split_array)
        lines_to_write.append(new_line)
    else:
        #print("no")
        lines_to_write.append(line.strip())

with open('phasstphase_dummy.vcf', 'w') as f:
    for line in lines_to_write:
        f.write(f"{line}\n")
