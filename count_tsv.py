# Using readlines()
file1 = open('whatshap_modified.vcf', 'r')
Lines = file1.readlines()
count = 0
phased_count = 0
all_count = 0
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
        ps_available = False
        for entry in split_array[8].split(":"):
            if entry == "PS":
                ps_available = True
        if haplotype == "0|1" or haplotype == "1|0" or haplotype == "0/1" or haplotype == "1/0":
            all_count += 1
            if ps_available:
                phased_count += 1

print("all count", all_count)
print("phased count", phased_count)