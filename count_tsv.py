# Using readlines()
file1 = open('hiphase_eval.tsv', 'r')
Lines = file1.readlines()

count = 0
switch_count = 0
flip_count = 0
total_variant = 0
phased_variant = 0
# Strips the newline character
for line in Lines:
    count += 1
    if count > 1:
        split_array = line.strip().split("\t")
        switch_count += int(split_array[11].split("/")[0])
        flip_count += int(split_array[11].split("/")[1])
        total_variant += int(split_array[7])
        phased_variant += int(split_array[8])
print(switch_count)
print(flip_count)
print(total_variant)
print(phased_variant)