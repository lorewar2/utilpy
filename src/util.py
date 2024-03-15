import os
import vcf
import string
import pyfaidx

def run_fastk_make_intermediate_files(k, fast_k_loc, intermediate_loc, ref_loc):
    # get the file name from path
    file_name = ref_loc.split("/")[-1]
    intermediate_path = "{}{}".format(intermediate_loc, file_name)
    command_to_run = "{} -k{} -t1 -N\"{}\" -T64 {}".format(fast_k_loc, k, intermediate_path, ref_loc)
    # run the command
    print(command_to_run)
    print(os.popen(command_to_run).read())
    return

def unique_kmers_for_parent_from_intermediates(logex_k_loc, intermediate_loc, parents):
    letters = string.ascii_uppercase
    number_of_parents = len(parents)
    letters = letters[0: number_of_parents]
    final_file_names = []
    ABC_commands = []
    input_files = ""
    # get the reference file names
    for (i, parent) in enumerate(parents):
        file_name = parent.split("/")[-1].split(".")[0]
        final_file_names.append("{}Unique.fa".format(file_name))
        temp_string = "{}".format(letters[i])
        input_files = "{} {}".format(input_files, parent)
        for (j, parent) in enumerate(parents):
            # the j index one is normal one
            if i == j:
                continue
            else:
                temp_string = "{} - {}".format(temp_string, letters[j])
        ABC_commands.append(temp_string)
    for (index, parent) in enumerate(parents):
        # make the command
        command = "{} -T64 '{} = {}' {}".format(logex_k_loc, final_file_names[index], ABC_commands[index], input_files)
        print(command)
        # run the command
        print(os.popen(command).read())
    return

def open_vcf_and_get_k_mer(k, vcf_loc, ref_loc):
    ref_alt_kmer_list = []
    haplotype_alleles_list = []
    ref_location_list = []
    phase_block_numbers = []
    variant_reader = vcf.Reader(filename = vcf_loc)
    ref_fasta = pyfaidx.Fasta(ref_loc)
    print("Gathering {}-mers from {} vcf and {} ref".format(k, vcf_loc, ref_loc))
    if k % 2 == 0:
        k_first_half_length = k // 2
        k_second_half_length = k // 2
    else:
        k_first_half_length = (k // 2) + 1
        k_second_half_length = k // 2
    for index, record in enumerate(variant_reader):
        if index % 10000 == 0:
            print("progress {:.2f}%".format(100 * variant_reader.read_bytes() / variant_reader.total_bytes()))
        if len(record.alleles) != 3:
            continue
        ref = record.alleles[0]
        alt = record.alleles[1]
        kmer_first_half = ref_fasta[record.CHROM][record.POS - k_first_half_length - 1 : record.POS - 1]
        kmer_second_half_ref = ref_fasta[record.CHROM][record.POS - 1 + len(ref) : record.POS + k_second_half_length - 1]
        kmer_second_half_alt = ref_fasta[record.CHROM][record.POS - 1 + len(ref) : record.POS + len(ref) - len(alt) + k_second_half_length - 1]
        if ((len(kmer_first_half) + len(kmer_second_half_ref) + len(ref)) == k) and ((len(kmer_first_half) + len(kmer_second_half_alt) + len(alt)) == k):
            ref_kmer = "{}{}{}".format(kmer_first_half, ref, kmer_second_half_ref).lower()
            alt_kmer = "{}{}{}".format(kmer_first_half, alt, kmer_second_half_alt).lower()
            try:
                phase_block = record.samples[0]["PS"]
            except:
                continue
            phase_block_numbers.append(phase_block)
            haplotype_alleles_list.append(record.samples[0]["GT"])
            ref_alt_kmer_list.append((ref_kmer, alt_kmer))
            ref_location_list.append((record.CHROM, record.POS))
            # test just 100 or specified iterations
            if index > 100:
                break
            #print(ref_kmer + "appended")
    return ref_alt_kmer_list, haplotype_alleles_list, ref_location_list, phase_block_numbers

def find_which_parent_contain_kstring(k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, tabex_loc, intermediate_loc, parent_ref_vec):
    ref_counter = [0, 0, 0, 0]
    alt_counter = [0, 0, 0, 0]
    for k_index, k_string in enumerate(k_string_vec):
        ref_k_string, alt_k_string = k_string
        temp_ref_count = []
        temp_alt_count = []
        for parent_id, parent in enumerate(parent_ref_vec):
            # check 4 parents for the ref_k_string
            result = search_for_kstring_in_intermediate(tabex_loc, intermediate_loc, parent, ref_k_string)
            temp_ref_count.append(result)
            # check 4 parents for the alt_k_string
            result = search_for_kstring_in_intermediate(tabex_loc, intermediate_loc, parent, alt_k_string)
            temp_alt_count.append(result)
        print("Reference location: {}".format(ref_loc_vec[k_index]))
        print("Haplotype: {}".format(haplotype_allele_vec[k_index]))
        print("Phase block number: {}".format(phase_blocks[k_index]))
        print("Ref kmer: {}".format(ref_k_string))
        print("Alt kmer: {}".format(alt_k_string))
        print("Ref kmer found: {}".format(temp_ref_count))
        print("Alt kmer found: {}".format(temp_alt_count))
        if ((temp_ref_count[0] == True) or (temp_ref_count[1] == True)) and ((temp_ref_count[2] == False) and (temp_ref_count[3] == False)):
            print("Hera exclusive ref kmer!")
        if ((temp_alt_count[0] == True) or (temp_alt_count[1] == True)) and ((temp_alt_count[2] == False) and (temp_alt_count[3] == False)):
            print("Hera exclusive alt kmer!")
        if ((temp_ref_count[2] == True) or (temp_ref_count[3] == True)) and ((temp_ref_count[0] == False) and (temp_ref_count[1] == False)):
            print("Steig exclusive ref kmer!")
        if ((temp_alt_count[2] == True) or (temp_alt_count[3] == True)) and ((temp_alt_count[0] == False) and (temp_alt_count[1] == False)):
            print("Steig exclusive alt kmer!")
        print("\n\n")
        
    return

def search_for_kstring_in_intermediate(tabex_loc, intermediate_loc, ref_loc, k_string):
    file_name = ref_loc.split("/")[-1]
    ktab_path = "{}{}.ktab".format(intermediate_loc, file_name)
    command_to_run = "{} {} {}".format(tabex_loc, ktab_path, k_string)
    #print(command_to_run)
    output = os.popen(command_to_run).read()
    print(output.splitlines()[-1])
    if output.find("Not found") == -1:
        exists = True
    else:
        exists = False
    return exists
