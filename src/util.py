import os
import vcf
import string
import pyfaidx
from multiprocessing import Process, Value, Array

TABEX_LOC = "/home/mweerakoon/fastk/FASTK/Tabex"
HERA1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera1/solTubHeraHap1.fa"
HERA2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/hera2/solTubHeraHap2.fa"
STIEG1_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg1/solTubStiegHap1.fa"
STIEG2_REF_LOC = "/data1/phasstphase_test/potato/k_mer_plot_stuff/parent_ref/stieg2/solTubStiegHap2.fa"

def find_specific_phaseblock_kmer (k, vcf_loc, ref_loc, phase_block_required):
    ref_alt_kmer_list = []
    haplotype_alleles_list = []
    ref_location_list = []
    variant_reader = vcf.Reader(filename = vcf_loc)
    ref_fasta = pyfaidx.Fasta(ref_loc)
    print("Gathering {}-mers from {} vcf and {} ref for phase block {}".format(k, vcf_loc, ref_loc, phase_block_required))
    if k % 2 == 0:
        k_first_half_length = k // 2
        k_second_half_length = k // 2
    else:
        k_first_half_length = (k // 2) + 1
        k_second_half_length = k // 2
    for index, record in enumerate(variant_reader):
        if index % 10000 == 0:
            print("Searching progress {:.2f}%".format(100 * variant_reader.read_bytes() / variant_reader.total_bytes()))
        try:
            phase_block = record.samples[0]["PS"]
            #print(phase_block)
            if phase_block != phase_block_required:
                continue
            else:
                print("FOUND")
        except:
            continue
        # processing the required phase block
        # get all the info
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
            haplotype_alleles_list.append(record.samples[0]["GT"])
            ref_alt_kmer_list.append((ref_kmer, alt_kmer))
            ref_location_list.append((record.CHROM, record.POS))
    # find in parent stuff outside variant reader
    for (index, ref_alt) in enumerate(ref_alt_kmer_list):
        ref_kmer, alt_kmer = ref_alt
        print(haplotype_alleles_list[index])
        # Hera1 search
        print("Searching for hera1 REF")
        (hera1_ref_result, hera1_ref_flags) = search_for_kstring_in_intermediate(TABEX_LOC, HERA1_REF_LOC, ref_kmer)
    # calculate all

    return 

def look_for_stieg_ref(thread_index, k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, tabex_loc, intermediate_loc, hera_ref, stieg_ref):
    for k_index, k_string in enumerate(k_string_vec):
        ref_k_string, alt_k_string = k_string
        stieg_ref_result, stieg_ref_flags = search_for_kstring_in_intermediate(tabex_loc, stieg_ref, k_string)
        print(stieg_ref)
    return

def thread_runner_kmer_search(number_of_threads, k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, TABEX_LOC, INTERMEDIATE_LOC, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC):
    # initialize variables
    threads = [None] * number_of_threads
    # get the length
    total_len = len(k_string_vec)
    one_thread_allocation = total_len / len(threads)
    print(total_len)
    for thread_index in range(len(threads)):
        # calculate the start and end of thread allocate data for each thread
        thread_start = int(one_thread_allocation * thread_index)
        thread_end = int(thread_start + one_thread_allocation)
        thread_k_string_vec = k_string_vec[thread_start: thread_end]
        thread_haplotype_allele_vec = haplotype_allele_vec[thread_start: thread_end]
        thread_ref_loc_vec = ref_loc_vec[thread_start: thread_end]
        thread_phase_blocks = phase_blocks[thread_start: thread_end]
        # run the thread
        threads[thread_index] = Process(target=find_which_parent_contain_kstring, args=(thread_index, thread_k_string_vec, thread_haplotype_allele_vec, thread_ref_loc_vec, thread_phase_blocks, TABEX_LOC, INTERMEDIATE_LOC, HERA_UNIQUE_LOC, STIEG_UNIQUE_LOC))
        threads[thread_index].start()

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
    # concancate two heras
    hera_concat_file = "{}hera_concat.fa".format(intermediate_loc)
    input_files = "{}{}.ktab {}{}.ktab".format(intermediate_loc, parents[0].split("/")[-1], intermediate_loc, parents[1].split("/")[-1])
    command = "{} -T64 '{} = A |. B' {}".format(logex_k_loc, hera_concat_file, input_files)
    print(os.popen(command).read())
    # concancate two steig
    stieg_concat_file = "{}stieg_concat.fa".format(intermediate_loc)
    input_files = "{}{}.ktab {}{}.ktab".format(intermediate_loc, parents[2].split("/")[-1], intermediate_loc, parents[3].split("/")[-1])
    command = "{} -T64 '{} = A |. B' {}".format(logex_k_loc, stieg_concat_file, input_files)
    print(os.popen(command).read())
    # make them unique A - B hera
    hera_output_file = "{}hera_unique.fa".format(intermediate_loc)
    input_files = "{}.ktab {}.ktab".format(hera_concat_file, stieg_concat_file)
    command = "{} -T64 '{} = A - B' {}".format(logex_k_loc, hera_output_file, input_files)
    print(os.popen(command).read())
    # make them unique B - A stieg
    stieg_output_file = "{}stieg_unique.fa".format(intermediate_loc)
    input_files = "{}.ktab {}.ktab".format(hera_concat_file, stieg_concat_file)
    command = "{} -T64 '{} = B - A' {}".format(logex_k_loc, stieg_output_file, input_files)
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
            if index > 1000:
                break
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
    return ref_alt_kmer_list, haplotype_alleles_list, ref_location_list, phase_block_numbers

def find_which_parent_contain_kstring(thread_index, k_string_vec, haplotype_allele_vec, ref_loc_vec, phase_blocks, tabex_loc, intermediate_loc, hera_ref, stieg_ref):
    final_result_blocks = [(phase_blocks[0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0])] # one block is phase block, haplotype counts (hera ref, hera alt, stieg ref, stieg alt)
    concancated_ref_k_string = ""
    concancated_alt_k_string = ""
    write_path = "{}final_result_{}.txt".format(intermediate_loc, thread_index)
    for k_index, k_string in enumerate(k_string_vec):
        ref_k_string, alt_k_string = k_string
        #print(k_index, phase_blocks[k_index])
        # run tabx when 100 k strings are collected
        if k_index % 1000 == 0:
            print("Thread {} progress {}%".format(thread_index, 100 * k_index / len(k_string_vec)))
        if (k_index % 20 == 0) and (k_index != 0):
            hera_ref_result, hera_ref_flags = search_for_kstring_in_intermediate(tabex_loc, hera_ref, concancated_ref_k_string)
            stieg_ref_result, stieg_ref_flags = search_for_kstring_in_intermediate(tabex_loc, stieg_ref, concancated_ref_k_string)
            hera_alt_result, hera_alt_flags = search_for_kstring_in_intermediate(tabex_loc, hera_ref, concancated_alt_k_string)
            stieg_alt_result, stieg_alt_flags = search_for_kstring_in_intermediate(tabex_loc, stieg_ref, concancated_alt_k_string)
            concancated_ref_k_string = ""
            concancated_alt_k_string = ""
            # process the stuff put in appropriate phase block and haplotype (increment)
            local_i = 0
            for global_i in range(k_index - 20, k_index):
                #print(local_i, global_i)
                #print("Reference location: {}".format(ref_loc_vec[global_i]))
                #print("Haplotype: {}".format(haplotype_allele_vec[global_i]))
                #print("Phase block number: {}".format(phase_blocks[global_i]))
                #print("Exclusive Ref kmer found: hera {} stieg {}".format(hera_ref_result[local_i], stieg_ref_result[local_i]))
                #print("Exclusive Alt kmer found: hera {} stieg {}".format(hera_alt_result[local_i], stieg_alt_result[local_i]))
                # make a new block if the phase block is different
                if phase_blocks[global_i] != final_result_blocks[-1][0]:
                    print("New phase block processing.. {}".format(phase_blocks[global_i]))
                    with open(write_path, 'a') as fw:
                        fw.write("{}\n".format(final_result_blocks[-1]))
                    final_result_blocks.append((phase_blocks[global_i], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]))
                # continue if not satisfying any of theese conditions # flag exclusive not phased
                if hera_ref_flags[local_i] or hera_alt_flags[local_i] or stieg_alt_flags[local_i] or stieg_ref_flags[local_i]:
                    local_i += 1
                    continue
                if not (hera_ref_result[local_i]and stieg_alt_result[local_i]) or (hera_alt_result[local_i] and stieg_ref_result[local_i]):
                    local_i += 1
                    continue
                if len(haplotype_allele_vec[global_i]) != 7:
                    local_i += 1
                    continue
                else:
                    haplotype_ref_alt = [0, 0, 0, 0]
                    if haplotype_allele_vec[global_i][0] == "1":
                        haplotype_ref_alt[0] = 1
                    if haplotype_allele_vec[global_i][2] == "1":
                        haplotype_ref_alt[1] = 1
                    if haplotype_allele_vec[global_i][4] == "1":
                        haplotype_ref_alt[2] = 1
                    if haplotype_allele_vec[global_i][6] == "1":
                        haplotype_ref_alt[3] = 1
                    #print(haplotype_ref_alt, haplotype_allele_vec[global_i])]
                print("processing {} {} {} {}".format(hera_ref_result[local_i], stieg_alt_result[local_i], hera_alt_result[local_i], stieg_ref_result[local_i]))
                if hera_ref_result[local_i]:
                    # increment the hera ref haplotypes
                    if haplotype_ref_alt[0] == 0:
                        final_result_blocks[-1][1][0] += 1
                    if haplotype_ref_alt[1] == 0:
                        final_result_blocks[-1][2][0] += 1
                    if haplotype_ref_alt[2] == 0:
                        final_result_blocks[-1][3][0] += 1
                    if haplotype_ref_alt[3] == 0:
                        final_result_blocks[-1][4][0] += 1
                if hera_alt_result[local_i]:
                    # increment the hera alt haplotypes
                    if haplotype_ref_alt[0] == 1:
                        final_result_blocks[-1][1][1] += 1
                    if haplotype_ref_alt[1] == 1:
                        final_result_blocks[-1][2][1] += 1
                    if haplotype_ref_alt[2] == 1:
                        final_result_blocks[-1][3][1] += 1
                    if haplotype_ref_alt[3] == 1:
                        final_result_blocks[-1][4][1] += 1
                if stieg_ref_result[local_i]:
                    # increment the stieg ref haplotypes
                    if haplotype_ref_alt[0] == 0:
                        final_result_blocks[-1][1][2] += 1
                    if haplotype_ref_alt[1] == 0:
                        final_result_blocks[-1][2][2] += 1
                    if haplotype_ref_alt[2] == 0:
                        final_result_blocks[-1][3][2] += 1
                    if haplotype_ref_alt[3] == 0:
                        final_result_blocks[-1][4][2] += 1
                if stieg_alt_result[local_i]:
                    # increment the stieg alt haplotypes
                    if haplotype_ref_alt[0] == 1:
                        final_result_blocks[-1][1][3] += 1
                    if haplotype_ref_alt[1] == 1:
                        final_result_blocks[-1][2][3] += 1
                    if haplotype_ref_alt[2] == 1:
                        final_result_blocks[-1][3][3] += 1
                    if haplotype_ref_alt[3] == 1:
                        final_result_blocks[-1][4][3] += 1
                local_i += 1
        concancated_ref_k_string = "{} {}".format(concancated_ref_k_string, ref_k_string)
        concancated_alt_k_string = "{} {}".format(concancated_alt_k_string, alt_k_string)
    with open(write_path, 'a') as fw:
        fw.write("{}\n".format(final_result_blocks[-1]))
    return

def search_for_kstring_in_intermediate(tabex_loc, ref_loc, k_string):
    array_for_results = []
    array_for_flags = []
    command_to_run = "{} {} {}".format(tabex_loc, ref_loc, k_string)
    output = os.popen(command_to_run).read()
    split_lines = output.splitlines()
    print(output)
    for (index, split_line) in enumerate(split_lines):
        if split_line.find("Not found") == -1:
            exists = True
        else:
            exists = False
        if index >= 1:
            if (split_line[27] != "1") and (exists == True):
                array_for_flags.append(True)
            else:
                array_for_flags.append(False)
            array_for_results.append(exists)
    return array_for_results, array_for_flags