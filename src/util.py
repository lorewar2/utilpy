import os
import vcf
import string
import pyfaidx
from multiprocessing import Process, Value, Array

TABEX_LOC = "/home/mweerakoon/fastk/FASTK/Tabex"
HERA1_REF_LOC = "intermediate/solTubHeraHap1.fa.ktab"
HERA2_REF_LOC = "intermediate/solTubHeraHap2.fa.ktab"
STIEG1_REF_LOC = "intermediate/solTubStiegHap1.fa.ktab"
STIEG2_REF_LOC = "intermediate/solTubStiegHap2.fa.ktab"
GROUND1_SAVE_LOC = "intermediate/hap1_ground.fa"
GROUND2_SAVE_LOC = "intermediate/hap2_ground.fa"

# modify this to only use the variant file and use the chr and region
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
            hera_ref_result = search_for_kstring_in_intermediate(tabex_loc, hera_ref, concancated_ref_k_string)
            stieg_ref_result = search_for_kstring_in_intermediate(tabex_loc, stieg_ref, concancated_ref_k_string)
            hera_alt_result = search_for_kstring_in_intermediate(tabex_loc, hera_ref, concancated_alt_k_string)
            stieg_alt_result = search_for_kstring_in_intermediate(tabex_loc, stieg_ref, concancated_alt_k_string)
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
                # only if both hera and stieg present
                #if not ((hera_ref_result[local_i] and stieg_alt_result[local_i]) or (hera_alt_result[local_i] and stieg_ref_result[local_i])):
                #    local_i += 1
                #    continue
                # only if all are unique
                #if ((hera_ref_result[local_i] == 2) or (stieg_alt_result[local_i] == 2) or (hera_alt_result[local_i] == 2) or (stieg_ref_result[local_i] == 2)):
                #    local_i += 1
                #    continue
                # only if processed correctly by the phasstphase
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
                if (hera_ref_result[local_i] == 1) and (stieg_ref_result[local_i] == 0):
                    # increment the hera ref haplotypes
                    if haplotype_ref_alt[0] == 0:
                        final_result_blocks[-1][1][0] += 1
                    if haplotype_ref_alt[1] == 0:
                        final_result_blocks[-1][2][0] += 1
                    if haplotype_ref_alt[2] == 0:
                        final_result_blocks[-1][3][0] += 1
                    if haplotype_ref_alt[3] == 0:
                        final_result_blocks[-1][4][0] += 1
                if (hera_alt_result[local_i] == 1) and (stieg_alt_result[local_i] == 0):
                    # increment the hera alt haplotypes
                    if haplotype_ref_alt[0] == 1:
                        final_result_blocks[-1][1][1] += 1
                    if haplotype_ref_alt[1] == 1:
                        final_result_blocks[-1][2][1] += 1
                    if haplotype_ref_alt[2] == 1:
                        final_result_blocks[-1][3][1] += 1
                    if haplotype_ref_alt[3] == 1:
                        final_result_blocks[-1][4][1] += 1
                if (stieg_ref_result[local_i] == 1) and (hera_ref_result[local_i] == 0):
                    # increment the stieg ref haplotypes
                    if haplotype_ref_alt[0] == 0:
                        final_result_blocks[-1][1][2] += 1
                    if haplotype_ref_alt[1] == 0:
                        final_result_blocks[-1][2][2] += 1
                    if haplotype_ref_alt[2] == 0:
                        final_result_blocks[-1][3][2] += 1
                    if haplotype_ref_alt[3] == 0:
                        final_result_blocks[-1][4][2] += 1
                if (stieg_alt_result[local_i] == 1) and (hera_alt_result[local_i] == 0):
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
    command_to_run = "{} {} {}".format(tabex_loc, ref_loc, k_string)
    output = os.popen(command_to_run).read()
    split_lines = output.splitlines()
    #print(output)
    for (index, split_line) in enumerate(split_lines):
        try:
            count = int(split_line.split()[1])
            if count < 13:
                exists = 1 # unique
            else:
                exists = 2 # exist but not unique, too many
        except:
            exists = 0 # not found
        if index >= 1:
            array_for_results.append(exists)
    return array_for_results

# make unique hera and stieg tabx files
def unique_kmers_for_parent_from_intermediates_from_reads(logex_k_loc, intermediate_loc, hera_read, steig_read):
    # make them unique A - B hera
    hera_output_file = "{}hera_unique.fa".format(intermediate_loc)
    input_files = "{}{}.ktab {}{}.ktab".format(intermediate_loc, hera_read.split("/")[-1],intermediate_loc, steig_read.split("/")[-1])
    command = "{} -T64 '{} = A - B' {}".format(logex_k_loc, hera_output_file, input_files)
    print(os.popen(command).read())
    # make them unique B - A stieg
    stieg_output_file = "{}stieg_unique.fa".format(intermediate_loc)
    command = "{} -T64 '{} = B - A' {}".format(logex_k_loc, stieg_output_file, input_files)
    print(os.popen(command).read())
    return

# kmers from reads with defined depth
def run_fastk_make_intermediate_files_from_reads(k, fast_k_loc, intermediate_loc, ref_loc, depth):
    # get the file name from path
    file_name = ref_loc.split("/")[-1]
    intermediate_path = "{}{}".format(intermediate_loc, file_name)
    command_to_run = "{} -k{} -t{} -N\"{}\" -T64 {}".format(fast_k_loc, k, depth, intermediate_path, ref_loc)
    # run the command
    print(command_to_run)
    print(os.popen(command_to_run).read())
    return

def modify_ref_and_save_hap1_hap2(vcf_loc, ref_loc):
    variant_reader = vcf.Reader(filename = vcf_loc)
    ref_fasta = pyfaidx.Fasta(ref_loc)
    # open the ref file and save chr 1 for testing
    string_array = []
    chr1_raw = ""
    chr1_head = ""
    with open(ref_loc, 'r') as rf:
        chr1_head = rf.readline()
        print(chr1_head)
        while True:
            line = rf.readline()
            if line[0] == ">":
                break
            string_array.append(line.strip())
    print("Getting stuff done")
    chr1_list = "".join(string_array)
    chr1_raw = list(chr1_list)
    print("done chr1 raw")
    print(len(chr1_raw))
    for index, record in enumerate(variant_reader):
        if record.CHROM != "chr1":
            break
        if index % 10000 == 0:
            print("progress {:.2f}%".format(100 * variant_reader.read_bytes() / variant_reader.total_bytes()))
        if len(record.alleles) != 2 or record.samples[0]["GT"] == "1|2" or record.samples[0]["GT"] == "1/2" or record.samples[0]["GT"] == "2/1" or record.samples[0]["GT"] == "2|1" or record.samples[0]["GT"] == "0/1" or record.samples[0]["GT"] == "1/0":
            continue
        ref = record.alleles[0]
        alt = record.alleles[1]
        if len(ref) == 1 and len(alt) == 1:
            # convert to haplotype 1\
            haplo_array = record.samples[0]["GT"].split("|")
            if haplo_array[1] == "1":
                if isinstance(alt, str):
                    chr1_raw[record.POS - 1] = alt
    # write the chr1_raw to file
    chr1_final = "".join(chr1_raw)
    with open(GROUND2_SAVE_LOC, 'a') as f1:
        f1.write("{}\n".format(chr1_head))
        f1.write("{}\n".format(chr1_final))
    return

def save_alt_and_ref_kmers_in_files(k, vcf_loc, ref_loc):
    variant_reader = vcf.Reader(filename = vcf_loc)
    ref_fasta = pyfaidx.Fasta(ref_loc)
    print("Gathering {}-mers from {} vcf and {} ref".format(k, vcf_loc, ref_loc))
    if k % 2 == 0:
        k_first_half_length = k // 2
        k_second_half_length = k // 2
    else:
        k_first_half_length = (k // 2) + 1
        k_second_half_length = k // 2
    with open(GROUND1_SAVE_LOC, 'a') as f1:
        with open(GROUND2_SAVE_LOC, 'a') as f2:
            for index, record in enumerate(variant_reader):
                if index % 10000 == 0:
                    print("progress {:.2f}%".format(100 * variant_reader.read_bytes() / variant_reader.total_bytes()))
                    if index > 1000:
                        break
                if len(record.alleles) != 3 or record.samples[0]["GT"] == "1|2" or record.samples[0]["GT"] == "1/2" or record.samples[0]["GT"] == "2/1" or record.samples[0]["GT"] == "2|1":
                    continue
                ref = record.alleles[0]
                alt = record.alleles[1]
                kmer_first_half = ref_fasta[record.CHROM][record.POS - k_first_half_length - 1 : record.POS - 1]
                kmer_second_half_ref = ref_fasta[record.CHROM][record.POS - 1 + len(ref) : record.POS + k_second_half_length - 1]
                kmer_second_half_alt = ref_fasta[record.CHROM][record.POS - 1 + len(ref) : record.POS + len(ref) - len(alt) + k_second_half_length - 1]
                if ((len(kmer_first_half) + len(kmer_second_half_ref) + len(ref)) == k) and ((len(kmer_first_half) + len(kmer_second_half_alt) + len(alt)) == k):
                    ref_kmer = "{}{}{}".format(kmer_first_half, ref, kmer_second_half_ref).lower()
                    alt_kmer = "{}{}{}".format(kmer_first_half, alt, kmer_second_half_alt).lower()
                    #print(index)
                    haplo_array = record.samples[0]["GT"].split("|")
                    f1.write(">{}\n".format(index))
                    f2.write(">{}\n".format(index))
                    if haplo_array[0] == "0":
                        f1.write("{}\n".format(ref_kmer))
                    if haplo_array[0] == "1":
                        f1.write("{}\n".format(alt_kmer))
                    if haplo_array[1] == "0":
                        f2.write("{}\n".format(ref_kmer))
                    if haplo_array[1] == "1":
                        f2.write("{}\n".format(alt_kmer))
    return

def make_result_file_from_dump (processing_folder):
    files = os.listdir(processing_folder)
    files = [f for f in files if os.path.isfile(processing_folder + '/' + f)]
    # only files which are .fa
    hap1_files = []
    hap2_files = []
    hap3_files = []
    hap4_files = []
    for file in files:
        if file.split("_")[3] == "0.fa":
            hap1_files.append(("{}/{}_hera.txt".format(processing_folder, file[0:-3]), "{}/{}_stieg.txt".format(processing_folder, file[0:-3])))
        if file.split("_")[3] == "1.fa":
            hap2_files.append(("{}/{}_hera.txt".format(processing_folder, file[0:-3]), "{}/{}_stieg.txt".format(processing_folder, file[0:-3])))
        if file.split("_")[3] == "2.fa":
            hap3_files.append(("{}/{}_hera.txt".format(processing_folder, file[0:-3]), "{}/{}_stieg.txt".format(processing_folder, file[0:-3])))
        if file.split("_")[3] == "3.fa":
            hap4_files.append(("{}/{}_hera.txt".format(processing_folder, file[0:-3]), "{}/{}_stieg.txt".format(processing_folder, file[0:-3])))
    # count number in hera and steig for the phase block
    write_path = "./haplotype1_result.csv"
    for hera_stieg_file in hap1_files:
        hera_file, stieg_file = hera_stieg_file
        # get the chr and location haplotype
        with open(hera_file, 'r') as fp:
            hera_lines = len(fp.readlines())
        with open(stieg_file, 'r') as fp:
            stieg_lines = len(fp.readlines())
        write_line = "{},{},{},{}\n".format(hera_file.split("/")[-1].split("_")[-5], hera_file.split("_")[-4], hera_lines, stieg_lines)
        with open(write_path, 'a') as fw:
            fw.write(write_line)
    write_path = "./haplotype2_result.csv"
    for hera_stieg_file in hap2_files:
        hera_file, stieg_file = hera_stieg_file
        # get the chr and location haplotype
        with open(hera_file, 'r') as fp:
            hera_lines = len(fp.readlines())
        with open(stieg_file, 'r') as fp:
            stieg_lines = len(fp.readlines())
        write_line = "{},{},{},{}\n".format(hera_file.split("/")[-1].split("_")[-5], hera_file.split("_")[-4], hera_lines, stieg_lines)
        with open(write_path, 'a') as fw:
            fw.write(write_line)
    write_path = "./haplotype3_result.csv"
    for hera_stieg_file in hap3_files:
        hera_file, stieg_file = hera_stieg_file
        # get the chr and location haplotype
        with open(hera_file, 'r') as fp:
            hera_lines = len(fp.readlines())
        with open(stieg_file, 'r') as fp:
            stieg_lines = len(fp.readlines())
        write_line = "{},{},{},{}\n".format(hera_file.split("/")[-1].split("_")[-5], hera_file.split("_")[-4], hera_lines, stieg_lines)
        with open(write_path, 'a') as fw:
            fw.write(write_line)
    write_path = "./haplotype4_result.csv"
    for hera_stieg_file in hap4_files:
        hera_file, stieg_file = hera_stieg_file
        # get the chr and location haplotype
        with open(hera_file, 'r') as fp:
            hera_lines = len(fp.readlines())
        with open(stieg_file, 'r') as fp:
            stieg_lines = len(fp.readlines())
        write_line = "{},{},{},{}\n".format(hera_file.split("/")[-1].split("_")[-5], hera_file.split("_")[-4], hera_lines, stieg_lines)
        with open(write_path, 'a') as fw:
            fw.write(write_line)
    return

def make_kmc_files_and_dump (processing_folder):
    # find all the files in the folder do today
    files = os.listdir(processing_folder)
    files = [f for f in files if os.path.isfile(processing_folder + '/' + f)]
    for file in files:
        if file != "CM034976.1_156_hap_0.fa":
            continue
        print("Processing {}".format(file))
        # make kmc files for all do today
        #~/kmc/bin/kmc -k21 -fa -ci1 -r 156_hap_1.fa 156_hap_1 .
        phase_fa = "{}/{}".format(processing_folder, file)
        phase_kmc = "{}/{}".format(processing_folder, file[0:-3])
        print(phase_fa)
        print(phase_kmc)
        command_to_run = "~/kmc/bin/kmc -k21 -fa -ci1 -r {} {} .".format(phase_fa, phase_kmc)
        print(command_to_run)
        os.system(command_to_run)
        # make intercept file for both do today
        # hera
        hera_intercept_kmc = "{}_hera".format(phase_kmc)
        stieg_intercept_kmc = "{}_stieg".format(phase_kmc)
        command_to_run = "~/kmc/bin/kmc_tools simple {} {}/../hera_unique intersect {}".format(phase_kmc, processing_folder, hera_intercept_kmc)
        print(command_to_run)
        os.system(command_to_run)
        command_to_run = "~/kmc/bin/kmc_tools simple {} {}/../stieg_unique intersect {}".format(phase_kmc, processing_folder, stieg_intercept_kmc)
        print(command_to_run)
        os.system(command_to_run)
        #~/kmc/bin/kmc_tools simple ./156_hap_1 ./hera_unique intersect test
        # make dump text file for both
        command_to_run = "~/kmc/bin/kmc_tools transform {} dump {}.txt".format(hera_intercept_kmc, hera_intercept_kmc)
        print(command_to_run)
        os.system(command_to_run)
        command_to_run = "~/kmc/bin/kmc_tools transform {} dump {}.txt".format(stieg_intercept_kmc, stieg_intercept_kmc)
        print(command_to_run)
        os.system(command_to_run)
        # find the number of lines in both
        command_to_run = "cat {}.txt | wc -l".format(hera_intercept_kmc)
        print(command_to_run)
        os.system(command_to_run)
        command_to_run = "cat {}.txt | wc -l".format(stieg_intercept_kmc)
        print(command_to_run)
        os.system(command_to_run)
    return

def make_fasta_file_for_each_phase_block_haplotype(k, vcf_loc, ref_loc, save_path):
    # for start just do 1 phase block 1 haplotype
    variant_reader = vcf.Reader(filename = vcf_loc)
    ref_fasta = pyfaidx.Fasta(ref_loc)
    print("Gathering {}-mers from {} vcf and {} ref and saving in fasta".format(k, vcf_loc, ref_loc))
    if k % 2 == 0:
        k_first_half_length = k // 2
        k_second_half_length = k // 2
    else:
        k_first_half_length = (k // 2) + 1
        k_second_half_length = k // 2
    for index, record in enumerate(variant_reader):
        if index % 10000 == 0:
            print("Searching progress {:.2f}%".format(100 * variant_reader.read_bytes() / variant_reader.total_bytes()))
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
            ref_chromosone = record.CHROM
            ref_location = record.POS
            try:
                phase_block = record.samples[0]["PS"]
            except:
                continue
            temp_haplo = record.samples[0]["GT"]
            if len(temp_haplo) != 7:
                continue
            haplotype_ref_alt = [0, 0, 0, 0]
            if temp_haplo[0] == "1":
                haplotype_ref_alt[0] = 1
            if temp_haplo[2] == "1":
                haplotype_ref_alt[1] = 1
            if temp_haplo[4] == "1":
                haplotype_ref_alt[2] = 1
            if temp_haplo[6] == "1":
                haplotype_ref_alt[3] = 1
            # write to fa file
            for (index, value) in enumerate(haplotype_ref_alt):
                write_path = "{}/{}_{}_hap_{}.fa".format(save_path, ref_chromosone, phase_block, index)
                with open(write_path, 'a') as fw:
                    if value == 1:
                        fw.write(">{}_{}\n{}\n".format(phase_block, ref_location, alt_kmer))
                    else:
                        fw.write(">{}_{}\n{}\n".format(phase_block, ref_location, ref_kmer))
    return

def find_specific_phaseblock_kmer (k, vcf_loc, ref_loc, phase_block_required, unique_hera, unique_stieg):
    ref_alt_kmer_list = []
    haplotype_alleles_list = []
    ref_location_list = []
    found_bool = False
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
                if found_bool == True: # encountered different block after found
                    break
                continue
            else:
                print("FOUND")
                found_bool = True
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
    # find in parent stuff outside variant reader haplotype 1
    haplotype1_stuff = []
    for (index, ref_alt) in enumerate(ref_alt_kmer_list):
        ref_kmer, alt_kmer = ref_alt
        print(haplotype_alleles_list[index])
        if len(haplotype_alleles_list[index]) != 7:
            continue
        haplotype_ref_alt = [0, 0, 0, 0]
        if haplotype_alleles_list[index][0] == "1":
            haplotype_ref_alt[0] = 1
        if haplotype_alleles_list[index][2] == "1":
            haplotype_ref_alt[1] = 1
        if haplotype_alleles_list[index][4] == "1":
            haplotype_ref_alt[2] = 1
        if haplotype_alleles_list[index][6] == "1":
            haplotype_ref_alt[3] = 1
        print("Searching for REF")
        (hera_ref_result) = search_for_kstring_in_intermediate(TABEX_LOC, unique_hera, ref_kmer)
        (stieg_ref_result) = search_for_kstring_in_intermediate(TABEX_LOC, unique_stieg, ref_kmer)
        print("Searching for ALT")
        (hera_alt_result) = search_for_kstring_in_intermediate(TABEX_LOC, unique_hera, alt_kmer)
        (stieg_alt_result) = search_for_kstring_in_intermediate(TABEX_LOC, unique_stieg, alt_kmer)
        # if not hera stieg together skip
        if not (((hera_ref_result[0]) and (stieg_alt_result[0])) or ((hera_alt_result[0]) and (stieg_ref_result[0]))):
            continue
        # if something is not unique skip
        if (hera_ref_result[0] == 2) or (hera_alt_result[0] == 2) or (stieg_ref_result[0] == 2) or (stieg_alt_result[0] == 2):
            continue
        if haplotype_ref_alt[0] == 0:
            if ((hera_ref_result[0]) and not (stieg_ref_result[0])):
                haplotype1_stuff.append((ref_location_list[index], 1, 0, 0, 0))
                print((ref_location_list[index], 1, 0, 0, 0))
            elif ((stieg_ref_result[0]) and not (hera_ref_result[0])):
                haplotype1_stuff.append((ref_location_list[index], 0, 0, 1, 0))
                print((ref_location_list[index], 0, 0, 1, 0))
        if haplotype_ref_alt[0] == 1:
            if ((hera_alt_result[0]) and not (stieg_alt_result[0])):
                haplotype1_stuff.append((ref_location_list[index], 0, 1, 0, 0))
                print((ref_location_list[index], 0, 1, 0, 0))
            elif ((stieg_alt_result[0]) and not (hera_alt_result[0])):
                haplotype1_stuff.append((ref_location_list[index], 0, 0, 0, 1))
                print((ref_location_list[index], 0, 0, 0, 1))
    write_path = "./intermediate/haplot1_result_block_{}.txt".format(phase_block_required)
    with open(write_path, 'a') as fw:
        for entry in haplotype1_stuff:
            fw.write("{}\n".format(entry))
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