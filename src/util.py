import os
# todo1
def run_fastk_make_intermediate_files(k, fast_k_loc, intermediate_loc, ref_loc):
    # get the file name from path
    file_name = ref_loc.split("/")[-1]
    intermediate_path = "{}{}".format(intermediate_loc, file_name)
    command_to_run = "{} -k{} -N\"{}\" -T64 {}".format(fast_k_loc, k, intermediate_path, ref_loc)
    # run the command
    print(command_to_run)
    print(os.popen(command_to_run).read())
    return

# todo2
def open_vcf_and_get_k_mer():

    return

# todo3
def search_for_kmer_in_intermediate():

    return

 #./FastK -k30 -N"./../" -T64 /data1/phasstphase_test/potato/reference/solTubHeraHap1.fa