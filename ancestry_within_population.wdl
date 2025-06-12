version 1.0 

#IMPORTS
## According to this: https://cromwell.readthedocs.io/en/stable/Imports/ we can import raw from github
## so we can make use of the already written WDLs provided by WARP/VUMC Biostatistics

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/genotype/Utils.wdl" as http_GenotypeUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/agd/AgdUtils.wdl" as http_AgdUtils


# WORKFLOW

workflow agd_ancestry_workflow{
    input{
        # workflow choices 

        Boolean external_spike_in = true 
        Boolean scope_supervised = true
        Boolean run_pca = true
        Boolean run_scope= true

        # required inputs original data as array of chromosomes 

        Array[File] source_pgen_files
        Array[File] source_pvar_files
        Array[File] source_psam_files
        Array[String] chromosomes
        String target_prefix

        File id_map_file

        # optional outputs for exporting 
        String? project_id
        String target_gcp_folder

        #inputs for subsetting ancestries
        File ancestry_id_file
        String ancestry_set
        String ancestry_column
        String iid_column = "GRID"

        # optional inputs for spike in data - required if merging spike in data 

        File? spike_in_pgen_file
        File? spike_in_pvar_file
        File? spike_in_psam_file
        File? spike_in_relatives_exclude 
    
        # optional inputs for supervised scope    - required if running supervised scope 

        File? supervised_scope_reference_pgen_file
        File? supervised_scope_reference_pvar_file
        File? supervised_scope_reference_psam_file
        File? supervised_scope_reference_superpop_file   
        File? supervised_scope_reference_relatives_exclude
        File? supervised_scope_reference_freq

        #optional inputs for unsupervised scope - required if running unsupervised scope
        String? scope_plink2_maf_filter = "--maf 0.01"

        String? scope_plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
        File? scope_long_range_ld_file

        String? pca_plink2_LD_filter_option = "--indep-pairwise 50 5 0.2"
        String? pca_plink2_maf_filter = "--maf 0.05"
        Int? K_supervised = 4
        Int? K_unsupervised =4
        Int? seed = 1234
    }

    # If the user chose to use the supervised scope and there is no precalculated reference allele frequency provided, allele frequency can be calculated from provided plink files instead
    if(scope_supervised && !defined(supervised_scope_reference_freq)){
        File required_supervised_scope_reference_pgen_file = select_first([supervised_scope_reference_pgen_file])
        File required_supervised_scope_reference_pvar_file = select_first([supervised_scope_reference_pvar_file])
        File required_supervised_scope_reference_psam_file = select_first([supervised_scope_reference_psam_file])
        File required_supervised_scope_reference_superpop_file = select_first([supervised_scope_reference_superpop_file])

        call CalculateFreq{
        input: 
            pgen_file = required_supervised_scope_reference_pgen_file,
            pvar_file = required_supervised_scope_reference_pvar_file,
            psam_file = required_supervised_scope_reference_psam_file,
            superpop_file = required_supervised_scope_reference_superpop_file,
            relatives_exclude = supervised_scope_reference_relatives_exclude
        }
        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFileFreq{
                input:
                    source_file = CalculateFreq.freq_file,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }
    }

    # If the user chose to use an external spike in, then first the spike-in data must be merged with the original data, and then the pipeline can proceed 
    if(external_spike_in){

        scatter (idx in range(length(chromosomes))) {
            String chromosome_for_spike_in = chromosomes[idx]
            File agd_pgen_file_for_spike_in = source_pgen_files[idx]
            File agd_pvar_file_for_spike_in = source_pvar_files[idx]
            File agd_psam_file_for_spike_in = source_psam_files[idx]

            call ConvertPgenToBed as ConvertPgenToBedForSpikeIn{
                input:
                    pgen = agd_pgen_file_for_spike_in,
                    pvar = agd_pvar_file_for_spike_in,
                    psam = agd_psam_file_for_spike_in, 
                    out_prefix = chromosome_for_spike_in
            }

            call SubsetChromosomeTGP{
                input: 
                    pgen_file = spike_in_pgen_file,
                    pvar_file = spike_in_pvar_file,
                    psam_file = spike_in_psam_file,
                    chromosome = chromosome_for_spike_in,
                    relatives_exclude = spike_in_relatives_exclude
            }

            call Merge1000genomesAGD {
                input:
                    agd_bed_file = ConvertPgenToBedForSpikeIn.convert_Pgen_out_bed,
                    agd_bim_file = ConvertPgenToBedForSpikeIn.convert_Pgen_out_bim,
                    agd_fam_file = ConvertPgenToBedForSpikeIn.convert_Pgen_out_fam,
                    TGP_bed_file = SubsetChromosomeTGP.subset_reference_out_bed_file,
                    TGP_bim_file = SubsetChromosomeTGP.subset_reference_out_bim_file,
                    TGP_fam_file = SubsetChromosomeTGP.subset_reference_out_fam_file,
                    chromosome = chromosome_for_spike_in
            }
        }
    }

    #now we can proceed with the pipeline, using select_first to get the first non empty array value to select either the merged or original input files 

    Array[File] my_pgen_files=select_first([Merge1000genomesAGD.merge_spike_out_pgen_file, source_pgen_files])
    Array[File] my_pvar_files=select_first([Merge1000genomesAGD.merge_spike_out_pvar_file, source_pvar_files])
    Array[File] my_psam_files=select_first([Merge1000genomesAGD.merge_spike_out_psam_file, source_psam_files])

    # now filter to just those individuals in the selected ancestry 
    scatter(idx in range(length(chromosomes))) {
        String chromosome_per_ancestry = chromosomes[idx]
        File pgen_file_per_ancestry = my_pgen_files[idx]
        File pvar_file_per_ancestry = my_pvar_files[idx]
        File psam_file_per_ancestry = my_psam_files[idx]

        #subset to just those individuals in the ancestry file using plink 
        call subset_pgen_by_ancestry{
            input:
                pgen_file = pgen_file_per_ancestry,
                pvar_file = pvar_file_per_ancestry,
                psam_file = psam_file_per_ancestry,
                ancestry_file = ancestry_id_file,
                target_ancestry = ancestry_set,
                ancestry_column = ancestry_column,
                iid_column = iid_column
        }
    }

    if(run_pca){
        scatter (idx in range(length(chromosomes))) {
            String chromosome_for_pca = chromosomes[idx]
            File pgen_file_for_pca = subset_pgen_by_ancestry.subset_pgen[idx]
            File pvar_file_for_pca = subset_pgen_by_ancestry.subset_pvar[idx]
            File psam_file_for_pca = subset_pgen_by_ancestry.subset_psam[idx]
            String replaced_sample_name_for_pca = "~{chromosome_for_pca}.psam"

            #remove high LD regions, LD prune, and maf filter 
             call PreparePlink as PreparePlinkPCA{
                input:
                    pgen_file = pgen_file_for_pca ,
                    pvar_file = pvar_file_for_pca,
                    psam_file = psam_file_for_pca,
                    long_range_ld_file = scope_long_range_ld_file,
                    plink2_maf_filter = pca_plink2_maf_filter,
                    plink2_LD_filter_option = pca_plink2_LD_filter_option,
                    chromosome = chromosome_for_pca
            }
        }

        call http_GenotypeUtils.MergePgenFiles as MergePgenFilesForPCA{
            input:
                pgen_files = PreparePlinkPCA.prepare_plink_unsupervised_output_pgen_file,
                pvar_files = PreparePlinkPCA.prepare_plink_unsupervised_output_pvar_file,
                psam_files = PreparePlinkPCA.prepare_plink_unsupervised_output_psam_file,
                output_prefix = target_prefix
        }

        #calculate IBD, remove related individuals, calculate PCs using plink & project related individuals on to PCs 
        call ibd_pca_project{
            input:
                pgen_file = MergePgenFilesForPCA.output_pgen, 
                pvar_file = MergePgenFilesForPCA.output_pvar,
                psam_file = MergePgenFilesForPCA.output_psam,  
                target_name = target_prefix   
        }

        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFile_PCAone {
                input:
                    source_file = ibd_pca_project.pca_merged,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }

        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFile_PCAtwo {
                input:
                    source_file = ibd_pca_project.eigenvalues,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }

         if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyOneFile as CopyFile_PCAthree {
                input:
                    source_file = ibd_pca_project.eigenvectors,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
        }
    }

    if(run_scope){
        scatter (idx in range(length(chromosomes))) {
            String chromosome_for_scope = chromosomes[idx]
            File pgen_file_for_scope =  subset_pgen_by_ancestry.subset_pgen[idx]
            File pvar_file_for_scope =  subset_pgen_by_ancestry.subset_pvar[idx]
            File psam_file_for_scope =  subset_pgen_by_ancestry.subset_psam[idx]
            String replaced_sample_name_for_scope = "~{chromosome_for_scope}.psam"

            #I think I need this to get the IDs correctly as GRIDS

            call http_AgdUtils.ReplaceICAIdWithGrid as ReplaceICAIdWithGridForScope {
                input:
                    input_psam = psam_file_for_scope,
                    id_map_file = id_map_file,
                    output_psam = replaced_sample_name_for_scope
            }
            call PreparePlink as PreparePlink{
                input:
                    pgen_file = pgen_file_for_scope,
                    pvar_file = pvar_file_for_scope,
                    psam_file = ReplaceICAIdWithGridForScope.output_psam,
                    long_range_ld_file = scope_long_range_ld_file,
                    plink2_maf_filter = scope_plink2_maf_filter,
                    plink2_LD_filter_option = scope_plink2_LD_filter_option,
                    chromosome = chromosome_for_scope 
            }
        }
        call http_GenotypeUtils.MergePgenFiles as MergePgenFilesForScope{
            input:
                pgen_files = PreparePlink.prepare_plink_unsupervised_output_pgen_file,
                pvar_files = PreparePlink.prepare_plink_unsupervised_output_pvar_file,
                psam_files = PreparePlink.prepare_plink_unsupervised_output_psam_file,
                output_prefix = target_prefix
        }

        call ConvertPgenToBed as ConvertPgenToBedForScope{
            input: 
                pgen = MergePgenFilesForScope.output_pgen, 
                pvar = MergePgenFilesForScope.output_pvar,
                psam = MergePgenFilesForScope.output_psam, 
        }

        call RunScopeUnsupervised{    
            input:
                bed_file = ConvertPgenToBedForScope.convert_Pgen_out_bed,
                bim_file = ConvertPgenToBedForScope.convert_Pgen_out_bim,
                fam_file = ConvertPgenToBedForScope.convert_Pgen_out_fam,
                K = K_unsupervised,
                output_string = target_prefix,
                seed = seed
        }

        if(scope_supervised){
            call QCAllelesBim{
                input:
                    bim_file = ConvertPgenToBedForScope.convert_Pgen_out_bim,
                    freq_file = select_first([CalculateFreq.freq_file, supervised_scope_reference_freq])
            }

            call PreparePlinkSupervised{
                input:
                    bed_file = ConvertPgenToBedForScope.convert_Pgen_out_bed,
                    bim_file = ConvertPgenToBedForScope.convert_Pgen_out_bim,
                    fam_file = ConvertPgenToBedForScope.convert_Pgen_out_fam,
                    variant_list = QCAllelesBim.out_variants
            }

            call RunScopeSupervised{
                input:
                    bed_file = PreparePlinkSupervised.prepare_plink_supervised_out_bed,
                    bim_file = PreparePlinkSupervised.prepare_plink_supervised_out_bim,
                    fam_file = PreparePlinkSupervised.prepare_plink_supervised_out_fam,
                    K = K_supervised,
                    output_string = target_prefix,
                    seed = seed,
                    topmed_freq = QCAllelesBim.out_frq
            }
        }
        
        if(defined(target_gcp_folder)){
            call http_GcpUtils.MoveOrCopyThreeFiles as CopyFiles_one{
                input:
                    source_file1 = RunScopeUnsupervised.outP,
                    source_file2 = RunScopeUnsupervised.outQ,
                    source_file3 = RunScopeUnsupervised.outV,
                    is_move_file = false,
                    project_id = project_id,
                    target_gcp_folder = select_first([target_gcp_folder])
            }
            if(scope_supervised){
                call http_GcpUtils.MoveOrCopyThreeFiles as CopyFiles_two {
                    input:
                        source_file1 = select_first([RunScopeSupervised.outP]),
                        source_file2 = select_first([RunScopeSupervised.outQ]),
                        source_file3 = select_first([RunScopeSupervised.outV]),
                        is_move_file = false,
                        project_id = project_id,
                        target_gcp_folder = select_first([target_gcp_folder])
                }
            }
        }
    }

 output {
        Array[File] ancestry_outputs = select_all([CopyFileFreq.output_file, CopyFile_PCAone.output_file, CopyFile_PCAtwo.output_file, CopyFiles_one.output_file1, CopyFiles_one.output_file2, CopyFiles_one.output_file3, CopyFiles_two.output_file1, CopyFiles_two.output_file2, CopyFiles_two.output_file3])
    }
}


# TASKS

## for frequency calculation 
task CalculateFreq{
    input{ 
        File pgen_file
        File pvar_file
        File psam_file
        File superpop_file
        File? relatives_exclude

        String? plink2_maf_filter = "--maf 0.001"

        Int memory_gb = 20
        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

    Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 4) + 20
    String out_name =  basename(pgen_file, ".pgen") + "_within_superpop_freqs"
    String out_file = out_name + ".frq.strat"

    command <<<
        # take the TGP data, remove duplicates, restrict to biallelic SNPs (necessary since cannot calculate frq within in plink2 in the same way as in plink1), and calculate within super populations 
        # include a MAF filter (default 0.001) to reduce size of output file
        plink2 --pgen ~{pgen_file} \
            --pvar ~{pvar_file} \
            --psam ~{psam_file} \
            --allow-extra-chr \
            --chr 1-22, X, Y \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 10 truncate \
            --max-alleles 2 \
            ~{plink2_maf_filter} \
            --rm-dup 'exclude-all' \
            --remove ~{relatives_exclude} \
            --make-bed \
            --out tgp_nodup  

        plink --bfile tgp_nodup \
            --allow-extra-chr \
            --freq --within ~{superpop_file} \
            --out ~{out_name}
    >>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output{
        File freq_file = out_file
    }
}

## for merging with spike-in data task 
task SubsetChromosomeTGP {
    input {
        File? pgen_file
        File? pvar_file
        File? psam_file
        String chromosome
        File? relatives_exclude
        Int? memory_gb = 20
        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }
    
    String out_string = "TGP_" + chromosome
    String in_chromosome = sub(chromosome, "chr", "")
    Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB") * 2) * 2 + 20
    
    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }
    
    command {
         plink2 \
             --pgen ~{pgen_file} --pvar ~{pvar_file} --psam ~{psam_file} \
             --allow-extra-chr \
             --chr ~{in_chromosome} \
             --remove ~{relatives_exclude} \
             --make-bed \
             --out ~{out_string}
    }
    
    output {
        File subset_reference_out_bed_file = out_string + ".bed"
        File subset_reference_out_bim_file = out_string + ".bim"
        File subset_reference_out_fam_file = out_string + ".fam"
    }
}


task ConvertPgenToBed{
    input {
        File pgen 
        File pvar 
        File psam 

        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

        String? out_prefix

        Int? memory_gb = 20

        Int? disk_size_multiplier = 4
        Int? disk_size_addition = 20
    }

    Int disk_size = ceil(size([pgen, pvar, psam], "GB"))*disk_size_multiplier + disk_size_addition

    String out_string = if defined(out_prefix) then out_prefix else basename(pgen, ".pgen")

    command {
        plink2 \
            --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            --make-bed \
            --out ~{out_string}
    }

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File convert_Pgen_out_bed = "${out_string}.bed"
        File convert_Pgen_out_bim = "${out_string}.bim"
        File convert_Pgen_out_fam = "${out_string}.fam"
    }
}

task Merge1000genomesAGD{
    input{
        File agd_bed_file
        File agd_bim_file
        File agd_fam_file

        File TGP_bed_file
        File TGP_bim_file
        File TGP_fam_file

        String? chromosome

        Int? memory_gb = 20


        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

    Int disk_size = ceil(size([agd_bed_file, TGP_bed_file], "GB")  * 2)*2 + 20


    String out_string = "AGD_TGP_" + chromosome
    String agd_prefix = basename(agd_bed_file, ".bed")
    String TGP_prefix = basename(TGP_bed_file, ".bed")
    String agd_prefix_rename= agd_prefix + "_renamed"
    String agd_prefix_2 = agd_prefix + "_2"
    String TGP_prefix_2 = TGP_prefix + "_2"

    String relocated_bed = agd_prefix + ".bed"
    String relocated_bim = agd_prefix + ".bim"
    String relocated_fam = agd_prefix + ".fam"

    String relocated_tgp_bed = TGP_prefix + ".bed"
    String relocated_tgp_bim = TGP_prefix + ".bim"
    String relocated_tgp_fam = TGP_prefix + ".fam"

    String agd_snp_list = agd_prefix + "_renamed.snplist"
    String TGP_snp_list = TGP_prefix + ".snplist"


    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    command{

        ln ~{agd_bed_file} ./~{relocated_bed}
        ln ~{agd_bim_file} ./~{relocated_bim}
        ln ~{agd_fam_file} ./~{relocated_fam}

        ln ~{TGP_bed_file} ./~{relocated_tgp_bed}
        ln ~{TGP_bim_file} ./~{relocated_tgp_bim}
        ln ~{TGP_fam_file} ./~{relocated_tgp_fam}

        plink2 \
            --bfile ~{agd_prefix} \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 20 truncate \
            --make-bed \
            --out ~{agd_prefix_rename}

        plink2 \
            --bfile ~{agd_prefix_rename} \
            --rm-dup force-first 'list' \
            --write-snplist \
            --out ~{agd_prefix_rename} 
        
        plink2 \
            --bfile ~{TGP_prefix} \
            --write-snplist \
            --out ~{TGP_prefix}

        plink \
            --bfile ~{agd_prefix_rename} \
            --bmerge ~{TGP_prefix} \
            --make-bed \
            --out merged_beds_files

        plink2 \
            --bfile ~{agd_prefix_rename} \
            --exclude merged_beds_files-merge.missnp \
            --extract-intersect ~{agd_snp_list} ~{TGP_snp_list} \
            --make-bed \
            --out ~{agd_prefix_2}

        plink2 \
            --bfile ~{TGP_prefix} \
            --exclude merged_beds_files-merge.missnp \
            --extract-intersect ~{agd_snp_list} ~{TGP_snp_list} \
            --make-bed \
            --out ~{TGP_prefix_2}

        plink \
            --bfile ~{agd_prefix_2} \
            --bmerge ~{TGP_prefix_2} \
            --make-bed \
            --out merged_beds_files2
        
        plink2 \
            --bfile merged_beds_files2 \
            --make-pgen \
            --out ~{out_string}
    }

    output{
        File merge_spike_out_pgen_file = out_string + ".pgen"
        File merge_spike_out_pvar_file = out_string + ".pvar"
        File merge_spike_out_psam_file = out_string + ".psam"
    }
}

## for PCA 
task ExtractVariants{
  input {
    File? pgen_file
    File? pvar_file
    File? psam_file 

    String? chromosome

    File? variants_extract_file

    Int? memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_pvar = chromosome + ".pvar"
  String new_psam = chromosome + ".psam"
  String intermediate_pgen = chromosome + "_varids.pgen"
  String intermediate_pvar = chromosome + "_varids.pvar"
  String intermediate_psam = chromosome + "_varids.psam"

  command {
    plink2 \
      --pgen ~{pgen_file} \
      --pvar ~{pvar_file} \
      --psam ~{psam_file} \
      --snps-only \
      --set-all-var-ids chr@:#:\$r:\$a \
      --new-id-max-allele-len 1000 \
      --make-pgen \
      --out ~{chromosome}_varids
    
    plink2 \
      --pgen ~{intermediate_pgen} \
      --pvar ~{intermediate_pvar} \
      --psam ~{intermediate_psam} \
      --extract ~{variants_extract_file} \
      --make-pgen \
      --out ~{chromosome}
  }

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }

  output {
    File extract_variants_output_pgen_file = new_pgen
    File extract_variants_output_pvar_file = new_pvar
    File extract_variants_output_psam_file = new_psam
  }

}

task ProjectPCA{
  input{
    File? pgen_file
    File? pvar_file
    File? psam_file
    File? PCA_loadings
    File? PCA_AF
    String? OUTNAME

    Int memory_gb = 20
    Int cpu = 8

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

  }

  Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB")  * 2) + 20

  String pca_file = OUTNAME + ".genotype.pca.sscore"
  String pca_variants = OUTNAME + "genotype.sscore.vars"

  command {
    plink2 --pgen ~{pgen_file} --pvar ~{pvar_file} --psam ~{psam_file} --score ~{PCA_loadings} \
    variance-standardize \
    cols=-scoreavgs,+scoresums \
    list-variants \
    header-read \
    --score-col-nums 3-22 \
    --read-freq ~{PCA_AF} \
    --out ~{OUTNAME}

    cp ~{OUTNAME}.sscore ~{pca_file}
    cp ~{OUTNAME}.sscore.vars ~{pca_variants}
    }

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  
  output {
    File output_pca_file = "~{pca_file}"
    File output_pca_variants="~{pca_variants}"
  }
}

## for scope
task PreparePlink{
  input {
    File? pgen_file
    File? pvar_file
    File? psam_file 

    String? chromosome

    String? plink2_maf_filter = "--maf 0.01"
    String? plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
    File? long_range_ld_file


    Int memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_pvar = chromosome + ".pvar"
  String new_psam = chromosome + ".psam"
  String out_prefix = chromosome 


  command {
    plink2 \
      --pgen ~{pgen_file} \
      --pvar ~{pvar_file} \
      --psam ~{psam_file} \
      ~{plink2_maf_filter} \
      --snps-only \
      --const-fid \
      --set-all-var-ids chr@:#:\$r:\$a \
      --new-id-max-allele-len 1000 \
      --make-pgen \
      --out maf_filtered
      
    plink2 \
        --pgen maf_filtered.pgen \
        --pvar maf_filtered.pvar \
        --psam maf_filtered.psam \
        --exclude range ~{long_range_ld_file} \
        --make-pgen \
        --out maf_filtered_longrange
    
    plink2 \
      --pgen maf_filtered_longrange.pgen \
      --pvar maf_filtered_longrange.pvar \
      --psam maf_filtered_longrange.psam \
      ~{plink2_LD_filter_option}

    plink2 \
        --pgen maf_filtered_longrange.pgen \
        --pvar maf_filtered_longrange.pvar \
        --psam maf_filtered_longrange.psam \
        --extract plink2.prune.in \
        --make-pgen \
        --out ~{out_prefix}
  }

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }

  output {
    File prepare_plink_unsupervised_output_pgen_file = new_pgen
    File prepare_plink_unsupervised_output_pvar_file = new_pvar
    File prepare_plink_unsupervised_output_psam_file = new_psam
  }
}

task QCAllelesBim{
    input {
        File? bim_file
        File? freq_file

        String docker = "blosteinf/r_utils_terra:0.1"
        Int memory_gb = 20
    }

    Int disk_size = ceil(size([bim_file, freq_file], "GB")  * 2) + 20

    command {
        ls /home/r-environment/
        Rscript /home/r-environment/allele_qc.R --in_freq ~{freq_file} --in_bim ~{bim_file} 
    }

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File out_frq = "corrected_freq.frq"
        File out_variants = "variants_to_extract.txt"
    }
}

task PreparePlinkSupervised{
    input { 
        File? bed_file
        File? bim_file
        File? fam_file 

        File? variant_list 
        String? out_string = "variant_filtered"

        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

  Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20
  Int memory_gb = 20

  String new_bed = out_string + ".bed"
  String new_bim = out_string + ".bim"
  String new_fam= out_string + ".fam"

  command { 
    plink2 \
        --bed ~{bed_file} \
        --bim ~{bim_file} \
        --fam ~{fam_file} \
        --extract ~{variant_list} \
        --make-bed \
        --out ~{out_string}
  }

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File prepare_plink_supervised_out_bed = new_bed
        File prepare_plink_supervised_out_bim = new_bim
        File prepare_plink_supervised_out_fam = new_fam
        String out_prefix = out_string
    }

}

task RunScopeUnsupervised{
    input{

        File bed_file
        File bim_file
        File fam_file

        Int? K
        String output_string
        Int? seed

        Int memory_gb = 60
        String docker = "blosteinf/scope:0.1"
    }

    String plink_binary_prefix =  basename(bed_file, ".bed")
    String relocated_bed = plink_binary_prefix + ".bed"
    String relocated_bim = plink_binary_prefix + ".bim"
    String relocated_fam = plink_binary_prefix + ".fam"

    String unsup_output = output_string + "_unsupervised_" 
    Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20

    command <<<
        ln ~{bed_file} ./~{relocated_bed}
        ln ~{bim_file} ./~{relocated_bim}
        ln ~{fam_file} ./~{relocated_fam}
        scope -g ~{plink_binary_prefix} -k ~{K} -seed ~{seed} -o ~{unsup_output}
        awk '{ for (i=1; i<=NF; i++) { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++) { str=str" "a[i,j]; } print str } }' ~{unsup_output}Qhat.txt > transposed_Qhat.txt
        cut -f2 ./~{relocated_fam} | paste - transposed_Qhat.txt > ~{unsup_output}Qhat.txt
    >>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
  }

    output {
        File outP= "${unsup_output}Phat.txt"
        File outQ= "${unsup_output}Qhat.txt"
        File outV= "${unsup_output}V.txt"
    }
}

task RunScopeSupervised{
    input{
       
        File bed_file
        File bim_file
        File fam_file

        Int? K
        String output_string
        Int? seed

        File? topmed_freq

        Int memory_gb = 60

        String docker = "blosteinf/scope:0.1"
    }

    String plink_binary_prefix = basename(bed_file, ".bed")
    String relocated_bed= plink_binary_prefix + ".bed"
    String relocated_bim= plink_binary_prefix + ".bim"
    String relocated_fam= plink_binary_prefix + ".fam"

    String sup_output = output_string + "_supervised_"

    Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20

    command <<<
        ln ~{bed_file} ./~{relocated_bed}
        ln ~{bim_file} ./~{relocated_bim}
        ln ~{fam_file} ./~{relocated_fam}
        scope -g ~{plink_binary_prefix} -freq ~{topmed_freq} -k ~{K} -seed ~{seed} -o ~{sup_output}
        ls
        awk '{ for (i=1; i<=NF; i++) { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++) { str=str" "a[i,j]; } print str } }' ~{sup_output}Qhat.txt > transposed_Qhat.txt
        cut -f2 ./~{relocated_fam} | paste - transposed_Qhat.txt > ~{sup_output}Qhat.txt
    >>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
  }

    output {
        File outP= "${sup_output}Phat.txt"
        File outQ= "${sup_output}Qhat.txt"
        File outV= "${sup_output}V.txt"
    }
}

task subset_pgen_by_ancestry {
  input {
    File pgen_file
    File psam_file
    File pvar_file
    File ancestry_file
    String target_ancestry
    String ancestry_column
    String iid_column

    Int memory_gb = 20
    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

    Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 2) + 20

    command <<<
    set -euo pipefail

    ln -s ~{pgen_file} base.pgen
    ln -s ~{psam_file} base.psam
    ln -s ~{pvar_file} base.pvar

    # Extract base filename without extension
    base_name=$(basename ~{pgen_file} .pgen)
    output_prefix="~{target_ancestry}_${base_name}"

    # Get column indices
    ancestry_col=$(head -n 1 ~{ancestry_file} | tr '\t' '\n' | awk -v col="~{ancestry_column}" '{ if ($1 == col) print NR; }')
    iid_col=$(head -n 1 ~{ancestry_file} | tr '\t' '\n' | awk -v col="~{iid_column}" '{ if ($1 == col) print NR; }')

    if [ -z "$ancestry_col" ] || [ -z "$iid_col" ]; then
      echo "Required column not found in ancestry_file" >&2
      exit 1
    fi

    # Create keep file with FID=0 and matching IID
    awk -v ac="$ancestry_col" -v ic="$iid_col" -v anc="~{target_ancestry}" \
      'BEGIN {FS=OFS="\t"} NR > 1 && $ac == anc {print $ic}' ~{ancestry_file} > keep_ids.txt

    awk 'NR==FNR { keep[$1]; next } $2 in keep { print $1, $2 }' OFS='\t' keep_ids.txt base.psam > keep_ids_fids.txt

    # Subset with plink2    
    plink2 \
      --pfile base \
      --keep keep_ids_fids.txt \
      --make-pgen \
      --out "$output_prefix"
  >>>

  output {
    File subset_pgen = "~{target_ancestry}_~{basename(pgen_file, '.pgen')}.pgen"
    File subset_psam = "~{target_ancestry}_~{basename(pgen_file, '.pgen')}.psam"
    File subset_pvar = "~{target_ancestry}_~{basename(pgen_file, '.pgen')}.pvar"
  }

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
}

task ibd_pca_project {
  input {
    File pgen_file
    File pvar_file
    File psam_file
    Int n_pcs = 20
    String target_name
    Float PI_HAT = 0.2
    Int my_n_splits = 100

    Int memory_gb = 20
    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    Int? disk_size_multiplier = 4
    Int? disk_size_addition = 20
  
  }
  
  Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB") * disk_size_multiplier) + disk_size_addition
  
  command <<<
    set -euo pipefail

    # Link input files
    ln -s ~{pgen_file} input.pgen
    ln -s ~{pvar_file} input.pvar
    ln -s ~{psam_file} input.psam

    target_name=~{target_name}
    threshold=~{PI_HAT}
    n_splits=~{my_n_splits}
    mem_for_plink=$(((~{memory_gb}-3) * 1024))  # Give buffer and convert GB to MB

    # Step 1: Convert to .bed format for PLINK1.9
    plink2 --pfile input --make-bed --out ${target_name}_step1

    # Step 2a: Estimate IBD
    for i in $(seq 1 $n_splits); do
    plink --bfile ${target_name}_step1 \
          --genome \
          --parallel ${i} ${n_splits} \
          --out ${target_name}_ibd_part${i} \
          --memory ${mem_for_plink} \
          --min ${threshold}
    done

    # Step 2b: Concatenate all *.genome files from the parallel runs
    head -n 1 ${target_name}_ibd_part1.genome > ${target_name}_ibd_all.genome  # Header
    for i in $(seq 1 $n_splits); do
        tail -n +2 ${target_name}_ibd_part${i}.genome >> ${target_name}_ibd_all.genome
    done


    # Step 3: create initial set of all individuals and extract related pairs
    genome_file="${target_name}_ibd_all.genome"
    awk 'NR > 1 { print $1, $2}' input.psam | sort | uniq > all_samples.txt
    awk -v thresh="$threshold" 'NR > 1 && $10 > thresh { print $1, $2, $3, $4 }' "$genome_file" > related_pairs.txt

    # Step 4: Greedily remove one idndividual from each related pair 
    > removed.txt  # empty file to store removed individuals

    while read fid1 iid1 fid2 iid2; do
        id1="$fid1 $iid1"
        id2="$fid2 $iid2"

        # If neither already removed, arbitrarily remove the second
        if ! grep -Fxq "$id1" removed.txt && ! grep -Fxq "$id2" removed.txt; then
            echo "$id2" >> removed.txt
        fi
    done < related_pairs.txt

    sort all_samples.txt > all_samples.sorted.txt
    sort removed.txt | uniq > removed.sorted.txt
    comm -23 all_samples.sorted.txt removed.sorted.txt > ${target_name}_unrelated.txt   
    echo "Count of unrelated samples:"
    wc -l ${target_name}_unrelated.txt

    # Step 5: PCA on unrelated
    plink2 \
      --pfile input \
      --keep ${target_name}_unrelated.txt \
      --freq counts \
      --ac-founders \
      --pca approx ~{n_pcs} allele-wts vcols=chrom,ref,alt  \
      --out ${target_name}_pca_unrelated \
      --memory ${mem_for_plink} 

    # Step 6: Project related
    start_col=3
    end_col=$((2 + ~{n_pcs}))

    #https://groups.google.com/g/plink2-users/c/W6DL5-hs_Q4/m/b_o3JMrxAwAJ
    #https://groups.google.com/g/plink2-users/c/ZO84YhMYabc

    #because it is suggested to compare 'apples to apples', I project all samples into the same space, rather than just projected unrelated samples
    end_col=$((5 + ~{n_pcs}))

    plink2 \
      --pfile input \
      --read-freq ${target_name}_pca_unrelated.acount \
      --score ${target_name}_pca_unrelated.eigenvec.allele 2 5 header-read variance-standardize \
      --score-col-nums 6-${end_col} \
      --out ${target_name}_projected_all_samples
      --memory ${mem_for_plink}

    # Step 7: Add a column indicating if the sample was included in the original PC 
    awk 'NR==FNR {ids[$1]; next} {print $0, ($1 in ids ? "related" : "unrelated")}' removed.sorted.txt ${target_name}_projected_all_samples.sscore > ${target_name}_pca_combined.tsv
 >>>

  output {
    File pca_merged = "${target_name}_pca_combined.tsv"
    File eigenvalues = "${target_name}_pca_unrelated.eigenval"
    File eigenvectors = "${target_name}_pca_unrelated.eigenvec"
  }

    runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
}