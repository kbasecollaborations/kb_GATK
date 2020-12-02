import subprocess
import os
import logging
import gzip

class GATKUtils:

    def __init__(self):
        self.path = "/kb/module/deps"
        pass

    def run_cmd(self, cmd):
        """
        This function runs a third party command line tool
        eg. bgzip etc.
        :param command: command to be run
        :return: success
        """
        command = " ".join(cmd)
        print(command)
        logging.info("Running command " + command)
        cmdProcess = subprocess.Popen(command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT,
                                      shell=True)
        for line in cmdProcess.stdout:
            logging.info(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            logging.info('return code: ' + str(cmdProcess.returncode))
            if cmdProcess.returncode != 0:
                raise ValueError('Error in running command with return code: '
                                 + command
                                 + str(cmdProcess.returncode) + '\n')
        logging.info("command " + command + " ran successfully")
        return "success"

    def validate_params(self, params):
        '''

        :param params:
        :return:
        '''
        if 'assembly_or_genome_ref' not in params:
            raise ValueError('required assembly_or_genome_ref field was not defined')
        elif 'variation_object_name' not in params:
            raise ValueError('required variation_object_name field was not defined')
        elif 'alignment_ref' not in params:
            raise ValueError('required alignment_ref field was not defined')
        elif 'snp_filter' not in params:
            raise ValueError('required snp_filter field was not defined')
        elif 'indel_filter' not in params:
            raise ValueError('required indel_filter field was not defined')

    def build_genome(self, assembly_file):
        '''
        build_genome: This function build genome from assembly data.
        :param assembly_file:
        :return:
        '''
        command = ["bwa"]
        command.append("index")
        #command.extend(["-a", "bwtsw"])   #Todo : need to figure out usage of bwtsw option. 
        command.append(assembly_file)
        self.run_cmd(command)

    def index_assembly(self, assembly_file):
        '''
        index_assembly: index assembly
        :param assembly_file:
        :return:
        '''
        command = ["samtools"]
        command.append("faidx")
        command .append(assembly_file)
        self.run_cmd(command)

    def generate_sequence_dictionary(self, assembly_file):
        '''
        generate_sequence_dictionary : It generates sequence dictionary for assembly input.
        :param assembly_file:
        :return:
        '''
        output_dict = assembly_file.replace("fa", "dict")
        command = ["java"]
        command.append("-jar")
        command.extend([os.path.join(self.path, "picard.jar"), "CreateSequenceDictionary"])
        command.extend(["REFERENCE=", assembly_file])
        command.extend(["Output=", output_dict])
        self.run_cmd(command)

    def mapping_genome(self, ref_genome, rev_fastq, fwd_fastq, output_dir, strain_info):
        '''
        mapping genome: it maps genome pair end fastq reads to reference genome.
        :param ref_genome:
        :param rev_fastq:
        :param fwd_fastq:
        :param output_dir:
        :param strain_info:
        :return:
        '''
        command = ["bwa"]
        command.append("mem")
        command.extend(["-t", "32"])
        command.extend(["-M", "-R"])
        command.append("\"@RG\\tID:"+ strain_info + "\\tLB:" + strain_info + "\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:" + strain_info + "\" ")
        command.append(ref_genome)
        command.append(rev_fastq)
        command.append(fwd_fastq)
        command.extend([">",os.path.join(output_dir, "aligned_reads.sam")])
        self.run_cmd(command)

    def duplicate_marking(self, output_dir, sam_file):
        '''
        duplicate_marking: locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA
        :param output_dir:
        :param sam_file:
        :return:
        '''
        command = ["java"]
        command.append("-Xmx4G")
        command.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("MarkDuplicatesSpark") 
        command.extend(["-I", sam_file])
        command.extend(["-M", os.path.join(output_dir, "dedup_metrics.txt")])
        command.extend(["-O", os.path.join(output_dir, "sorted_dedup_reads.bam")])
        self.run_cmd(command)

    def sort_bam_index(self, output_dir):
        '''
        sort_bam_index: sort bam index file
        :param output_dir:
        :return:
        '''
        command = ["samtools"]
        command.append("index")
        command.append(os.path.join(output_dir, "aligned_reads.bam"))
        self.run_cmd(command)

    def collect_alignment_and_insert_size_metrics(self, assembly_file, output_dir):
        '''
        collect_alignment_and_insert_size_metrics :
        :param assembly_file:
        :param output_dir:
        :return:
        '''
        command1 = ["java"]
        command1.extend(["-jar", os.path.join(self.path, "picard.jar")])
        command1.append("CollectAlignmentSummaryMetrics")
        command1.extend(["R=", assembly_file])
        command1.extend(["I=", os.path.join(output_dir,"sorted_dedup_reads.bam")])
        command1.extend(["O=", os.path.join(output_dir,"alignment_metrics.txt")])
        self.run_cmd(command1)

        command2 = ["java"]
        command2.extend(["-jar", os.path.join(self.path, "picard.jar")])
        command2.extend(["CollectInsertSizeMetrics INPUT=", os.path.join(output_dir, "sorted_dedup_reads.bam")])
        command2.extend(["OUTPUT=", os.path.join(output_dir,"insert_metrics.txt" )])
        command2.extend(["HISTOGRAM_FILE=" + os.path.join(output_dir, "insert_size_histogram.pdf")])
        self.run_cmd(command2)

        command3 = ["samtools"]
        command3.append("depth")
        command3.extend(["-a", os.path.join(output_dir, "sorted_dedup_reads.bam")])
        command3.extend([">", os.path.join(output_dir, "depth_out.txt")])
        self.run_cmd(command3)

    def variant_calling(self, assembly_file, output_dir):
        '''
        variant_calling: The variants identified in this step will be filtered and provided as input for Base Quality Score Recalibration (BQSR)
        :param assembly_file:
        :param output_dir:
        :return:
        '''
        command = ["java"]
        command.append("-Xmx4G")
        command.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("HaplotypeCaller")
        command.extend(["-R", assembly_file])
        command.extend(["-I", os.path.join(output_dir, "sorted_dedup_reads.bam")])
        command.extend(["-O", os.path.join(output_dir, "raw_variants.vcf")])
        self.run_cmd(command)

    def extract_variants(self, assembly_file, output_dir):
        '''
        extract_variants: This step separates SNPs and Indels so they can be processed and used independently
        :param assembly_file:
        :param output_dir:
        :return:
        '''
        command1 = ["java"]
        command1.append("-Xmx4G")
        command1.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command1.append("SelectVariants")
        command1.extend(["-R", assembly_file])
        command1.extend(["-V", os.path.join(output_dir, "raw_variants.vcf")])
        command1.extend(["--select-type", "SNP"])
        command1.extend(["-O", os.path.join(output_dir, "raw_snps.vcf")])
        self.run_cmd(command1)

        command2 = ["java"]
        command2.append("-Xmx4G")
        command2.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command2.append("SelectVariants")
        command2.extend(["-R", assembly_file])
        command2.extend(["-V", os.path.join(output_dir, "raw_variants.vcf")])
        command2.extend(["--select-type", "INDEL"])
        command2.extend(["-O", os.path.join(output_dir, "raw_indels.vcf")])
        self.run_cmd(command2)


    def filter_SNPs(self, assembly_file, output_file, output_dir, params):
        '''
        filter_SNPs: QD < 2.0: This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.
                     FS > 60.0: This is the Phred-scaled probability that there is strand bias at the site.
                     MQ < 40.0: This is the root mean square mapping quality over all the reads at the site.
                     SOR > 4.0: This is another way to estimate strand bias using a test similar to the symmetric odds ratio test.
                     MQRankSum < -8.0: Compares the mapping qualities of the reads supporting the reference allele and the alternate allele.
                     ReadPosRankSum < -8.0: Compares whether the positions of the reference and alternate alleles are different within the reads.
        :param assembly_file:
        :param output_file:
        :param output_dir:
        :param params:
        :return:
        '''
        command  = ["java"]
        command.append("-Xmx4G")
        command.extend(["-jar", os.path.join(self.path,"gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("VariantFiltration")
        command.extend(["-R", assembly_file])
        command.extend(["-V", os.path.join(output_dir, "raw_snps.vcf")])
        command.extend(["-O", os.path.join(output_dir, output_file)])
        command.extend(["-filter-name", "\"QD_filter\"", "-filter", "\"QD", "<", params['snp_filter']['snp_qd_filter'] + "\""])
        command.extend(["-filter-name", "\"FS_filter\"", "-filter", "\"FS", "<", params['snp_filter']['snp_fs_filter'] + "\""])
        command.extend(["-filter-name", "\"MQ_filter\"", "-filter", "\"MQ", "<", params['snp_filter']['snp_mq_filter'] + "\""])
        command.extend(["-filter-name", "\"SOR_filter\"", "-filter", "\"SOR", "<", params['snp_filter']['snp_sor_filter'] + "\""])
        command.extend(["-filter-name", "\"MQRankSum_filter\"", "-filter", "\"MQRankSum", "<", params['snp_filter']['snp_mqrankSum_filter'] + "\""])
        command.extend(["-filter-name", "\"ReadPosRankSum_filter\"", "-filter", "\"ReadPosRankSum", "<", params['snp_filter']['snp_readposranksum_filter'] + "\""])
        self.run_cmd(command)

    def filter_Indels(self, assembly_file, output_file, output_dir, params):
        '''
        filter_Indels:
                    QD < 2.0: This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples.
                    FS > 200.0: This is the Phred-scaled probability that there is strand bias at the site.
                    SOR > 10.0: This is another way to estimate strand bias using a test similar to the symmetric odds ratio test.
        :param assembly_file:
        :param output_file:
        :param output_dir:
        :param params:
        :return:
        '''
        command  = ["java"]
        command.append("-Xmx4G")
        command.extend(["-jar", os.path.join(self.path,"gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("VariantFiltration")
        command.extend(["-R", assembly_file])
        command.extend(["-V", os.path.join(output_dir, "raw_indels.vcf")])
        command.extend(["-O", os.path.join(output_dir, output_file)])
        command.extend(["-filter-name", "QD_filter", "-filter", "\"QD", "<", params['indel_filter']['indel_qd_filter'] + "\""])
        command.extend(["-filter-name", "FS_filter", "-filter", "\"FS", "<", params['indel_filter']['indel_fs_filter'] + "\""])
        command.extend(["-filter-name", "SOR_filter", "-filter", "\"SOR", "<", params['indel_filter']['indel_sor_filter'] + "\""])
        self.run_cmd(command)

    def exclude_filtered_variants(self, output_dir):
        '''
        exclude_filtered_variants: extract only the passing variants and provide this as input to BQSR
        :param output_dir:
        :return:
        '''
        command1 = ["java"]
        command1.append("-Xmx4G")
        command1.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command1.extend(["SelectVariants", "--exclude-filtered"])
        command1.extend(["-V", os.path.join(output_dir, "filtered_snps.vcf")])
        command1.extend(["-O", os.path.join(output_dir, "bqsr_snps.vcf")])
        self.run_cmd(command1)

        command2 = ["java"]
        command2.append("-Xmx4G")
        command2.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command2.extend(["SelectVariants", "--exclude-filtered"])
        command2.extend(["-V", os.path.join(output_dir, "filtered_indels.vcf")])
        command2.extend(["-O", os.path.join(output_dir, "bqsr_indels.vcf")])
        self.run_cmd(command2)

    def base_quality_score_recalibration(self, assembly_file, data_table, output_dir):
        '''
        base_quality_score_recalibration: BQSR is performed twice. The second pass is optional, only required to produce a recalibration report.
        :param assembly_file:
        :param data_table:
        :param output_dir:
        :return:
        '''
        command = ["java"]
        command.append("-Xmx4G")
        command.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("BaseRecalibrator")
        command.extend(["-R", assembly_file])
        command.extend(["-I", os.path.join(output_dir, "sorted_dedup_reads.bam")])
        command.extend(["-known-sites", os.path.join(output_dir, "bqsr_snps.vcf")])
        command.extend(["-known-sites", os.path.join(output_dir, "bqsr_indels.vcf")])
        command.extend(["-O", os.path.join(output_dir, data_table)])
        self.run_cmd(command)
    
    def apply_BQSR(self, assembly_file, data_table, output_dir):
        '''
        apply_BQSR: This step applies the recalibration computed in the first BQSR step to the bam file.
        :param assembly_file:
        :param data_table:
        :param output_dir:
        :return:
        '''
        command = ["java"]
        command.append("-Xmx4G")
        command.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("ApplyBQSR")
        command.extend(["-R", assembly_file])
        command.extend(["-I", os.path.join(output_dir, "sorted_dedup_reads.bam")])
        command.extend(["-bqsr",os.path.join(output_dir, data_table)])
        command.extend(["-O", os.path.join(output_dir, "recal_reads.bam")])
        self.run_cmd(command)

    def analyze_covariates(self, output_dir):
        '''
        analyze_covariates: This step produces a recalibration report based on the output from the two BQSR runs.
        :param output_dir:
        :return:
        '''
        command = ["java"]
        command.append("-Xmx4G")
        command.extend(["-jar", os.path.join(self.path,"gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("AnalyzeCovariates")
        command.extend(["-before", os.path.join(output_dir, "recal_data.table")])
        command.extend(["-after", os.path.join(output_dir, "post_recal_data.table")])
        command.extend(["-plots", os.path.join(output_dir, "recalibration_plots.pdf")])
        self.run_cmd(command)

    def bgzip_vcf_file(self, filepath):
        '''
        bgzip_vcf_file: This step zip the vcf file.
        :param filepath:
        :return:
        '''
        bzfilepath = filepath + ".gz"
        command = ["bgzip", filepath]
        self.run_cmd(command)
        return bzfilepath        

    def index_vcf_file(self, filepath):
        '''
        index_vcf_file: this step index vcf file.
        :param filepath:
        :return:
        '''
        bzfilepath = self.bgzip_vcf_file(filepath)
        command  = ["tabix", "-p", "vcf", bzfilepath]
        self.run_cmd(command)
        return bzfilepath       

    def reheader(self, filepath, strain_info):
        '''
        :param filepath:
        :param strain_info:
        :return:
        '''
        reheader_vcf_path = filepath.replace(".vcf.gz", "") + "_reheader.vcf.gz"
        new_header_path = filepath.replace(".vcf.gz", "") + "_newheader.vcf"
        header_path = filepath.replace(".vcf.gz", "") + "_header.vcf"
        header_cmd = ["tabix -H", filepath, ">", header_path]
        self.run_cmd(header_cmd)

        pattern = 's/sample_1/' + strain_info + '/'
        replace_cmd = ['sed' , pattern, header_path, ">",  new_header_path]
        self.run_cmd(replace_cmd)

        command = ["tabix -r", new_header_path, filepath, ">", reheader_vcf_path]
        self.run_cmd(command)

        return reheader_vcf_path
