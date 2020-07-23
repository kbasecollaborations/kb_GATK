import subprocess
import os

class GATKUtils:

    def __init__(self):
        self.path = "/kb/module/deps"
        pass

    def run_cmd(self, command):
        cmd = " ".join(command)

        try:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if stdout:
                print("ret> ", process.returncode)
                print("OK> output ", stdout)
            if stderr:
                print("ret> ", process.returncode)
                print("Error> error ", stderr.strip())

        except OSError as e:
            print("OSError > ", e.errno)
            print("OSError > ", e.strerror)
            print("OSError > ", e.filename)

    def validate_params(self, params):
        if 'assembly_or_genome_reff' not in params:
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
        command = ["bwa"]
        command.append("index")
        command.extend(["-a", "bwts"])
        command.append(assembly_file)
        self.run_cmd(command)

    def index_assembly(self, assembly_file):
        command = ["samtools"]
        command.append("faidx")
        command .append(assembly_file)
        self.run_cmd(command)

    def generate_sequence_dictionary(self, assembly_file):
        output_dict = assembly_file.replace("fa", "dict")
        command = ["java"]
        command.append("-jar")
        command.extend([os.path.join(self.path, "picard.jar"), "CreateSequenceDictionary"])
        command.extend(["REFERENCE=", assembly_file])
        command.extend(["Output=", output_dict])
        self.run_cmd(command)

    def mapping_genome(self, ref_genome, rev_fastq, fwd_fastq, output_dir):
        command = ["bwa"]
        command.append("mem")
        command.extend(["-t", "32"])
        command.extend(["-M", "-R"])
        command.append("\"@RG\\tID:sample_1\\tLB:sample_1\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:sample_1\" ")
        command.append(ref_genome)
        command.append(rev_fastq)
        command.append(fwd_fastq)
        command.extend([">",os.path.join(output_dir, "aligned_reads.sam")])
        self.run_cmd(command)

    def duplicate_marking(self, output_dir, sam_file):
        command = ["java"]
        command.extend(["-jar", os.path.join(self.path, "picard.jar")])
        command.append("SortSam")
        command.extend(["INPUT=", sam_file])
        command.extend(["OUTPUT=", os.path.join(output_dir, "aligned_reads.bam")])
        command.append("SORT_ORDER=coordinate")
        self.run_cmd(command)

    def sort_bam_index(self, output_dir):
        command = ["samtools"]
        command.append("index")
        command.append(os.path.join(output_dir, "aligned_reads.bam"))
        self.run_cmd(command)

    def collect_alignment_and_insert_size_metrics(self, assembly_file, output_dir):
        command1 = ["java"]
        command1.extend(["-jar", os.path.join(self.path, "picard.jar")])
        command1.append("CollectInsertSizeMetrics INPUT=")
        command1.append(os.path.join(os.path.join(output_dir, "aligned_reads.bam")))
        command1.extend(["OUTPUT=",os.path.join(output_dir,"insert_metrics.txt")])
        command1.extend(["HISTOGRAM_FILE=" + os.path.join(output_dir, "insert_size_histogram.pdf")])
        self.run_cmd(command1)

        command2 = ["java"]
        command2.extend(["-jar", os.path.join(self.path, "picard.jar")])
        command2.extend(["CollectInsertSizeMetrics INPUT=", os.path.join(output_dir, "aligned_reads.bam")])
        command2.extend(["OUTPUT=", os.path.join(output_dir,"insert_metrics.txt" )])
        command2.extend(["HISTOGRAM_FILE=" + os.path.join(output_dir, "insert_size_histogram.pdf")])
        self.run_cmd(command2)

        command3 = ["samtools"]
        command3.append("depth")
        command3.extend(["-a", os.path.join(output_dir, "aligned_reads.bam")])
        command3.extend([">", os.path.join(output_dir, "depth_out.txt")])
        self.run_cmd(command3)

    def variant_calling(self, assembly_file, output_dir):
        command = ["java"]
        command.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("HaplotypeCaller")
        command.extend(["-R", assembly_file])
        command.extend(["-I", os.path.join(output_dir, "aligned_reads.bam")])
        command.extend(["-O", os.path.join(output_dir, "raw_variants.vcf")])
        self.run_cmd(command)

    def extract_variants(self, assembly_file, output_dir):
        command1 = ["java"]
        command1.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command1.append("SelectVariants")
        command1.extend(["-R", assembly_file])
        command1.extend(["-V", os.path.join(output_dir, "raw_variants.vcf")])
        command1.extend(["--select-type", "SNP"])
        command1.extend(["-O", os.path.join(output_dir, "raw_snps.vcf")])
        self.run_cmd(command1)

        command2 = ["java"]
        command2.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command2.append("SelectVariants")
        command2.extend(["-R", assembly_file])
        command2.extend(["-V", os.path.join(output_dir, "raw_variants.vcf")])
        command2.extend(["--select-type", "INDEL"])
        command2.extend(["-O", os.path.join(output_dir, "raw_indels.vcf")])
        self.run_cmd(command2)


    def filter_SNPs(self, assembly_file, output_file, output_dir, params):
        command  = ["java"]
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
        command  = ["java"]
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
        command1 = ["java"]
        command1.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command1.extend(["SelectVariants", "--exclude-filtered"])
        command1.extend(["-V", os.path.join(output_dir, "filtered_snps.vcf")])
        command1.extend(["-O", os.path.join(output_dir, "bqsr_snps.vcf")])
        self.run_cmd(command1)

        command2 = ["java"]
        command2.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command2.extend(["SelectVariants", "--exclude-filtered"])
        command2.extend(["-V", os.path.join(output_dir, "filtered_indels.vcf")])
        command2.extend(["-O", os.path.join(output_dir, "bqsr_indels.vcf")])
        self.run_cmd(command2)

    def base_quality_score_recalibration(self, assembly_file, data_table, output_dir):
        command = ["java"]
        command.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("BaseRecalibrator")
        command.extend(["-R", assembly_file])
        command.extend(["-I", os.path.join(output_dir, "aligned_reads.bam")])
        command.extend(["-known-sites", os.path.join(output_dir, "bqsr_snps.vcf")])
        command.extend(["-known-sites", os.path.join(output_dir, "bqsr_indels.vcf")])
        command.extend(["-O", os.path.join(output_dir, data_table)])
        self.run_cmd(command)
    
    def apply_BQSR(self, assembly_file, data_table, output_dir):
        command = ["java"]
        command.extend(["-jar", os.path.join(self.path, "gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("ApplyBQSR")
        command.extend(["-R", assembly_file])
        command.extend(["-I", os.path.join(output_dir,"aligned_reads.bam")])
        command.extend(["-bqsr",os.path.join(output_dir, data_table)])
        command.extend(["-O", os.path.join(output_dir, "recal_reads.bam")])
        self.run_cmd(command)

    def analyze_covariates(self, output_dir):
        command = ["java"]
        command.extend(["-jar", os.path.join(self.path,"gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar")])
        command.append("AnalyzeCovariates")
        command.extend(["-before", os.path.join(output_dir,"recal_data.table")])
        command.extend(["-after", os.path.join(output_dir, "post_recal_data.table")])
        command.extend(["-plots", os.path.join(output_dir, "recalibration_plots.pdf")])
        self.run_cmd(command)
