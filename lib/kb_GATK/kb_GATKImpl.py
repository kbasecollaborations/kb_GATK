# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
from installed_clients.VariationUtilClient import VariationUtil
from installed_clients.KBaseReportClient import KBaseReport
from kb_GATK.Utils.GATKUtils import GATKUtils
from kb_GATK.Utils.DownloadAlignmentUtils import DownloadAlignmentUtils


import shutil
#END_HEADER


class kb_GATK:
    '''
    Module Name:
    kb_GATK

    Module Description:
    A KBase module: kb_GATK
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbasecollaborations/kb_GATK.git"
    GIT_COMMIT_HASH = "5e6e4bdca9a7749bba0abab081736c56007212ed"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.gu = GATKUtils()
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.vu = VariationUtil(self.callback_url)
        self.du = DownloadAlignmentUtils(self.callback_url)
        #END_CONSTRUCTOR
        pass


    def run_kb_GATK(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_GATK

        

        source_ref1 = '43745/111/1'
        out = self.du.downloadreadalignment(source_ref1, params, self.callback_url)
        #sam_file = os.path.join(out['destination_dir'],"reads_alignment.sam")

        sam_file = "/kb/module/test/reads_alignment.sam"
        src = "/kb/module/test/genome"
        output_dir = self.shared_folder
        #output_dir = os.path.join(self.shared_folder, str(uuid.uuid4()))
        #os.mkdir(output_dir)
        dest = shutil.copytree(src, os.path.join(output_dir, "genome"))

        assembly_file = os.path.join(src, "reference/test_assembly.fa")

        fwd_fastq = "/kb/module/test/bt_test_data/reads_1.fq"
        rev_fastq = "/kb/module/test/bt_test_data/reads_2.fq"

        output_dir = self.shared_folder + "/"  #TODO need to use uuid for storing all the intermediate results rather that using shared_folder.

        self.gu.build_genome(assembly_file)

        self.gu.index_assembly(assembly_file)

        self.gu.generate_sequence_dictionary(assembly_file)

        #self.gu.mapping_genome(assembly_file, fwd_fastq, rev_fastq, output_dir )

        self.gu.duplicate_marking(output_dir, sam_file)

        self.gu.sort_bam_index(output_dir)

        self.gu.collect_alignment_and_insert_size_metrics(assembly_file, output_dir)
   
        self.gu.analyze_covariates( output_dir)
   
        self.gu.variant_calling(assembly_file, output_dir)
   
        self.gu.extract_variants(assembly_file, output_dir)
   
        self.gu.filter_SNPs(assembly_file, "filtered_snps.vcf", output_dir, params)
   
        self.gu.filter_Indels(assembly_file, "filtered_indels.vcf", output_dir, params)
   
        self.gu.exclude_filtered_variants(output_dir)
   
        self.gu.base_quality_score_recalibration(assembly_file, "recal_data.table", output_dir)
   
        self.gu.apply_BQSR(assembly_file, "recal_data.table", output_dir)
   
        self.gu.base_quality_score_recalibration(assembly_file, "post_recal_data.table", output_dir)
   
        self.gu.apply_BQSR(assembly_file,  "post_recal_data.table", output_dir)
   
        self.gu.filter_SNPs(assembly_file, "filtered_snps_final.vcf", output_dir, params)
   
        self.gu.filter_Indels(assembly_file, "filtered_indels_final.vcf", output_dir, params)

        os.system("grep   '##fileformat' " + output_dir + "filtered_snps_final.vcf > " + output_dir + "sample.vcf")
        cmd = "grep -v  '##' " + output_dir + "filtered_snps_final.vcf >> " + output_dir + "sample.vcf"
        
        os.system(cmd)

        params['vcf_staging_file_path'] = output_dir + "sample.vcf"

        self.vu.save_variation_from_vcf(params)

        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': 'Success'},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_kb_GATK

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_GATK return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
