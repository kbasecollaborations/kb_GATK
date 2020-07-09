import os
from installed_clients.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from installed_clients.AssemblyUtilClient import AssemblyUtil

class DownloadAlignmentUtils:

    def __init__(self, callback_url):
       self.callbackURL = os.environ['SDK_CALLBACK_URL']
       self.au = AssemblyUtil(self.callbackURL)
       pass 

    def downloadreadalignment(self, source_ref, params, callback_url):
        self.callback_url = callback_url
        self.ru = ReadsAlignmentUtils(self.callback_url)
        params['source_ref'] = source_ref
        params['downloadSAM'] = 1

        params['destination_dir'] = '/kb/module/work/tmp'
        params['stats'] = {
                    "properly_paired":1,
                    "multiple_alignments":1,
                    "singletons":1,
                    "alignment_rate":1,
                    "unmapped_reads":1,
                    "mapped_reads":1,
                    "total_reads":1
                 }
        return self.ru.download_alignment(params)
  
    def download_genome(self, genomeref, output_dir):
        file = self.au.get_assembly_as_fasta({
          'ref': genomeref,
          'filename': os.path.join(output_dir, "ref_genome.fa")
        })
        return file

