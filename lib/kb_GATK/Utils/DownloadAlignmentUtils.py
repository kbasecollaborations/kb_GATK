import os
from installed_clients.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from installed_clients.AssemblyUtilClient import AssemblyUtil

class DownloadAlignmentUtils:

    def __init__(self):
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.au = AssemblyUtil(self.callbackURL)
        self.ru = ReadsAlignmentUtils(self.callbackURL)
        pass

    def download_genome(self, assembly_ref, output_dir):
        file = self.au.get_assembly_as_fasta({
            'ref': assembly_ref,
            'filename': os.path.join(output_dir, "ref_genome.fna")
        })
        return file

    def downloadreadalignment(self, params):
        params['downloadSAM'] = 1
        return self.ru.download_alignment(params)
  


