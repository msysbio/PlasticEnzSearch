import unittest
from abundance import annotation
from args import Args
import translate_search
from quantify_hmm import quantify_hmm
from annotation import blast_search
from pathmanager import PathManager
from output import create_html, remove_temps
import PlasticTools
import os
import database_operations

data_folder = "/home/jasper/Thesis/"

args1 = Args(
        output = f'{data_folder}Test_data/output',
        contigs = f'{data_folder}Test_data/contigs-fixed.fa',
        plastic = 'all',
        mappings=f'{data_folder}Test_data/GC125633.bam,{data_folder}Test_data/GC125648.bam,{data_folder}Test_data/GC125657.bam,{data_folder}Test_data/GC125668.bam')
       
args2 = Args(     
        output = f'{data_folder}Test_data2/output',
        contigs = f'{data_folder}Test_data2/Citadell.contigs.fa',
        plastic = 'all',
        mappings=f'{data_folder}Test_data2/GC125618.sorted.bam,{data_folder}Test_data2/GC127864.sorted.bam')


class TestPlasticTools(unittest.TestCase):
    def setUp(self):
        self.args = args2
        self.p = PathManager(self.args)
        
        #database_operations.database_fetch(self.args.plastic, self.args.output)
    
    def test1_run_prodigal(self):
        # Try to run the function
        try:
            translate_search.run_prodigal(self.p)
        except Exception as e:
            print(e)
    
        # Check that the output directories were created
        self.assertTrue(os.path.exists(self.p.temps))
    
        # Check that the .faa and .ffn files were created
        contigs_base = os.path.basename(self.p.contigs).split(".")[0]
        aa_file = os.path.join(self.p.output, "temps", f"{contigs_base}.faa")
        nt_file = os.path.join(self.p.output, "temps", f"{contigs_base}.ffn")
        self.assertTrue(os.path.exists(aa_file))
        self.assertTrue(os.path.exists(nt_file))
    
        # Check that the prodigal log file was created
        prodigal_log_file = os.path.join(self.p.output, "temps", "prodigal.log")
        self.assertTrue(os.path.exists(prodigal_log_file))
    
    def test2_run_hmmer(self):
        # Try to run the function
        try:
            translate_search.run_hmmer(self.p)
        except Exception as e:
            print(e)
        
        # Define contigs_base
        contigs_base = os.path.basename(self.p.contigs).split(".")[0]
    
        # Check that the hmmsearch output and log files were created for each plastic type
        if self.p.plastic == "all":
            plastic_names = self.p.all_plastics
        else:
            plastic_names = self.p.plastic.split(',')
        for plastic_name in plastic_names:
            temp_dir = os.path.join(self.p.temps, plastic_name.lower())
            hmm_output = os.path.join(temp_dir, f"{contigs_base}_{plastic_name}_HMMER.out")
            log_file = os.path.join(temp_dir, f"{plastic_name}_hmmsearch.log")
            self.assertTrue(os.path.exists(hmm_output), f"{plastic_name}: {hmm_output}")
            self.assertTrue(os.path.exists(log_file), f"{plastic_name}: {log_file}")
    
        # Check that the hmmsearch output fasta file was created for each plastic type
        for plastic_name in plastic_names:
            temp_dir = os.path.join(self.p.temps, plastic_name.lower())
            output = os.path.join(temp_dir, f"{contigs_base}_{plastic_name}_hmm_output.fasta")
            self.assertTrue(os.path.exists(output), f"{plastic_name}: {output}")


    def test4_blast_search(self):
        # Try to run the function
        try:
            blast_search.blast_search(self.p)
        except Exception as e:
            print(e)
    
        # Check that the output directories were created
        self.assertTrue(os.path.exists(self.p.temps))
    
        # Check that the annotation log file was created
        log_file = os.path.join(self.p.temps, "annotation.log")
        self.assertTrue(os.path.exists(log_file), log_file)
    
        # Check that the annotation Excel files were created for each plastic type
        if self.p.plastic == "all":
            plastic_names = self.p.all_plastics
        else:
            plastic_names = self.p.plastic.split(',')
        for plastic_name in plastic_names:
            temp_dir = os.path.join(self.p.temps, plastic_name.lower())
            output_file = os.path.join(temp_dir, f"{self.p.contigs_base}_{plastic_name}_hmm_output_annotation.xlsx")
            self.assertTrue(os.path.exists(output_file), f"{plastic_name}: {output_file}")

    def test5_annotation(self):
        # Try to run the function
        try:
            annotation(self.p)
        except Exception as e:
            print(e)
    
        # Check that the output directories were created
        self.assertTrue(os.path.exists(self.p.temps))
        
        # Define contigs_base
        contigs_base = os.path.basename(self.p.contigs).split(".")[0]
    
        # Check that the corrected .ffn and .saf files were created for each plastic type
        if self.p.plastic == "all":
            plastic_names = self.p.all_plastics
        else:
            plastic_names = self.p.plastic.split(',')
        for plastic_name in plastic_names:
            temp_dir = os.path.join(self.p.temps, plastic_name.lower())
            corrected_ffn_file = os.path.join(temp_dir, f"{contigs_base}_{plastic_name}_corrected.ffn")
            saf_file = os.path.join(temp_dir, f"{contigs_base}_{plastic_name}.saf")
            self.assertTrue(os.path.exists(corrected_ffn_file))
            self.assertTrue(os.path.exists(saf_file))
    
        # Check that the featureCounts log and output files were created for each plastic type
        if self.p.plastic == "all":
            plastic_names = self.p.all_plastics
        else:
            plastic_names = self.p.plastic.split(',')
        for plastic_name in plastic_names:
            temp_dir = os.path.join(self.p.temps, plastic_name.lower())
            mapping_files = self.p.mappings.split(',')
            for mapping_file in mapping_files:
                mapping_file = mapping_file.strip()  # Remove any leading/trailing whitespace
                log_file = os.path.join(temp_dir, f"{os.path.basename(mapping_file)}_featureCounts.log")
                program_output_file = os.path.join(temp_dir, f"{os.path.basename(mapping_file)}_featureCounts.out")
                self.assertTrue(os.path.exists(log_file), log_file)
                self.assertTrue(os.path.exists(program_output_file), program_output_file)
    
        # Check that the mapping_summary.tsv file was created for each plastic type
        if self.p.plastic == "all":
            plastic_names = self.p.all_plastics
        else:
            plastic_names = self.p.plastic.split(',')
        for plastic_name in plastic_names:
            temp_dir = os.path.join(self.p.temps, plastic_name.lower())
            tsv_file = os.path.join(temp_dir, 'mapping_summary.tsv')
            self.assertTrue(os.path.exists(tsv_file), tsv_file)
        
    def test6_html(self):
        # Try to run the function
        try:
            create_html(self.p)
            remove_temps(self.p)
        except Exception as e:
            print(e)

        html_files = [file for file in os.listdir(self.p.output) if file.endswith('.html')]
        self.assertTrue(html_files != [], html_files)

def run_tests():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestPlasticTools)
    suite._tests.sort(key=lambda x: x._testMethodName)
    unittest.TextTestRunner().run(suite)

if __name__ == '__main__':
    run_tests()
    #PlasticTools.main(args1, debug=True, blast=False)
    #PlasticTools.main(args2, debug=True, blast=False)

