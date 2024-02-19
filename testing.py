#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 19:04:17 2023

@author: u0145079
"""

import unittest
from abundance import annotation

class TestAnnotation(unittest.TestCase):
    def test_annotation(self):
        # Set up a mock args object
        class Args:
            """
            output = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Plastic_tool/Plastic_tool'
            contigs = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/Daphnia_microbiome/contigs-new/st.donatus/contigs-fixed.fa'
            plastic = 'all'
            mappings='/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/Daphnia_microbiome/bams/daphnia_microbiome/st.donatus/GC127855-ready.bam,/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/Daphnia_microbiome/bams/daphnia_microbiome/st.donatus/GC127867-ready.bam'
            """

            output = '/home/jasper/OneDrive/Thesis/Plastic_tool'
            contigs = '/home/jasper/Test_data/contigs-fixed.fa'
            plastic = 'all'
            mappings='/home/jasper/Test_data/GC125633.bam,/home/jasper/Test_data/GC125648.bam,/home/jasper/Test_data/GC125657.bam,/home/jasper/Test_data/GC125668.bam'


        args = Args()

        # Run the function
        annotation(args)

        # Add assertions here to check the results
        # For example, check if certain files were created, or if the function
        # returned the expected result

if __name__ == '__main__':
    unittest.main()
