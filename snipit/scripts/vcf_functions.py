#!/usr/bin/env python3

# imports of built-ins
import os
import sys
import argparse
import collections
from itertools import cycle, chain
import csv
import math
import re

# imports from other modules
from Bio import SeqIO
import vcf

aa_letter_to_abbrv = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "*": "Ter"
}

aa_abbrv_to_letter = {v:k for k,v in aa_letter_to_abbrv.items()}

def get_nonsynonymous_sub(alt_base, ann_list):
    strs_to_replace_re = re.compile(r'|'.join(re.escape(key) for key in aa_abbrv_to_letter.keys()))

    for ann in ann_list:
        base,var_type,effect_size,gene,*var_properties = ann.split("|")
        if var_type == "synonymous_variant":
            continue
        if base.upper() not in ["A","C","T","G"]:
            continue
        if base==alt_base:
            aa_sub = var_properties[6]
            if aa_sub.startswith("p."):
                aa_sub = aa_sub.replace("p.","")
                aa_sub = strs_to_replace_re.sub(lambda x: aa_abbrv_to_letter.get(x.group(), "_"), aa_sub)
                return aa_sub
            else:
                Exception("Does not look like an aa sub")

def parse_vcf(vcf_file):
    """
        Parse snpEff annotations from a vcf file
    """
    substitution_annot = {}
    with open(vcf_file, 'r') as infile:
        vcf_reader = vcf.Reader(infile)
        for record in vcf_reader:
            record_alt_annots = {}
            print(record.INFO["ANN"])
            for alt in record.ALT:
                print(alt)
                record_alt_annots[str(alt)] = get_nonsynonymous_sub(str(alt), record.INFO["ANN"])
            substitution_annot[record.POS] = record_alt_annots
    return substitution_annot

def get_vcf_sub_for_position_base(pos,base):
    return None

if __name__ == "__main__":
    # call parse_vcf(), passing the first argument. Otherwise print help text
    if len(sys.argv) < 2:
        print("Usage: python vcf_functions.py <vcf_file>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    if not os.path.isfile(vcf_file):
        print(f"Error: {vcf_file} is not a valid file.")
        sys.exit(1)
    if not vcf_file.endswith(".vcf"):
        print(f"Error: {vcf_file} does not have a .vcf extension.")
        sys.exit(1)
    # Parse the VCF file and print the substitutions
    print(f"Parsing VCF file: {vcf_file}")
    vars = parse_vcf(vcf_file)
    print(vars)