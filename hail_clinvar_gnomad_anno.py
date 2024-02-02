import hail as hl
import xml.etree.ElementTree as ET
import pandas as pd
import csv
import json
from collections import Counter
from Bio import Entrez
import re
import os
import logging
import argparse

Entrez.api_key = '8ffb436a03ddf2e6b41406517e8562b57e08'
Entrez.email = 'phonegor2@gmail.com'

temp_dir = '/mnt/beegfs/tmp'
hl.init(
  master='local[56]',
  tmp_dir = temp_dir,
  local_tmpdir = temp_dir,
  spark_conf={"spark.local.dir": temp_dir}
)

# Load your Hail VariantDataset first
# hl.import_vcf('/mnt/beegfs/hongyf/SLURM/results/variant_calling/haplotypecaller/joint_variant_calling/joint_germline_recalibrated.vcf.gz', reference_genome = 'GRCh38', force_bgz = True, min_partitions = 28).write('/mnt/beegfs/hongyf/SLURM/results/variant_calling/haplotypecaller/joint_variant_calling/joint_germline_recalibrated.mt', overwrite = True)
# hl.import_vcf('/mnt/beegfs/hongyf/data/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.noINFO.vcf.gz', reference_genome = 'GRCh38', force_bgz = True, min_partitions = 28).write('/mnt/beegfs/hongyf/data/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.noINFO.mt', overwrite = True)
# hl.import_vcf('/mnt/beegfs/hongyf/data/gnomad.genomes.v3.1.2.hgdp_tgp.chr21.vcf.bgz', reference_genome = 'GRCh38', force_bgz = True, min_partitions = 28).write('/mnt/beegfs/hongyf/data/gnomad.genomes.v3.1.2.hgdp_tgp.chr21.mt', overwrite = True)

# Define a function that translates GT calls to simplified zygosity strings using Hail expressions
def get_zygosity(call):
    return hl.if_else(
        call.is_het(),
        'het',  # Heterozygous
        hl.if_else(
            call.is_hom_ref() | call.is_hom_var(),
            'hom',  # Homozygous (either reference or variant)
            'unknown'  # For no-call or multi-allelic cases
        )
    )

def convert_review_status(row):
    #condition_data = json.loads(row)
    #for rcv_key, rcv_dict in condition_data.items():
    #if 'review_status' in row:
    if row == 'criteria provided, multiple submitters, no conflicts':
        return '\u2605\u2605'
    elif row == 'reviewed by expert panel':
        return '\u2605\u2605\u2605'
    elif row == 'practice guideline':
        return '\u2605\u2605\u2605\u2605'
    else:
        return None
    
# Define the function to retrieve the information
def get_medgen_info(CUI_list):
    disease_list = []
    inheritance_list = []
    disease_detail_list = []
    disease_source_list = []
    
    for CUI in CUI_list:
        disease = None
        inheritance = None
        disease_detail = None
        disease_source = None
        
        # Perform the Entrez search
        handle = Entrez.esearch(db="medgen", retmax=1, term=CUI + "[Accession]", idtype="acc")
        records = Entrez.read(handle)
        handle.close()
        # Check if any IDs were found
        if records['IdList']:
            medgen_UID = records['IdList'][0]
            
            # Perform the Entrez summary
            handle = Entrez.esummary(db="medgen", id=medgen_UID, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            document_summary = records['DocumentSummarySet']['DocumentSummary'][0]
            # Retrieve 'Disease'
            if document_summary.get('Title'):
                disease = document_summary['Title']
                
            concept_meta = document_summary['ConceptMeta']
            # Retrieve 'ModesOfInheritance'
            if concept_meta.get('ModesOfInheritance'):
                inheritance = concept_meta['ModesOfInheritance'][0]['Name'].rstrip('inheritance').rstrip()
                
            # Retrieve 'Definitions'
            if concept_meta.get('Definitions'):
                disease_detail = re.sub(r'(\\n)+', '\n', str(concept_meta['Definitions'][0]))
                disease_source = concept_meta['Definitions'][0].attributes['source']
            
        disease_list.append(disease)
        inheritance_list.append(inheritance)
        disease_detail_list.append(disease_detail)
        disease_source_list.append(disease_source)
    
    return disease_list, inheritance_list, disease_detail_list, disease_source_list

def multiallelelics_split(joint_mt):
    # Filter to biallelic variants and annotate
    bi = joint_mt.filter_rows(hl.len(joint_mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False)
    
    # Filter to multiallelic variants and split
    multi = joint_mt.filter_rows(hl.len(joint_mt.alleles) > 2)
    split = hl.split_multi_hts(multi)
    
    # Combine the biallelic and split multiallelic variants
    return split.union_rows(bi)

def process_joint_matrix_table(joint_mt, clinvar_ht):
    """
    Processes the joint matrix table by splitting multiallelic variants,
    annotating with ClinVar data, filtering, and updating filters.

    :param joint_mt: The Hail MatrixTable of the joint dataset.
    :param clinvar_ht: The Hail Table containing ClinVar annotations.
    :return: The processed Hail MatrixTable.
    """
    # Split multiallelic variants if necessary
    if joint_mt.aggregate_rows(hl.agg.any(hl.len(joint_mt.alleles) > 2)):
        joint_mt = multiallelelics_split(joint_mt)

    # Annotate rows with ClinVar data
    joint_mt = joint_mt.annotate_rows(**clinvar_ht[joint_mt.row_key])

    # Filter rows to those defined in the ClinVar dataset
    joint_mt = joint_mt.filter_rows(hl.is_defined(joint_mt.Name))

    # Update the filters field: if empty, assign {"PASS"}
    joint_mt = joint_mt.annotate_rows(filters=hl.if_else(
        hl.len(joint_mt.filters) == 0, hl.set(["PASS"]), joint_mt.filters))

    joint_mt = joint_mt.drop(joint_mt.info)
    return joint_mt

def process_sample_matrix_table(sample_mt, gnomad_mt, pop_code):
    """
    Annotates the sample matrix table with zygosity, genotype, MedGen, and population-specific allele frequency.

    :param sample_mt: The Hail MatrixTable for a specific sample.
    :param gnomad_mt: The Hail MatrixTable containing gnomAD data.
    :param pop_code: Population code to use for allele frequency ('global' or 'eas').
    :return: A pandas DataFrame with the necessary annotations and allele frequencies.
    """
    # Annotate entries with zygosity and genotype
    sample_mt = sample_mt.annotate_entries(
        Zygosity=get_zygosity(sample_mt.GT),
        Genotype=hl.str(sample_mt.alleles[sample_mt.GT[0]]) + "/" + hl.str(sample_mt.alleles[sample_mt.GT[1]])
    )

    # Annotate rows with MedGen information
    sample_mt = sample_mt.annotate_rows(
        MedGen=sample_mt.ReportedPhenotypeInfo.map(lambda pheno: pheno.split(':')[0])
    )

    # Annotate rows with the appropriate allele frequency based on the population code
    if pop_code == 'global':
        sample_mt = sample_mt.annotate_rows(
            AF_global=gnomad_mt.rows()[sample_mt.row_key].info.AF
        )
        sample_pd = sample_mt.make_table().to_pandas()
        sample_pd['Allele frequency (global)'] = sample_pd['AF_global'].apply(lambda x: f"{x[0]*100:.5f}%")
    elif pop_code == 'eas':
        sample_mt = sample_mt.annotate_rows(
            AF_eas=gnomad_mt.rows()[sample_mt.row_key].info.AF_eas
        )
        sample_pd = sample_mt.make_table().to_pandas()
        sample_pd['Allele frequency (EAS)'] = sample_pd['AF_eas'].apply(lambda x: f"{x[0]*100:.5f}%")

    return sample_pd

def process_sample_pd(sample_pd, sample):
    """
    Processes a pandas DataFrame with various annotations.

    :param sample_pd: The pandas DataFrame containing sample data.
    :param sample: The name of the sample, used for stripping from column names.
    :return: The processed pandas DataFrame.
    """
    # Annotate the DataFrame with disease information
    disease_info_cols = ['Disease associated', 'Inheritance', 'Disease detail', 'Disease source']
    sample_pd[disease_info_cols] = sample_pd.apply(
        lambda row: pd.Series(get_medgen_info(row['MedGen'])),
        axis=1
    )

    # Convert review status using the provided function
    sample_pd['Review status'] = sample_pd['ReviewStatus'].apply(convert_review_status)

    # Clean up the variant names
    sample_pd['Variant'] = sample_pd['Name'].apply(lambda x: re.sub(r'\([^)]+\)', '', x))
    sample_pd['Gene'] = sample_pd['Name'].apply(lambda x: re.search(r'\((.*?)\)', x).group(1))

    # Rename columns
    rename_columns = {'rsid': 'rsID', 'OriginSimple': 'Origin', 'ClinicalSignificance': 'Clinical significance'}
    sample_pd.rename(columns=rename_columns, inplace=True)

    # Strip sample prefix from column names
    pattern_lstrip = r'^' + re.escape(sample + '.')
    sample_pd.columns = sample_pd.columns.str.replace(pattern_lstrip, '', regex=True)

    sample_pd.drop(columns = ['RS# (dbSNP)', 'RCVaccession'], inplace=True)
    return sample_pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the VCF matrix table to pandas df and annotate it with ClinVar and gnomAD data.')
    parser.add_argument('--clinvar_ht_path', default='/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/ClinVar/clinvar_summary_20240127.ht', help='Path to the ClinVar Hail table.')
    parser.add_argument('--gnomad_mt_path', default='/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/eas_af-only-gnomad.genomes.v4.0.sites.mt', help='Path to the gnomAD Hail matrix table.')
    parser.add_argument('--vcf_mt_path', default='/mnt/beegfs/hongyf/SLURM/results/variant_calling/haplotypecaller/joint_variant_calling/joint_germline_recalibrated.mt', help='Path to the VCF Hail matrix table.')
    parser.add_argument('--output_path', default='/mnt/beegfs/hongyf/SLURM/results/annotation/Hail/', help='Output directory for annotated results.')
    parser.add_argument('--pop_code', choices=['global', 'eas'], default='eas', help='Population code to use for allele frequency.')
    parser.add_argument('--samples', help='Comma-separated list of sample IDs to process. If not provided, all samples will be processed.')
    
    # Parse the arguments
    args = parser.parse_args()

    clinvar_ht_path = args.clinvar_ht_path
    gnomad_mt_path = args.gnomad_mt_path
    vcf_mt_path = args.vcf_mt_path

    # Configure logging
    logging.basicConfig(filename = args.output_path + '.command.log',
                    level = logging.INFO,
                    format = '%(asctime)s %(levelname)s:%(message)s')
    
    clinvar_ht = hl.read_table(clinvar_ht_path)
    gnomad_mt = hl.read_matrix_table(gnomad_mt_path)
    joint_mt = hl.read_matrix_table(vcf_mt_path)
    joint_mt = process_joint_matrix_table(joint_mt, clinvar_ht)

    samples = args.samples.split(',') if args.samples else joint_mt.s.collect()

    for sample in samples:
        logging.info(f'Processing sample: {sample}')
        
        sample_dir = os.path.join(args.output_path, sample)
        output = os.path.join(sample_dir, 'var_details.tsv')
        # Create the directory if it doesn't exist
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        else:
            if os.path.exists(output) and  os.stat(output).st_size == 0:
                logging.info(f'No valid pathogenic variants found for sample: {sample}')
                continue


        sample_mt = joint_mt.filter_cols(joint_mt.s == sample)
        # Filter out entries where the genotype is homozygous reference or missing
        sample_mt = sample_mt.filter_rows(hl.agg.any(~sample_mt.GT.is_hom_ref() & hl.is_defined(sample_mt.GT)))
        
        valid_variant_count = sample_mt.count_rows()
        if valid_variant_count > 0:
            logging.info(f'Found {valid_variant_count} valid pathogenic variants for sample: {sample}')
            sample_pd = process_sample_matrix_table(sample_mt, gnomad_mt, args.pop_code)
            sample_pd = process_sample_pd(sample_pd, sample)
            sample_pd.to_csv(output, sep='\t', index=False, quoting = csv.QUOTE_MINIMAL)
        else:
            logging.info(f'No valid pathogenic variants found for sample: {sample}')
            # Create an empty TSV file
            with open(output, 'w') as file:
                pass

hl.stop()

#sample_mt = sample_mt.select_rows(sample_mt.rsid, sample_mt.filters, 
                                  #sample_mt.VariationID, sample_mt.Name, sample_mt.Type, sample_mt.GeneSymbol, sample_mt.OriginSimple, sample_mt.ClinicalSignificance_y, sample_mt.ReviewStatus_y,
                                  #sample_mt.ReportedPhenotypeInfo, sample_mt.Submitter, sample_mt.Description, sample_mt.ExplanationOfInterpretation, sample_mt.AF_eas)

#vcf_mt = vcf_mt.select_rows(vcf_mt.rsid, vcf_mt.filters, vcf_mt.VariationID, vcf_mt.Name, vcf_mt.Type, vcf_mt.GeneSymbol, vcf_mt.OriginSimple, vcf_mt.ClinicalSignificance_y, vcf_mt.ReviewStatus_y, 
                            #vcf_mt.ReportedPhenotypeInfo, vcf_mt.Submitter, vcf_mt.Description, vcf_mt.ExplanationOfInterpretation, vcf_mt.AF_global)

#sample_pd['condition'] = sample_pd.apply(lambda row: supportive_explanation_most_common(row, 'origin'), axis=1)
#sample_pd['condition'] = sample_pd.apply(lambda row: supportive_explanation_most_common(row, 'inheritance'), axis=1)
#sample_pd['condition'] = sample_pd['condition'].apply(filter_condition)
#sample_pd['condition'] = sample_pd['condition'].apply(convert_review_status)
#sample_pd['condition'] = sample_pd['condition'].apply(lambda row: rename_rcv_header(row, 'disease', 'Disease'))
#sample_pd['condition'] = sample_pd['condition'].apply(lambda row: rename_rcv_header(row, 'origin', 'Origin'))
#sample_pd['condition'] = sample_pd['condition'].apply(lambda row: rename_rcv_header(row, 'inheritance', 'Inheritance'))
#sample_pd['Gene'] = sample_pd['variation_name'].apply(lambda x: re.search(r'\((.*?)\)', x).group(1))
#sample_pd['Strand'] = sample_pd['info.STRAND'].apply(modify_strand)
#sample_pd['Zygosity'] =  sample_pd[sample_name +'.GT'].apply(determine_zygosity)
#sample_pd['Review status'] = sample_pd['ReviewStatus_y'].apply(convert_review_status)
#sample_pd['Allele frequency (EAS)'] = sample_pd['AF_eas'].apply(lambda x: f"{x[0]*100:.5f}%")
#sample_pd['Clinical significance'] = sample_pd['ClinicalSignificance_y']
#sample_pd['Disease associated'] = sample_pd['ReportedPhenotypeInfo'].apply(lambda pheno_list: [pheno.split(':')[1] for pheno in pheno_list])
#sample_pd['MedGen'] = sample_pd['ReportedPhenotypeInfo'].apply(lambda pheno_list: [pheno.split(':')[0] for pheno in pheno_list])
#sample_pd[['Inheritance', 'Disease detail', 'Disease detail source']] = sample_pd.apply(lambda row: pd.Series(get_medgen_info(row['MedGen'])), axis=1)
#sample_pd['Genotype'] = sample_pd.apply(lambda row: row['alleles'][row[sample_name +'.GT'].alleles[0]] + '/' + row['alleles'][row[sample_name +'.GT'].alleles[1]], axis=1)

#sample_pd['Variant'] = sample_pd['Name'].apply(lambda x: re.sub(r'\([^)]+\)', '', x))



#sample_pd[['rsID', 'Genotype', 'Zygosity', 'Gene', 'Variant', 'Clinical significance', 'Origin', 'Allele frequency (EAS)', 
            #'Disease associated', 'Review status', 'Inheritance', 'Disease detail source', 'Disease detail', 
            #'Submitter', 'Description', 'ExplanationOfInterpretation']].to_csv('/home/methylation/report-generator/test1_var_details_0128.tsv', sep = '\t', index = False)

# NA12878_chr1_pd[['rsID', 'Genotype', 'Zygosity', 'Gene', 'Variant', 'Clinical significance', 'Allele frequency (global)', 'Review status', 'Origin', 'ReportedPhenotypeInfo', 'Submitter', 'Description', 'ExplanationOfInterpretation']].to_csv('/home/methylation/report-generator/NA12878_var_details_0128.tsv', sep='\t',index=False,  quoting=csv.QUOTE_NONE)

""" def most_common_in_json(explanation_dict, key):
    values = [item[key] for item in explanation_dict.values() if key in item]
    # Find the most common non-None value
    for value, count in Counter(values).most_common():
        if value is not None:
            return value
    return "NO info in ClinVar, please check manually" """

""" def supportive_explanation_most_common(row, col_key):
    condition_data = json.loads(row['condition'])
    explanation_data = json.loads(row['explanation'])
    for rcv_key in condition_data.keys():
        keys_to_extract = condition_data[rcv_key].get('supportive_clinicalassertionID', '').split(',')
        # If keys_to_extract is empty, get every key from explanation_data
        if keys_to_extract == ['']:
            keys_to_extract = list(explanation_data.keys()) 
        # Extract values based on the keys
        most_common_value_in_col = most_common_in_json({k: explanation_data.get(k) for k in keys_to_extract}, col_key)
        # Update condition_data with the extracted value
        condition_data[rcv_key][col_key] = most_common_value_in_col
    # Convert the modified condition_data back to a string
    updated_condition_string = json.dumps(condition_data)
    # Return the updated string
    return updated_condition_string """

""" def modify_strand(value):
    if value == 1:
        return '+'
    elif value == -1:
        return '-'
    else:
        return value  # Return the original value if it's not 1 or -1 """

""" def rename_rcv_header(row, old, new):
    condition_data = json.loads(row)
    for rcv_key, rcv_dict in condition_data.items():
        if old in rcv_dict:
            rcv_dict[new] = rcv_dict.pop(old)
    return json.dumps(condition_data) """

""" def filter_condition(json_str):
    # Parse the JSON string
    condition_dict = json.loads(json_str)
    # Access the 'rcv_index' list
    rcv_index_list = list(condition_dict.keys())
    # Check conditions for 'practice guideline' and 'reviewed by expert panel'
    has_practice_guideline = any(condition_dict[rcv]['review_status'] == 'practice guideline' for rcv in rcv_index_list)
    if has_practice_guideline:
        filtered_rcv_index = [rcv for rcv in rcv_index_list if condition_dict[rcv]['review_status'] == 'practice guideline']
        filtered_condition = {key: condition_dict[key] for key in filtered_rcv_index if key in condition_dict}
        return json.dumps(filtered_condition)
    has_expert_panel = any(condition_dict[rcv]['review_status'] == 'reviewed by expert panel' for rcv in rcv_index_list)
    if has_expert_panel:
        filtered_rcv_index = [rcv for rcv in rcv_index_list if condition_dict[rcv]['review_status'] == 'reviewed by expert panel']
        filtered_condition = {key: condition_dict[key] for key in filtered_rcv_index if key in condition_dict}
        return json.dumps(filtered_condition)
    
    return json_str """

