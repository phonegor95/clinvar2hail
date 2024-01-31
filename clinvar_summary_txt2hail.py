import pandas as pd
import hail as hl
import gzip
import argparse
import re
import sys

temp_dir = '/mnt/beegfs/tmp'
hl.init(
  master='local[28]',
  tmp_dir = temp_dir,
  local_tmpdir = temp_dir,
  spark_conf={"spark.local.dir": temp_dir}
)

# Identify the header (last comment line)
def get_header_line(gzfile):
    with gzip.open(gzfile, 'rt') as file:
        header_line = None
        skip_rows = 0
        for line in file:
            if line.startswith('#'):
                header_line = line.strip()
                header_line = header_line.lstrip('#').split('\t')
            else:
                return header_line, skip_rows
            skip_rows += 1
        # Ensure that a header line is found
        if header_line is None:
            raise ValueError("No header line found in the file")
    
# Define a custom aggregation function
""" def custom_agg_review_status(series):
    if 'practice guideline' in series.values:
        return ['practice guideline']
    elif 'reviewed by expert panel' in series.values:
        return ['reviewed by expert panel']
    else:
        return list(series) """

def filter_review_status_to_priority(df):
    # Define the priority for each review status
    priority = {
        'practice guideline': 1,
        'reviewed by expert panel': 2
    }
    
    # Map the review statuses to priorities
    df['SortPriority'] = df['ReviewStatus'].map(lambda x: priority.get(x, 3))
    
    # Sort the DataFrame based on VariationID and SortPriority
    df.sort_values(by=['VariationID', 'SortPriority'], inplace=True)
    
    # Drop duplicates based on VariationID for rows with SortPriority < 3
    df_filtered = df[df['SortPriority'] < 3].drop_duplicates(subset='VariationID', keep='first')
    
    # Find VariationIDs that have been kept with SortPriority < 3
    kept_variation_ids = df_filtered['VariationID'].unique()
    
    # Filter out rows with SortPriority == 3 that have the same VariationID as those kept with SortPriority < 3
    df_no_priority = df[(df['SortPriority'] == 3) & (~df['VariationID'].isin(kept_variation_ids))]
    
    # Concatenate the filtered priority rows with the no_priority rows
    df = pd.concat([df_filtered, df_no_priority], ignore_index=True)
    
    # Drop the temporary SortPriority column
    df.drop('SortPriority', axis=1, inplace=True)
    
    # Return the modified DataFrame
    return df

def meets_review_criteria(review_status_list):
    # Count the occurrences of each review status
    review_status_counts = {
        'criteria provided, single submitter': review_status_list.count('criteria provided, single submitter'),
        'practice guideline': review_status_list.count('practice guideline'),
        'reviewed by expert panel': review_status_list.count('reviewed by expert panel')
    }
    
    # Check if the criteria are met
    has_single_submitter = review_status_counts['criteria provided, single submitter'] > 1
    has_one_practice_guideline = review_status_counts['practice guideline'] > 0
    has_one_reviewed_by_expert_panel = review_status_counts['reviewed by expert panel'] > 0
    
    # Return True if at least one of the criteria is met
    return has_single_submitter or has_one_practice_guideline or has_one_reviewed_by_expert_panel

def read_and_filter_variant_summary(variant_summary_path, variant_allele_relat_path, chromosome_criteria, review_status_criteria, clinical_significance_pathogenic, allele_type_filter, variant_type_filter):
    header_line, skip_rows = get_header_line(variant_summary_path)
    df = pd.read_csv(variant_summary_path, sep='\t', comment=None, skiprows=skip_rows, names=header_line, low_memory=False)
    df = df[df['Assembly'] == 'GRCh38']
    df = df[df['Chromosome'].isin(chromosome_criteria)]
    df.replace({'Chromosome': {'MT': 'M'}}, inplace=True)
    df['Chromosome'] = df['Chromosome'].apply(lambda x: 'chr' + str(x))
    df = df[(df['PositionVCF'] != -1) & (df['ReferenceAlleleVCF'] != 'na') &  (df['AlternateAlleleVCF'] != 'na')]
    df = df[~df['Type'].isin(allele_type_filter)]

    header_line, skip_rows = get_header_line(variant_allele_relat_path)
    variant_allele_relat_df = pd.read_csv(variant_allele_relat_path, sep='\t', comment=None, skiprows=skip_rows, names=header_line, low_memory=False)
    df = df[df['VariationID'].isin(variant_allele_relat_df[~variant_allele_relat_df['Type'].isin(variant_type_filter)]['VariationID'])]
    print(f"df after variant_type_filter {df}")
    df = df[df['ReviewStatus'].isin(review_status_criteria)]
    df = df[df['ClinicalSignificance'].str.contains('|'.join(clinical_significance_pathogenic), case=True, na=False)]

    if df.groupby('VariationID').agg({'AlleleID': set}).apply(lambda x: len(x) > 1).duplicated().any():
        print(f"Error: groupby VariationID Column: AlleleID has duplicate values.")
        sys.exit(1)  # Exit the script with an error code
    
    # until now, no duplicated alleles (locus+alleles) of every variation (locus) after above condition
    df = df.drop_duplicates(subset='VariationID', keep='first')
    return df

def read_and_filter_submission_summary(file_path, phenotype_filter):
    header_line, skip_rows = get_header_line(file_path)
    df = pd.read_csv(file_path, sep='\t', comment=None, skiprows=skip_rows, names=header_line, low_memory=False)
    df = df[~(df['ReportedPhenotypeInfo'].isin(phenotype_filter) | df['ReportedPhenotypeInfo'].str.startswith('na:'))]
    df = filter_review_status_to_priority(df)
    return df

def group_and_filter_submission_summary(df, clinical_significance_pathogenic, clinical_significance_benign):
    df = df.groupby(['VariationID', 'ReportedPhenotypeInfo']).agg({
        'ClinicalSignificance': list,
        'Description': list,
        'ReviewStatus': list,
        'Submitter': list,
        'SCV': list,
        'ExplanationOfInterpretation': list
    }).reset_index()

    # Subset the DataFrame based on the criteria
    df = df[df['ClinicalSignificance'].apply(
        lambda x: any(value in clinical_significance_pathogenic for value in x) and not any(value in clinical_significance_benign for value in x)
    )].reset_index(drop=True)

    df = df[df['ReviewStatus'].apply(meets_review_criteria)].reset_index(drop=True)

    # Group by 'VariationID' to collect all different 'ReportedPhenotypeInfo' entries
    df = df.groupby('VariationID').agg({
        'ReportedPhenotypeInfo': list,
        'ClinicalSignificance': list,
        'Description': list,
        'ReviewStatus': list,
        'Submitter': list,
        'SCV': list,
        'ExplanationOfInterpretation': list
    }).reset_index()
    return df

def merge_and_prepare_for_hail(submission_df, variant_df, drop_col):
    merged_df = submission_df.merge(variant_df, on='VariationID', suffixes=('_list', ''))
    merged_df = merged_df.drop(columns=drop_col)
    print(f"merged_df: {merged_df}")
    ht = hl.Table.from_pandas(merged_df)
    ht = ht.transmute(locus=hl.locus(ht.Chromosome, ht.PositionVCF, reference_genome='GRCh38'))
    ht = ht.transmute(alleles=[ht.ReferenceAlleleVCF, ht.AlternateAlleleVCF])
    ht = ht.key_by(ht.locus, ht.alleles)
    return ht

""" def annotate_with_gnomad(ht, gnomad_path):
    gnomad_ht = hl.read_matrix_table(gnomad_path)
    ht = ht.annotate(AF_global=gnomad_ht.rows()[ht.key].info.AF)
    ht = ht.annotate(AF_eas=gnomad_ht.rows()[ht.key].info.AF_eas)
    return ht """

def main():

    variant_df = read_and_filter_variant_summary(
        variant_summary,
        variant_allele_relat,
        Chromosome_criteria,
        ReviewStatus_criteria,
        ClinicalSignificance_Pathogenic,
        AlleleType_filter,
        VariantType_filter
    )
    print(f"variant_summary_df: {variant_df}")
    variant_df['VariationID'].duplicated().any()

    submission_df = read_and_filter_submission_summary(submission_summary, Phenotype_filter)
    submission_df = group_and_filter_submission_summary(
        submission_df, 
        ClinicalSignificance_Pathogenic, 
        ClinicalSignificance_Benign
    )
    print(f"submission_summary_df: {submission_df}")

    variant_submission_ht = merge_and_prepare_for_hail(submission_df, variant_df, drop_col)
    del submission_df 
    del variant_df
    # variant_submission_ht = annotate_with_gnomad(variant_submission_ht, gnomad)
    variant_submission_ht.describe()
    variant_submission_ht.show()

    variant_submission_ht.write(output_path, overwrite=True)

# Run the main function
if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description='Process genomic data with specified filters and criteria.')

    # Define arguments that the script can accept
    parser.add_argument("--work_dir", default='/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/ClinVar/', help="Working directory path")
    parser.add_argument("--submission_summary", default='submission_summary_20240127.txt.gz', help="Submission summary file")
    parser.add_argument("--variant_summary", default='variant_summary_20240127.txt.gz', help="Variant summary file")
    parser.add_argument("--variant_allele_relat", default='variation_allele_20240127.txt.gz', help="Variant Allele relationship file")
    parser.add_argument("--date", default='20240127', help="Date of ClinVar file")
    parser.add_argument("--gnomad", default='/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/eas_af-only-gnomad.genomes.v4.0.sites.mt', help="GnomAD file path")
    parser.add_argument("--ClinicalSignificance_Pathogenic", type=str, default='Pathogenic,Likely pathogenic,Pathogenic low penetrance,Likely pathogenic low penetrance,Established risk allele,Likely risk allele', help="Clinical significance pathogenic categories")
    parser.add_argument("--ClinicalSignificance_Benign", type=str, default='Uncertain significance,Uncertain risk allele,Benign,Likely benign', help="Clinical significance benign categories")
    parser.add_argument("--ReviewStatus_criteria", type=str, default='practice guideline,reviewed by expert panel,criteria provided, multiple submitters, no conflicts', help="Review status criteria")
    parser.add_argument("--AlleleType_filter", type=str, default='Complex,copy number loss,copy number gain,fusion,Inversion,Tandem duplication,Translocation', help="Allele type filters")
    parser.add_argument("--VariantType_filter", type=str, default='Haplotype,CompoundHeterozygote,Distinct chromosomes,Phase unknown,Diplotype', help="Variant type filters")
    parser.add_argument("--Phenotype_filter", type=str, default='C3661900:not provided,CN169374:not specified,na:See cases', help="Phenotype filters")
    parser.add_argument("--Chromosome_criteria", type=str, default=','.join([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']), help="Chromosome criteria")
    parser.add_argument("--drop_col", type=str, default='AlleleID,GeneID,HGNC_ID,ClinSigSimple,nsv/esv (dbVar),Origin,TestedInGTR,OtherIDs,Assembly,ChromosomeAccession,Start,Stop,Cytogenetic,ReferenceAllele,AlternateAllele,NumberSubmitters,SomaticClinicalImpact,SomaticClinicalImpactLastEvaluated,ReviewStatusClinicalImpact,Oncogenicity,OncogenicityLastEvaluated,ReviewStatusOncogenicity,ClinicalSignificance_list,ReviewStatus_list,PhenotypeList,SCV,ExplanationOfInterpretation', help="Drop variant submission merged df col")

    # Define and parse command-line arguments here as shown previously
    args = parser.parse_args()

    work_dir = args.work_dir
    submission_summary = args.work_dir + args.submission_summary
    variant_summary = args.work_dir + args.variant_summary
    variant_allele_relat = args.work_dir + args.variant_allele_relat
    gnomad = args.gnomad
    date = args.date
    output_path = work_dir + 'clinvar_summary_' + date + '.ht'
    
    # Convert the comma-separated strings into lists
    Chromosome_criteria = args.Chromosome_criteria.split(',')
    ClinicalSignificance_Pathogenic = args.ClinicalSignificance_Pathogenic.split(',')
    ClinicalSignificance_Benign = args.ClinicalSignificance_Benign.split(',')
    ReviewStatus_criteria = re.split(r',(?=\S)', args.ReviewStatus_criteria)
    AlleleType_filter = args.AlleleType_filter.split(',')
    VariantType_filter = args.VariantType_filter.split(',')
    Phenotype_filter = args.Phenotype_filter.split(',')
    drop_col = args.drop_col.split(',')
    main()