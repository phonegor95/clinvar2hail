import pandas as pd
import hail as hl
import gzip
import argparse

temp_dir = '/mnt/beegfs/tmp'
hl.init(
  master='local[26]',
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
def custom_agg_review_status(series):
    if 'practice guideline' in series.values:
        return ['practice guideline']
    elif 'reviewed by expert panel' in series.values:
        return ['reviewed by expert panel']
    else:
        return list(series)

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

def process_genomic_data(work_dir, submission_summary, variant_summary, gnomad,
                         ClinicalSignificance_Pathogenic, ClinicalSignificance_Benign,
                         ReviewStatus_criteria, VariantType_filter, Phenotype_filter,
                         Chromosome_criteria):
    # Function implementation goes here
    # You can add your data processing logic here, using the provided parameters
    print("Processing genomic data with the following parameters:")
    print(f"Work directory: {work_dir}")
    # Add more print statements or actual processing logic as needed
    pass

#work_dir = '/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/ClinVar/'
#submission_summary = work_dir + 'submission_summary_20240127.txt.gz'
#variant_summary = work_dir + 'variation_allele_20240127.txt.gz'
#gnomad = '/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/eas_af-only-gnomad.genomes.v4.0.sites.mt'
#ClinicalSignificance_Pathogenic = ['Pathogenic', 'Likely pathogenic', 'Pathogenic, low penetrance', 'Likely pathogenic, low penetrance', 'Established risk allele', 'Likely risk allele']
#ClinicalSignificance_Benign = ['Uncertain significance', 'Uncertain risk allele', 'Benign', 'Likely benign']
#ReviewStatus_criteria = ['practice guideline', 'reviewed by expert panel', 'criteria provided, multiple submitters, no conflicts']
#VariantType_filter = ['Complex', 'copy number loss', 'copy number gain', 'fusion', 'Inversion', 'Tandem duplication', 'Translocation']
#Phenotype_filter = ['C3661900:not provided', 'CN169374:not specified', 'na:See cases']
#Chromosome_criteria = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']

# Create the parser
parser = argparse.ArgumentParser(description='Process genomic data with specified filters and criteria.')

# Define arguments that the script can accept
parser.add_argument("--work_dir", default='/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/ClinVar/', help="Working directory path")
parser.add_argument("--submission_summary", default='submission_summary_20240127.txt.gz', help="Submission summary file")
parser.add_argument("--variant_summary", default='variant_summary_20240127.txt.gz', help="Variant summary file")
parser.add_argument("--date", default='20240127', help="Date of ClinVar file")
parser.add_argument("--gnomad", default='/mnt/beegfs/hongyf/igenomes_base/Homo_sapiens/GATK/GRCh38/Annotation/GermlineResource/eas_af-only-gnomad.genomes.v4.0.sites.mt', help="GnomAD file path")
parser.add_argument("--ClinicalSignificance_Pathogenic", type=str, default='Pathogenic,Likely pathogenic,Pathogenic low penetrance,Likely pathogenic low penetrance,Established risk allele,Likely risk allele', help="Clinical significance pathogenic categories")
parser.add_argument("--ClinicalSignificance_Benign", type=str, default='Uncertain significance,Uncertain risk allele,Benign,Likely benign', help="Clinical significance benign categories")
parser.add_argument("--ReviewStatus_criteria", type=str, default='practice guideline,reviewed by expert panel,criteria provided multiple submitters no conflicts', help="Review status criteria")
parser.add_argument("--VariantType_filter", type=str, default='Complex,copy number loss,copy number gain,fusion,Inversion,Tandem duplication,Translocation', help="Variant type filters")
parser.add_argument("--Phenotype_filter", type=str, default='C3661900:not provided,CN169374:not specified,na:See cases', help="Phenotype filters")
parser.add_argument("--Chromosome_criteria", type=str, default=','.join([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']), help="Chromosome criteria")

# Parse the arguments
args = parser.parse_args()

work_dir = args.work_dir
submission_summary = args.work_dir + args.submission_summary
variant_summary = args.work_dir + args.variant_summary
gnomad = args.gnomad
date = args.date

# Split the comma-separated strings into lists
ClinicalSignificance_Pathogenic = args.ClinicalSignificance_Pathogenic.split(',')
ClinicalSignificance_Benign = args.ClinicalSignificance_Benign.split(',')
ReviewStatus_criteria = args.ReviewStatus_criteria.split(',')
VariantType_filter = args.VariantType_filter.split(',')
Phenotype_filter = args.Phenotype_filter.split(',')
Chromosome_criteria = args.Chromosome_criteria.split(',')


header_line_variant, skip_rows_variant = get_header_line(variant_summary)
variant_summary_df = pd.read_csv(variant_summary, sep='\t', comment=None, skiprows = skip_rows_variant, names = header_line_variant, low_memory=False)
variant_summary_df = variant_summary_df[variant_summary_df['Assembly'] == 'GRCh38']
variant_summary_df = variant_summary_df[variant_summary_df['Chromosome'].isin(Chromosome_criteria)]
variant_summary_df.replace({'Chromosome': {'MT': 'M'}}, inplace=True)
variant_summary_df['Chromosome'] = variant_summary_df['Chromosome'].apply(lambda x: 'chr' + str(x))
variant_summary_df = variant_summary_df[variant_summary_df['PositionVCF'] != -1]
variant_summary_df = variant_summary_df[variant_summary_df['ReferenceAlleleVCF'] != 'na']
variant_summary_df = variant_summary_df[variant_summary_df['AlternateAlleleVCF'] != 'na']
variant_summary_df = variant_summary_df[variant_summary_df['ReviewStatus'].isin(ReviewStatus_criteria)]
variant_summary_df = variant_summary_df[variant_summary_df['ClinicalSignificance'].str.contains('|'.join(ClinicalSignificance_Pathogenic), case=True, na=False)]
variant_summary_df = variant_summary_df[~variant_summary_df['Type'].isin(VariantType_filter)]
variant_summary_df = variant_summary_df.drop_duplicates(subset='VariationID', keep='first')
print(f"variant_summary_df: {variant_summary_df}")

header_line_submission, skip_rows_submission = get_header_line(submission_summary)
submission_summary_df = pd.read_csv(submission_summary, sep='\t', comment=None, skiprows = skip_rows_submission, names = header_line_submission, low_memory=False)
submission_summary_df = submission_summary_df[~submission_summary_df['ReportedPhenotypeInfo'].isin(Phenotype_filter)]
submission_summary_df = filter_review_status_to_priority(submission_summary_df)

# Step 1: Group by 'VariationID' and 'ReportedPhenotypeInfo' and aggregate the other columns
submission_summary_df = submission_summary_df.groupby(['VariationID', 'ReportedPhenotypeInfo']).agg({
    'ClinicalSignificance': lambda x: list(x),
    #'DateLastEvaluated': lambda x: list(x),
    'Description': lambda x: list(x),
    #'SubmittedPhenotypeInfo': lambda x: x.value_counts().to_dict(),
    'ReviewStatus': lambda x: list(x),
    #'CollectionMethod': lambda x: x.value_counts().to_dict(),
    #'OriginCounts': lambda x: x.value_counts().to_dict(),
    'Submitter': lambda x: list(x),
    'SCV': lambda x: list(x),
    #'SubmittedGeneSymbol': lambda x: x.value_counts().to_dict(),
    'ExplanationOfInterpretation': lambda x: list(x)
}).reset_index()


# Subset the DataFrame based on the criteria
submission_summary_df = submission_summary_df[
    submission_summary_df['ClinicalSignificance'].apply(
    lambda x: any(value in ClinicalSignificance_Pathogenic for value in x) and not any(value in ClinicalSignificance_Benign for value in x)
    )
].reset_index(drop=True)

submission_summary_df = submission_summary_df[
    submission_summary_df['ReviewStatus'].apply(meets_review_criteria)
].reset_index(drop=True)

# Step 2: Group by 'VariationID' to collect all different 'ReportedPhenotypeInfo' entries
submission_summary_df = submission_summary_df.groupby('VariationID').agg({
    'ReportedPhenotypeInfo': lambda x: list(x),
    'ClinicalSignificance': lambda x: list(x),
    #'DateLastEvaluated': lambda x: list(x),
    'Description': lambda x: list(x),
    #'SubmittedPhenotypeInfo': lambda x: list(x),
    'ReviewStatus': lambda x: list(x),
    #'CollectionMethod': lambda x: list(x),
    #'OriginCounts': lambda x: list(x),
    'Submitter': lambda x: list(x),
    'SCV': lambda x: list(x),
    #'SubmittedGeneSymbol': lambda x: list(x),
    'ExplanationOfInterpretation': lambda x: list(x)
}).reset_index()
print(f"variant_summary_df: {variant_summary_df}")

variant_submission_summary_df = submission_summary_df.merge(variant_summary_df, on='VariationID')
print(f"variant_submission_summary_df: {variant_submission_summary_df}")

# Step 3: Convert the final grouped DataFrame to a Hail Table
variant_submission_summary_ht = hl.Table.from_pandas(variant_submission_summary_df[['VariationID', 'Name', 'Type', 'GeneSymbol', 'OriginSimple', 'Chromosome', 'PositionVCF', 'ReferenceAlleleVCF', 'AlternateAlleleVCF', 'ClinicalSignificance_y', 'ReviewStatus_y', 'ClinSigSimple', 'LastEvaluated', 'ReportedPhenotypeInfo', 'SCV', 'Submitter', 'Description', 'ExplanationOfInterpretation']])
del variant_summary_df
del submission_summary_df
del variant_submission_summary_df
variant_submission_summary_ht = variant_submission_summary_ht.transmute(locus = hl.locus(variant_submission_summary_ht.Chromosome, variant_submission_summary_ht.PositionVCF, reference_genome='GRCh38'))
variant_submission_summary_ht = variant_submission_summary_ht.transmute(alleles = [variant_submission_summary_ht.ReferenceAlleleVCF, variant_submission_summary_ht.AlternateAlleleVCF])
variant_submission_summary_ht = variant_submission_summary_ht.key_by(variant_submission_summary_ht.locus, variant_submission_summary_ht.alleles)
gnomad_ht = hl.read_table(gnomad)
variant_submission_summary_ht = variant_submission_summary_ht.annotate(AF_global = gnomad_ht.rows()[variant_submission_summary_ht.key].info.AF)
variant_submission_summary_ht = variant_submission_summary_ht.annotate(AF_eas = gnomad_ht.rows()[variant_submission_summary_ht.key].info.AF_eas)
variant_submission_summary_ht.describe()
variant_submission_summary_ht.show()
variant_submission_summary_ht.write(work_dir + 'clinvar_summary_' + date +'_w_gnomad.ht', overwrite=True)