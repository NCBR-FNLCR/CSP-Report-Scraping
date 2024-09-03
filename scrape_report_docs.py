import sys
import os
import re
import datetime
import glob
import pandas as pd
import numpy as np
from labkey.api_wrapper import APIWrapper
from IPython.display import display
# import docx
from docx2python import docx2python
import argparse
from argparse import RawTextHelpFormatter
import tika
from tika import parser
tika.initVM()

# Returns sentence containing NCF1 (if there)
def parse_ncf1(file_name):
##### SSL issue from tika

    neg_rpt = parser.from_file(file_name)
    text = str(neg_rpt['content'])
    text = pd.Series(re.split('\n|\.', text))

    match = text[text.str.contains('NCF1', case = False, na = False)]
    if not match.empty:
        match = '; '.join(match)
    else:
        match = None

    return match


def is_empty(test):
    return pd.isna(test) or str(test).isspace() or len(str(test)) == 0


# Replaces A with T and C with G and vice versa to get flipped strand
def get_flipped_strand(dna_change):

    if pd.isna(dna_change):
        return None

    dna_change = str(dna_change)
    flipped_dna_change = ''

    for letter in dna_change:
        if letter == 'A':
            letter = 'T'
        elif letter == 'T':
            letter = 'A'
        elif letter == 'C':
            letter = 'G'
        elif letter == 'G':
            letter = 'C'

        flipped_dna_change = flipped_dna_change + letter

    return flipped_dna_change


# takes in final table and checks scraped variant classification w/ Seqr classification from Tags
def check_variant_classification(variants):

    variants.reset_index(inplace = True, drop = True)
    for i in range(0, variants.shape[0]):
        report_class = variants.at[i, 'Classification']
        seqr_class = variants.at[i, 'Tags']

        # get variant classification tag from seqr if available.
        if not is_empty(seqr_class):
            seqr_class = pd.Series(str(seqr_class).split('|'))
            seqr_class = seqr_class[seqr_class.str.contains('pathogenic|likely pathogenic|uncertain significance|likely benign|benign',
                                                            na=False, case=False)]

            # if able to compare, do so. Otherwise just continue to the next comparison
            if not seqr_class.empty and not is_empty(report_class):
                seqr_class = seqr_class.str.replace('VUS.*', 'VUS', regex=True)
                report_class = str(report_class)
                if 'uncertain' in report_class.lower():
                    report_class = 'VUS'

                if report_class.lower() not in seqr_class.str.lower().tolist():
                    new_note = str(variants.at[i, 'Comparison Notes']) + ' | Variant classification does not match'
                    new_note = new_note.strip(' |,')
                    variants.at[i, 'Comparison Notes'] = new_note

    return variants


# seqr_df passes in all variants that are shared across families.
# Returns subject level variant table. Each row is a phenotips ID and a variant that the person has
def build_subject_variant_tab(seqr_df, id_dict, fam_dict, mrn_dict, name_dict):
    subj_variants = pd.DataFrame(
        columns=['MRN', 'Patient_Name', 'Phenotips_ID', 'SEQR_ID', 'Phenotips_Fam_ID', 'SEQR_Fam_ID',
                 'hgvsc', 'hgvsp', 'Seqr_Gene', 'Current_HGNC_Gene', 'Chrom', 'Pos', 'Ref', 'Alt',
                 'Genotype', 'in_HGMD', 'Zygosity', 'Tags', 'Notes', 'Effect', 'sift', 'polyphen'])

    # loop thru all rows of SEQR which contains variants shared across families
    for ind, row in seqr_df.iterrows():
        # go thru all samples in seqr row
        for i in range(0, 10):
            colname = 'GT_genotype_' + str(i)
            samplename = 'sample_id_' + str(i)
            genotype = str(row[colname])

            seqr_id = str(row[samplename])
            seqr_id = seqr_id[:seqr_id.find('_')]
            phen_id = id_dict.get(seqr_id)
            mrn = mrn_dict.get(seqr_id)
            name = name_dict.get(seqr_id)

            # Don't allow blanks
            if '/' in genotype and genotype != './.':

                seqr_fam = str(row['family'])
                seqr_fam = seqr_fam[:seqr_fam.find('_')]
                phen_fam = fam_dict.get(seqr_fam)

                hom_alt = row['alt'] + '/' + row['alt']
                htz = row['ref'] + '/' + row['alt']
                hom_ref = row['ref'] + '/' + row['ref']

                zygosity = 'Unknown'

                # don't allow any hom ref subjects in the variants output
                if genotype == hom_ref:
                    continue
                elif genotype == hom_alt:
                    zygosity = 'Homozygous Alt'
                elif genotype == htz:
                    zygosity = 'Heterozygous'
                elif genotype[:genotype.find('/')] != row['ref'] and genotype[genotype.find('/') + 1:] != row['ref']:
                    zygosity = 'Compound Het'

                # populate individual level row from SEQR data into subj_variants
                subj_variants = subj_variants.append({
                    'MRN': mrn,
                    'Patient_Name': name,
                    'Phenotips_ID': phen_id,
                    'SEQR_ID': seqr_id,
                    'Phenotips_Fam_ID': phen_fam,
                    'SEQR_Fam_ID': seqr_fam,
                    'hgvsc': row['hgvsc'],
                    'hgvsp': row['hgvsp'],
                    'Seqr_Gene': row['gene'],
                    'Current_HGNC_Gene': row['latest_gene_name'],
                    'Chrom': row['chrom'],
                    'Pos': row['pos'],
                    'Ref': row['ref'],
                    'Alt': row['alt'],
                    'Genotype': genotype,
                    'in_HGMD': row['in_hgmd'],
                    'Zygosity': zygosity,
                    'Tags': row['tags'],
                    'Notes': row['notes'],
                    'Effect': row['effect'],
                    'sift': row['sift'],
                    'polyphen': row['polyphen']}, ignore_index=True, sort = False)

    subj_variants.sort_values(['Phenotips_ID'], inplace=True)

    return subj_variants


# Returns all rows from Seqr that are tagged with Inconclusive Negative, but are missing actual Negative CRIS reports
def check_inconclusive_neg_variants(variants, seqr, mrn_df):

    # makes dicts out of our mrn dictionary
    id_dict = dict(zip(mrn_df.Seqr_ID, mrn_df.Phenotips_ID))
    fam_dict = dict(zip(mrn_df.Seqr_Family_ID, mrn_df.Phenotips_Family_ID))
    mrn_dict = dict(zip(mrn_df.Seqr_ID, mrn_df.MRN))
    name_dict = dict(zip(mrn_df.Seqr_ID, mrn_df.Patient_Name))

    # keep only Inconclusive negative report tagged variants
    seqr = seqr[seqr['tags'].str.contains('Inconclusive negative report', na=False)]

    # Get subject level table for easier comparison
    seqr_individual = build_subject_variant_tab(seqr, id_dict, fam_dict, mrn_dict, name_dict)

    no_reports = pd.DataFrame(columns=['MRN', 'Phenotips_ID', 'SEQR_ID', 'Phenotips_Fam_ID', 'SEQR_Fam_ID', 'hgvsc',
                                       'hgvsp', 'Seqr_Gene', 'Current_HGNC_Gene', 'Chrom', 'Pos', 'Ref', 'Alt', 'Genotype',
                                       'Zygosity', 'Tags', 'Notes', 'Effect', 'sift', 'polyphen', 'in_HGMD'])

    # Checks each SEQR row against list of Negative reports
    for ind, row in seqr_individual.iterrows():
        mrn = row['MRN']

        # cut down variants to those that match the MRN
        have_neg_reports = variants[variants['MRN'] == mrn]

        # if no reports to be found, save the seqr row
        if have_neg_reports.empty:
            no_reports = no_reports.append(row, sort=False)

    return no_reports


# Returns all rows from Seqr that are tagged with CRIS Report, but are missing actual Positive CRIS reports
def check_cris_report_variants(variants, seqr, mrn_df):

    # makes dicts out of our mrn dictionary
    id_dict = dict(zip(mrn_df.Seqr_ID, mrn_df.Phenotips_ID))
    fam_dict = dict(zip(mrn_df.Seqr_Family_ID, mrn_df.Phenotips_Family_ID))
    mrn_dict = dict(zip(mrn_df.Seqr_ID, mrn_df.MRN))
    name_dict = dict(zip(mrn_df.Seqr_ID, mrn_df.Patient_Name))

    # Limit to variants with CRIS Report Tag
    seqr = seqr[seqr['tags'].str.contains('CRIS Report', na=False)]

    # Get subject level table for easier comparison
    seqr_individual = build_subject_variant_tab(seqr, id_dict, fam_dict, mrn_dict, name_dict)
    seqr_individual.to_csv('CRIS_Report_Seqr_Subject_Level.csv', index = False)

    no_reports = pd.DataFrame(columns=['MRN', 'Phenotips_ID', 'SEQR_ID', 'Phenotips_Fam_ID', 'SEQR_Fam_ID',
                                       'hgvsc', 'hgvsp', 'Seqr_Gene', 'Current_HGNC_Gene', 'Chrom', 'Pos', 'Ref', 'Alt',
                                       'Genotype', 'Zygosity', 'Tags', 'Notes', 'Effect', 'sift', 'polyphen', 'in_HGMD'])

    # Checks each SEQR row against list of CRIS Reports
    for ind, row in seqr_individual.iterrows():
        mrn = row['MRN']
        dna_change = row['hgvsc']

        # cut down variants to those that match the MRN and the DNA Change
        have_mrn = variants[variants['MRN'] == mrn]
        have_cris_reports = have_mrn[have_mrn['hgvsc'] == dna_change]

        # if no reports to be found, save the seqr row
        if have_cris_reports.empty:
            flipped_dna_change = get_flipped_strand(dna_change)
            have_cris_reports = have_mrn[have_mrn['hgvsc'] == flipped_dna_change]

        if have_cris_reports.empty:
            no_reports = no_reports.append(row, sort=False)

    return no_reports


# Copies info from SEQR into the variants table:
#   variants = df with all scraped variants
#   variant_ind = row in variants table to look for ID / assign variant data
#   seqr_row = SEQR row that has matched dna/protein change
def copy_info(variants, variant_ind, seqr_row):
    matched_id = False

    # Check all IDs in Seqr row for a match with the variant row ID
    for col in range(0, 10):
        sample_name = 'sample_id_' + str(col)
        seqr_id = str(seqr_row[sample_name])

        # Found a match! Adds seqr variant info to df
        if str(variants.at[variant_ind, 'Seqr_ID']) in seqr_id:
            matched_id = True
            variants.at[variant_ind, 'Seqr_Gene'] = seqr_row['gene']
            variants.at[variant_ind, 'Current_HGNC_Gene'] = seqr_row['latest_gene_name']
            variants.at[variant_ind, 'Chrom'] = seqr_row['chrom']
            variants.at[variant_ind, 'Pos'] = seqr_row['pos']
            variants.at[variant_ind, 'Ref'] = seqr_row['ref']
            variants.at[variant_ind, 'Alt'] = seqr_row['alt']
            variants.at[variant_ind, 'hgvsc'] = seqr_row['hgvsc']
            variants.at[variant_ind, 'hgvsp'] = seqr_row['hgvsp']
            variants.at[variant_ind, 'in_HGMD'] = seqr_row['in_hgmd']
            variants.at[variant_ind, 'Tags'] = seqr_row['tags']
            variants.at[variant_ind, 'Notes'] = seqr_row['notes']
            break

    return variants, matched_id


# Checks HGNC alias gene names. Returns True if the genes are aliases of each other
def are_alias_genes(gene1, gene2):

    if str(gene1) == str(gene2):
        return True

    # check if gene1 is the approved symbol
    gene_info = hgnc[hgnc['Approved symbol'] == str(gene1)]
    if gene_info.shape[0] == 1:

        # check if gene2 is in the older symbols
        gene2 = '|' + str(gene2) + '|'
        previous = gene_info[gene_info['Previous symbols'].str.contains(gene2, regex = False, na = False)]
        if not previous.empty:
            return True
        # check if gene2 is an alias
        else:
            alias = gene_info[gene_info['Alias symbols'].str.contains(gene2, regex = False, na = False)]
            if not alias.empty:
                return True
    # check if gene2 is the approved symbol
    else:
        gene_info = hgnc[hgnc['Approved symbol'] == str(gene2)]
        if gene_info.shape[0] == 1:

            # check if gene1 is in the older symbols
            gene1 = '|' + str(gene1) + '|'
            previous = gene_info[gene_info['Previous symbols'].str.contains(gene1, regex=False, na=False)]
            if not previous.empty:
                return True
            # check if gene1 is an alias
            else:
                alias = gene_info[gene_info['Alias symbols'].str.contains(gene1, regex=False, na=False)]
                if not alias.empty:
                    return True

    return False


# Adds previous gene names to Seqr spreadsheet
def add_previous_gene_names(seqr):

    # This column will contain the gene listed in seqr + all previous names
    seqr['latest_gene_name'] = None
    seqr.reset_index(inplace = True, drop = True)
    for i in range(0, seqr.shape[0]):
        latest_name = hgnc.get(str(seqr.at[i, 'gene']))
        if not pd.isna(latest_name):
            seqr.at[i, 'latest_gene_name'] = latest_name
        else:
            seqr.at[i, 'latest_gene_name'] = seqr.at[i, 'gene']


    return seqr


# Adds variants that match DNA Change from SEQR
def add_seqr_variants(variants, seqr):

    # Add seqr fields
    variants.insert(9, 'Seqr_Gene', None)
    variants.insert(10, 'Current_HGNC_Gene', None)
    variants['Chrom'] = None
    variants['Pos'] = None
    variants['Ref'] = None
    variants['Alt'] = None
    variants['hgvsc'] = None
    variants['hgvsp'] = None
    variants['in_HGMD'] = None
    variants['Tags'] = None
    variants['Notes'] = None
    variants['Comparison Notes'] = ''

    # Loop through each row in scraped variants
    for i in range(0, variants.shape[0]):

        # stores if family id exists in the seqr download (for diagnostic purposes)
        family_exists = seqr['family'].str.contains(str(variants.at[i, 'Seqr_Family_ID']), na=False).any()

        # Removes any spaces from rows that look like c.2164 G>A for example
        dna_change = str(variants.at[i, 'DNA Change'])
        if dna_change.startswith('c.'):
            dna_change = dna_change.replace(' ', '')
            variants.at[i, 'DNA Change'] = dna_change

        # Removes spaces from protein change field
        protein_change = str(variants.at[i, 'Protein Change'])
        if protein_change.startswith('p.'):
            protein_change = protein_change.replace(' ', '')
            variants.at[i, 'Protein Change'] = protein_change

        matched_id = False
        matched_dna_change = False
        matched_notes = False
        matched_protein = False
        fam_id = variants.at[i, 'Seqr_Family_ID']

        # get most recent hgnc name for variant
        report_gene = hgnc.get(str(variants.at[i, 'Gene']))
        if pd.isna(report_gene):
            report_gene = str(variants.at[i, 'Gene'])

        # Checks for Gene field match with Seqr (or previous hgnc gene name)
        # seqr_rows = seqr[seqr['possible_gene_names'].str.contains('|' + report_gene + '|', na = False, case = False, regex = False)]
        seqr_rows = seqr[seqr['latest_gene_name'] == report_gene]
        if seqr_rows.empty:
            variants.at[i, 'Comparison Notes'] = 'No Gene match'
            continue

        # At this point, seqr_rows contains all Seqr rows with a matching gene
        # Loop through each gene match seqr_row and check dna_change, notes, protein_change
        for ind, row in seqr_rows.iterrows():

            # checks IDs if the DNA Change is actually there
            if dna_change in str(row['hgvsc']):
                matched_dna_change = True
                variants, matched_id = copy_info(variants, i, row)
            else:
                flipped_dna_change = get_flipped_strand(dna_change)
                # see if flipped DNA change will match
                if flipped_dna_change in str(row['hgvsc']):
                    matched_dna_change = True
                    variants, matched_id = copy_info(variants, i, row)
                else:
                    # Try to see if the variant nomenclature shows up in the notes
                    if 'nomenclature' in str(row['notes']) or 'Nomenclature' in str(row['notes']):
                        if dna_change in str(row['notes']) or flipped_dna_change in str(row['notes']):
                            matched_notes = True
                            variants, matched_id = copy_info(variants, i, row)

                    # if not in the notes, try looking at protein change
                    else:
                        # if there's a protein change match, check if the DNA match field has same number of base deletion
                        if protein_change in str(row['hgvsp']):

                            num_bases = dna_change.split('del')

                            # get number of bases from DNA Change field and compare to seqr
                            if len(num_bases) > 1 and num_bases[1].isdigit():
                                num_bases = int(num_bases[1])
                                num_bases_seqr = row['hgvsc'].split('del')

                                # compares dna change base deletion length to seqr base deletion length
                                if len(num_bases_seqr) > 1:
                                    num_bases_seqr = len(num_bases_seqr[1])
                                    if num_bases == num_bases_seqr:
                                        matched_protein = True
                                        variants, matched_id = copy_info(variants, i, row)

            if matched_id:
                if matched_notes:
                    variants.at[i, 'Comparison Notes'] = 'Variant nomenclature was found in notes. Not an exact match.'
                elif matched_protein:
                    variants.at[i, 'Comparison Notes'] = 'Protein Change match found in Seqr. DNA Change was not an exact match.'
                break

        if not matched_id:
            if matched_dna_change:
                variants.at[i, 'Comparison Notes'] = 'Other subjects in SEQR have this variant, but not this subject.'
            elif matched_notes:
                variants.at[i, 'Comparison Notes'] = 'Variant nomenclature was found in notes in other subjects, but not this subject.'
            elif matched_protein:
                variants.at[i, 'Comparison Notes'] = 'Protein Change match found in Seqr for other subjects, but not this subject.'
            else:
                variants.at[i, 'Comparison Notes'] = 'No DNA Change match'

            if not family_exists:
                variants.at[i, 'Comparison Notes'] = 'Seqr Family ID not found in Seqr download. ' + str(
                    variants.at[i, 'Comparison Notes'])

    variants.rename(columns = {'Gene': 'Report Gene'}, inplace = True)

    return variants


# Returns MRN, SEQR, Phenotips ID dictionary
def make_id_dict():

    labkey_server = 'hgrepo.niaid.nih.gov'
    project_name = 'CSP'  
    contextPath = 'labkey'
    schema = 'lists'
    table = 'CSP Master'
    api = APIWrapper(labkey_server, project_name, contextPath, api_key=,use_ssl=True)

    result = api.query.select_rows(schema, table)
    master_table = pd.DataFrame(result["rows"])
    gen_master_table = master_table[['GRISClinicalID','CRIMSON_MRN','CRIS_Patient_Name','PT_Phenotips_Family_ID','Batch_Received_Genome','Genome_ID','Missing_Active_Order_Combined']]
    exo_master_table = master_table[['GRISClinicalID','CRIMSON_MRN','CRIS_Patient_Name','PT_Phenotips_Family_ID','Batch_Received_Exome','Exome_ID','Missing_Active_Order_Combined']]
    gen_master_table = gen_master_table[~gen_master_table['Batch_Received_Genome'].isna()]
    exo_master_table = exo_master_table[~exo_master_table['Batch_Received_Exome'].isna()]
    gen_master_table = gen_master_table[gen_master_table['Missing_Active_Order_Combined']=='No']
    exo_master_table = exo_master_table[exo_master_table['Missing_Active_Order_Combined']=='No']


    gen_master_table.rename(columns={'GRISClinicalID':'Phenotips_ID','CRIMSON_MRN':'MRN','CRIS_Patient_Name':'Patient_Name','PT_Phenotips_Family_ID':'Phenotips_Family_ID',
                                 'Batch_Received_Genome':'Batch_Received' , 'Genome_ID':'Genome ID'}, inplace=True)
    exo_master_table.rename(columns={'GRISClinicalID':'Phenotips_ID','CRIMSON_MRN':'MRN','CRIS_Patient_Name':'Patient_Name','PT_Phenotips_Family_ID':'Phenotips_Family_ID',
                                 'Batch_Received_Exome':'Batch_Received','Exome_ID':'Exome ID' }, inplace=True)
    gen_master_table['Exome ID']=''
    exo_master_table['Genome ID']=''
    mrn_dict = pd.concat([gen_master_table, exo_master_table])
                         
    mrn_dict = mrn_dict.dropna(subset=['Exome ID','Genome ID'])
    mrn_dict = mrn_dict[['MRN', 'Patient_Name', 'Phenotips_ID', 'Phenotips_Family_ID', 'Batch_Received']]
    mrn_dict.columns = ['MRN', 'Patient_Name', 'Phenotips_ID', 'Phenotips_Family_ID', 'Batch_Received']
    mrn_dict.drop_duplicates(inplace=True)


    api = APIWrapper(labkey_server, project_name, contextPath, api_key=,use_ssl=True)
    Seqr_results = api.query.select_rows(schema_name='lists',query_name='Seqr Map')
    Seqr_table = pd.DataFrame(result["rows"])
    Seqr_table = Seqr_table[['GRISClinicalID', 'seqr_Individual_ID','seqr_Family_ID']]
    Seqr_table.rename(columns={'GRISClinicalID':'Phenotips_ID','seqr_Individual_ID':'Seqr_ID',
                                 'seqr_Family_ID':'Seqr_Family_ID' }, inplace=True)
    Seqr_table.drop_duplicates(inplace=True)
    mrn_dict = mrn_dict.merge(Seqr_table, on='Phenotips_ID', how='left')
    mrn_dict.to_csv('dict.csv')

    return mrn_dict


def get_secondary_findings(variants):

    ### here Nov 5
    

    variants= pd.read_csv('Variants.csv')
    variants  = variants.dropna(subset=['Batch_Received'])
    secondary = variants[variants['Tags'].str.contains('Secondary finding', na = False)]
    secondary = secondary[secondary['Tags'].str.contains('CRIS Report', na = False)]
    secondary['to_remove'] = ''
   
    # Loop through all secondary findings reports
    # and check if the words "Secondary findings: Not detected." are in the report
    # if so, remove from the secondary findings list, and put a note in Comparison Notes that
    # the report has been incorrectly tagged in Seqr for this person
    for i, row in secondary.iterrows():

        received = secondary.at[i, 'Batch_Received']
        mrn = secondary.at[i, 'MRN']
        mrn = mrn.replace("-", "")
        batch = int(re.search(r'\d+', received).group())
        
        
        if batch > 121:     ####need to update this everytime!!!!
           continue
        batch_name = 'Batch ' + str(batch)
        print(batch_name)

        report_types = pd.Series(glob.glob('/Volumes/NIAID_Centralized_Sequencing/REPORTS/FINAL/' + batch_name + '/*'))
        #report_dir = report_types[report_types.str.contains('Positive', case=False)]  
        ### NOV 2021
        #report_types = pd.DataFrame(glob.glob('/Volumes/NIAID_Centralized_Sequencing/REPORTS/FINAL/' + batch_name + '/*'))
        if len(report_types) == 0:
            continue
        for sub_dir in report_types:
            if 'POSITIVE' in sub_dir.upper():
                report_dir = sub_dir
        
        #report_dir =report_types[report_types.iloc[:,0].str.contains('POSITIVE', case=False)].loc[:,0].to_list()
        #report_dir = report_dir[0]

        report_path = report_dir + "/" + mrn + "*.docx"

        fnames = glob.glob(report_path)
        fnames = [x for x in fnames if "~$" not in x]
        if len(fnames) == 0:
            continue
            
        doc = docx2python(fnames[0])
        text = doc.text.lower()

        if "secondary findings: not detected" in text:
            if variants.at[i, 'Comparison Notes'] is None:
                variants.at[i, 'Comparison Notes'] = "Tagged in Seqr with secondary finding, but not detected in actual patient report."
            else:
                variants.at[i, 'Comparison Notes'] = str(variants.at[i, 'Comparison Notes']) + ", Tagged in Seqr with secondary finding, but not detected in actual patient report."
            secondary.at[i, 'to_remove'] = 'Remove'

    secondary = secondary[~secondary['to_remove'].str.match('Remove', na = False)]
    secondary.drop(columns = ['to_remove'], inplace = True)

    return secondary, variants


def main():

    global hgnc


    hgnc = pd.read_csv('hgnc_name_dict.csv')  
    hgnc = dict(zip(hgnc['Previous Symbol'], hgnc['Current Symbol']))

    print('REMEMBER: Add colors to the header for Report vs. SEQR columns!')

    #### Usage statement ####
    parsestr = 'Scrapes all positive reports for variant information.\n\
            Includes Negative reports if flagged. Outputs a CSV.\n\n\
            Usage:\n\
                scrape_report_docs.py -o outdir -s seqr_files_dir \n\n'

    #parser = argparse.ArgumentParser(description=parsestr, formatter_class=RawTextHelpFormatter)
    #parser.add_argument('-o', '--outdir', required=True, type=str, help='Output file directory')
    #parser.add_argument('-s', '--seqr_dir', required=True, type=str, help='Directory where Seqr downloads are stored')

    #args = parser.parse_args()
    outdir = 
    seqr_dir = 

    #add / to paths to use directories
    if seqr_dir.isspace() or seqr_dir == '':
        seqr_dir = ''
    elif not seqr_dir.endswith("/"):
        seqr_dir = seqr_dir + "/"

    if outdir.isspace() or outdir == '':
        outdir = ''
    elif not outdir.endswith("/"):
        outdir = outdir + "/"

    ####### Make a combined seqr file from all the files in seqr_dir #######
    print('Combining Seqr files')
    seqr_files = glob.glob(seqr_dir + '*.csv')

    seqr_combined = pd.DataFrame()
    for file in seqr_files:
        curr_table = pd.read_csv(file, dtype = str)
        seqr_combined = seqr_combined.append(curr_table)

    seqr_combined.drop_duplicates(inplace = True)
    seqr_combined.to_csv('seqr_combined.csv', index = False)

    seqr_combined = add_previous_gene_names(seqr_combined)
    seqr_combined.to_csv('seqr_combined_hgnc.csv', index = False)


    flip_reports = 3
    neg_patients = pd.DataFrame(columns=['MRN', 'Batch'])
    print('Including negative reports.')

    #mrn_dict = make_id_dict()
    #mrn_dict.to_csv(outdir + 'dict.csv', index=False)
    mrn_dict = pd.read_csv('dict.csv')
    ##### NOv 2021 mrn dict doesnt have info for WGS 
    # mrn_dict = pd.read_csv(outdir + 'dict.csv', dtype=str)

    variants = pd.DataFrame(columns = ['MRN', 'Zygosity', 'DNA Change', 'RefSeq_ID', 'Disease Inheritance',
                                       'Gene', 'Classification', 'Associated Disease', 'Protein Change', 'OMIM'])
    cma = pd.DataFrame(columns = ['MRN', 'CMA Result', 'Associated Disease', 'Classification'])

    # Make sure we're scraping all available batches in the FINAL directory
    batch_dirs = glob.glob('/Volumes/NIAID_Centralized_Sequencing/REPORTS/FINAL/Batch*')

    # add the extra categories
    batch_dirs = batch_dirs + ['/Volumes/NIAID_Centralized_Sequencing/REPORTS/FINAL/RAPIDTAT',
                               '/Volumes/NIAID_Centralized_Sequencing/REPORTS/FINAL/Sanger only'] #was named GeneDx folder, now is RAPITAT



    print('Scraping Reports...')

    num_positive = 0
    num_negative = 0

    unscraped = pd.DataFrame(columns = ['Filename', 'MRN', 'Batch'])


    # Loop through all report directories
    for batch_dir in batch_dirs:

        batch_name = batch_dir.replace('/Volumes/NIAID_Centralized_Sequencing/REPORTS/FINAL/', '')
        
 
        print('Scraping ' + batch_name)

        report_types = pd.Series(glob.glob(batch_dir + '/*'))

        # Flips between positive and negative reports
        for i in range(1, flip_reports):

            # Switch between positive and negative reports
            if i == 1:
                report_dir = report_types[report_types.str.contains('Positive', case = False)]
                if report_dir.empty:
                    continue
                report_dir = report_dir.reset_index(drop = True)[0]

                dir_name_doc = report_dir + '/*.docx'
                fnames = glob.glob(dir_name_doc)
                fnames = [x for x in fnames if "~$" not in x]
                num_positive = num_positive + len(fnames)
            else:
                report_dir = report_types[report_types.str.contains('Negative', case = False)]
                if report_dir.empty:
                    continue
                report_dir = report_dir.reset_index(drop = True)[0]
                dir_name_pdf = report_dir + '/*.pdf'
                fnames = glob.glob(dir_name_pdf)
                num_negative = num_negative + len(fnames)

            # Loop through all files in the directory
            for file in fnames:

                # Get MRN from the file name
                mrn = file.replace(report_dir + '/', '')
                mrn = mrn[0:7]
                mrn = "-".join([mrn[a:a + 2] for a in range(0, len(mrn), 2)])

                # Set to true, unless the report is actually scraped
                unable_to_scrape = True
                doc = None

                should_print = mrn == ''

                # If we're doing Positive Reports
                if i == 1:
                    try:
                        doc = docx2python(file)
                    except:
                        # can't scrape so add to list of unscraped
                        print('Cannot scrape file! ' + str(file))
                        unscraped = unscraped.append(pd.DataFrame({'Filename': file, 'MRN': mrn, 'Batch': batch_name}, index = [0]))
                        continue

                    # Split up doc text by periods or new lines
                    text_list = re.split('\. |\n', doc.text)
                    text_list = pd.Series(text_list)
                    text_list = text_list[~text_list.str.isspace()]
                    text_list = text_list[text_list != '']
                    text_list = text_list.reset_index(drop=True)

                    # Check to see if this is a CMA Positive report
                    if not text_list[text_list.str.contains('COPY NUMBER VARIANT')].empty:

                        unable_to_scrape = False

                        dna_change = text_list[text_list.str.contains('CMA', case = False)]
                        dna_change = dna_change[dna_change.str.contains('showed')
                                                | dna_change.str.contains('identified')
                                                | dna_change.str.contains('detected')]
                        dna_change = dna_change.str.strip()
                        dna_change = ';'.join(dna_change)
                        dna_change = dna_change.replace('CMA detected copy number variant', '')
                        dna_change = dna_change.strip(' ;:')

                        cma_df = pd.DataFrame({'DNA_Change': dna_change, 'MRN': mrn}, index = [0])


                        disease = text_list[text_list.str.contains('OMIM', case = False)]
                        disease = disease[~disease.str.contains('OMIM \*')]
                        disease = disease.str.strip()
                        disease.rename('Associated_Disease')
                        disease = disease.to_frame()
                        disease['MRN'] = mrn


                        pathogenicity = ''


                        benign = text_list[text_list == 'Benign']
                        for index, word in benign.iteritems():
                            if index + 5 < len(text_list):
                                test = text_list[index:index+5].tolist()
                                if test == ['Benign', 'Likely Benign', 'Uncertain Significance', 'Likely Pathogenic', 'Pathogenic']:
                                    text_list.drop(index = list(range(index, index+5)), inplace = True)
                                    break

                        classification = text_list[text_list.str.contains('pathogenic|likely pathogenic|uncertain significance|likely benign|benign',
                                                                          na = False, case = False)]
                        classification = classification.str.strip()
                        classification.reset_index(drop = True, inplace = True)

                        if not classification.empty:


                            line = classification[0].lower()
                            if 'likely pathogenic' in line:
                                pathogenicity = 'Likely Pathogenic'
                            elif 'pathogenic' in line:
                                pathogenicity = 'Pathogenic'
                            elif 'likely benign' in line:
                                pathogenicity = 'Likely Benign'
                            elif 'benign' in line:
                                pathogenicity = 'Benign'
                            elif 'uncertain significance' in line:
                                pathogenicity = 'Uncertain Significance'
                            elif 'inconclusive' in line:
                                pathogenicity = 'Inconclusive'

                        cma_df = cma_df.merge(disease, on = 'MRN', how = 'left')
                        cma_df['Classification'] = pathogenicity
                        cma_df.rename(columns = {'0_x': 'CMA Result', '0_y': 'Associated Disease'}, inplace = True)
                        cma_df = cma_df.reindex(columns = ['MRN', 'CMA Result', 'Associated Disease', 'Classification'])
                        cma = cma.append(cma_df, sort = False)


                    counter = 1

                    for table in doc.body:
                        df = pd.DataFrame(table)

                        if 'Gene' in df.iloc[0, 0]:

                            unable_to_scrape = False
                            counter = counter + 1

                            num_rows = df.shape[0]
                            num_cols = df.shape[1]

                            for row in range(0, num_rows):
                                for col in range(0, num_cols):
                                    df.at[row, col] = ' '.join(df.at[row, col])

                            df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
                            df = df.apply(lambda x: x.str.replace('- ', '-') if x.dtype == "object" else x)
                            df = df.apply(lambda x: x.str.replace('  ', '-') if x.dtype == "object" else x)

                            if should_print:
                                df.to_csv(outdir + mrn + '_raw' + str(counter) + '.csv', index=False)

                            # Assign first row to column names
                            df.columns = df.loc[0]
                            df = df.loc[1:]
                            df.reset_index(drop=True, inplace=True)


                            if 'DNA Change' not in df.columns:

                                if 'Benign' in df.columns:
                                    continue
                                else:
                                    df.insert(0, 'MRN', mrn)
                                    df.rename(columns = {'Variant': 'DNA Change', 'Inheritance': 'Disease Inheritance'}, inplace = True)
                                    df['Protein Change'] = None

                                    for j in range(0, df.shape[0]):

                                        dna_change = df.at[j, 'DNA Change'].split(' ')
                                        if df.at[j, 'DNA Change'].startswith('chr'):
                                            df.at[j, 'DNA Change'] = dna_change[3]
                                            df.at[j, 'Protein Change'] = dna_change[4].strip('()')

                            else:
                                if 'Gene' in df.columns:

                                    df_cols = pd.Series(df.columns)
                                    if df_cols.str.contains('OMIM', case = False, na = False).any():
                                        omim_col = df_cols[df_cols.str.contains('OMIM', case = False, na = False)]
                                        omim_col = omim_col.tolist()[0]
                                        df.rename(columns={omim_col: 'OMIM'}, inplace=True)
                                    else:
                                        df['OMIM'] = ''

                                    if 'Associated disease' in df.columns:
                                        df.rename(columns={'Associated disease': 'Associated Disease'}, inplace=True)
                                    if 'Disease inheritance' in df.columns:
                                        df.rename(columns={'Disease inheritance': 'Disease Inheritance'}, inplace=True)
                                    if 'Inheritance' in df.columns:
                                        df.rename(columns={'Inheritance': 'Disease Inheritance'}, inplace=True)

                                    if df.shape[0] > 1:


                                        for row_ind in range(1, df.shape[0]):

                                            curr_gene = df.at[row_ind, 'Gene']
                                            if should_print:
                                                print('GENE:')
                                                print('|' + str(curr_gene) + '|')
                                                print(pd.isna(curr_gene))
                                                print(str(curr_gene).isspace())
                                                print(len(str(curr_gene)) == 0)

                                            if pd.isna(curr_gene) or str(curr_gene).isspace() or len(str(curr_gene)) == 0:
                                                df.at[row_ind, 'Gene'] = df.at[row_ind - 1, 'Gene']

                                                # if gene was missing, there's a good chance this info is also missing. Use row above
                                                if is_empty(df.at[row_ind, 'DNA Change']):
                                                    df.at[row_ind, 'DNA Change'] = df.at[row_ind - 1, 'DNA Change']
                                                if is_empty(df.at[row_ind, 'Protein Change']):
                                                    df.at[row_ind, 'Protein Change'] = df.at[row_ind - 1, 'Protein Change']
                                                if is_empty(df.at[row_ind, 'Zygosity']):
                                                    df.at[row_ind, 'Zygosity'] = df.at[row_ind - 1, 'Zygosity']
                                                if is_empty(df.at[row_ind, 'Classification']):
                                                    df.at[row_ind, 'Classification'] = df.at[row_ind - 1, 'Classification']

                                            if 'Associated Disease' not in df.columns.to_list():
                                            	print("File: "+file)
                                            	continue
                                            	
                                            if 'Associated Disease' in df.columns.to_list():
                                                if df.at[row_ind, 'Associated Disease'] == '':

                                                    prev_genes = df.at[row_ind - 1, 'Gene'].split(', ')
                                                    curr_genes = df.at[row_ind, 'Gene'].split(', ')
                                                    shared_gene = list(set(prev_genes) & set(curr_genes))

                                                    if len(shared_gene) > 0:
                                                        if is_empty(df.at[row_ind, 'Associated Disease']):
                                                            df.at[row_ind, 'Associated Disease'] = df.at[row_ind - 1, 'Associated Disease']
                                                        if is_empty(df.at[row_ind, 'OMIM']):
                                                            df.at[row_ind, 'OMIM'] = df.at[row_ind - 1, 'OMIM']
                                                        if is_empty(df.at[row_ind, 'Disease Inheritance']):
                                                            df.at[row_ind, 'Disease Inheritance'] = df.at[row_ind - 1, 'Disease Inheritance']
                                                        if is_empty(df.at[row_ind, 'Classification']):
                                                            df.at[row_ind, 'Classification'] = df.at[row_ind - 1, 'Classification']

                                        
                                        variant_groups = df.groupby(['Gene', 'DNA Change', 'Zygosity', 'Classification'])
                                        group_sizes = variant_groups.size()
                                        group_sizes = group_sizes[group_sizes > 1]

                                        if not group_sizes.empty:

                                            new_df = pd.DataFrame()
                                            for group in variant_groups:
                                                new_row = group[1]

                                                if new_row.shape[0] > 1:
                                                    new_row.reset_index(inplace=True, drop=True)

                                                    disease_info = new_row[['Associated Disease', 'OMIM',
                                                                            'Disease Inheritance']].drop_duplicates()
                                                    disease_info.fillna('', inplace = True)

                                                    assoc_disease = '; '.join(disease_info['Associated Disease'].tolist()).strip('; ')
                                                    omim_id = '; '.join(disease_info['OMIM'].tolist()).strip('; ')
                                                    disease_inher = '; '.join(disease_info['Disease Inheritance'].tolist()).strip('; ')

                                                    new_row = new_row.loc[0]
                                                    new_row['Associated Disease'] = assoc_disease
                                                    new_row['OMIM'] = omim_id
                                                    new_row['Disease Inheritance'] = disease_inher

                                                new_df = new_df.append(new_row, sort = False, ignore_index = True)
                                            df = new_df

                                        if should_print:
                                            df.to_csv(outdir + mrn + '_final.csv', index=False)

                                df.insert(0, 'MRN', mrn)
                                # df.to_csv(outdir + mrn + '_' + str(counter) + '.csv', index=False)


                            df['RefSeq_ID'] = None

                            for ind, row in df.iterrows():

                                fixed_dna_change = row['DNA Change']
                                refseq_id = text_list[text_list.str.contains(row['Gene'], na = False, regex = False)]
                                refseq_id = refseq_id[refseq_id.str.contains(row['DNA Change'], na = False, regex = False)]
                                refseq_id = refseq_id[refseq_id.str.contains('NM_', na = False, regex = False)]

                        
                                if refseq_id.empty:
                                    fixed_dna_change = ''.join(str(row['DNA Change']).split())
                                    refseq_id = text_list[text_list.str.contains(row['Gene'], na=False, regex = False)]
                                    refseq_id = refseq_id[refseq_id.str.contains(fixed_dna_change, na=False, regex = False)]
                                    refseq_id = refseq_id[refseq_id.str.contains('NM_', na=False, regex = False)]

                                if not refseq_id.empty:
                                    refseq_id = refseq_id.str.findall('NM_\d+')
                                    refseq_id = pd.Series([element for list_ in refseq_id for element in list_]).unique().tolist()
                                    refseq_id = ','.join(refseq_id)
                                    df.at[ind, 'RefSeq_ID'] = refseq_id
                                    df.at[ind, 'DNA Change'] = fixed_dna_change

                            variants = variants.append(df, sort = False, ignore_index = True)

                else:
                    unable_to_scrape = False
                    neg_rpt = pd.DataFrame(columns = ['MRN', 'Batch'])
                    neg_rpt['MRN'] = [mrn]
                    neg_rpt['Batch'] = [batch_name]

                    neg_date = file.replace(report_dir + '/', '')
                    neg_date = re.findall('\d{8}', neg_date)
                    if len(neg_date) == 1:
                        neg_date = str(neg_date[0])

                   
                        if neg_date.startswith('20'):
                            neg_date = pd.to_datetime(neg_date, format = '%Y%m%d')
                  
                        else:
                            neg_date = pd.to_datetime(neg_date, format = '%m%d%Y')
                    else:
                        neg_date = None

                    neg_rpt['Date'] = [neg_date]

 

                    neg_patients = neg_patients.append(neg_rpt, sort = False)

                if unable_to_scrape:
                    unscraped = unscraped.append(pd.DataFrame({'Filename': file, 'MRN': mrn, 'Batch': batch_name}, index = [0]))

    print('Scraped all reports.')
    print(num_positive)
    print(num_negative)

    variants.to_csv(outdir + 'all_variants.csv', index = False)
    neg_patients.to_csv(outdir + 'neg_rpt_patients.csv', index = False)
    unscraped.to_csv(outdir + 'unscraped.csv', index = False)

    variants.drop_duplicates(inplace = True)
    variants = variants.merge(mrn_dict, on='MRN', how='left')
    variants = variants[['Phenotips_ID', 'Seqr_ID', 'Phenotips_Family_ID', 'Seqr_Family_ID', 'Batch_Received', 'MRN', 'Patient_Name', 'RefSeq_ID',
                         'Gene', 'DNA Change', 'Protein Change', 'Zygosity', 'Classification', 'Associated Disease', 'OMIM', 'Disease Inheritance']]
    variants.to_csv(outdir + 'all_variants_ids.csv', index=False)
    neg_patients = neg_patients.merge(mrn_dict, on = 'MRN', how = 'left')

    cma = cma.merge(mrn_dict, on = 'MRN', how = 'left')
    cma.drop_duplicates(inplace = True)
    cma = cma[['Phenotips_ID', 'Seqr_ID', 'Phenotips_Family_ID', 'Seqr_Family_ID',
               'MRN', 'Patient_Name', 'CMA Result', 'Associated Disease', 'Classification']]
    cma.to_csv(outdir + 'all_cma.csv', index = False)
    # cma = pd.read_csv(outdir + 'all_cma.csv', dtype = str)

    print('Matching SEQR variants with CRIS Report scraped variants...')
    # variants_merged = pd.read_csv(outdir + 'all_variants_seqr.csv', dtype = str)
    variants_merged = add_seqr_variants(variants, seqr_combined)
    variants_merged.to_csv(outdir + 'all_variants_seqr.csv', index = False)

    print('Finding "CRIS Report" tagged variants that are missing Positive Reports...')
    # variants_merged = pd.read_csv(outdir + 'all_variants_seqr.csv', dtype=str)
    missing_positive_rpt = check_cris_report_variants(variants_merged, seqr_combined, mrn_dict)

    print('Finding "Inconclusive negative report" tagged variants that are missing Negative Reports...')
    missing_negative_rpt = check_inconclusive_neg_variants(neg_patients, seqr_combined, mrn_dict)

    print('Double checking secondary findings...')
    secondary_findings, variants_merged = get_secondary_findings(variants_merged)

    print('Checking variant classification between reports and Seqr...')
    variants_merged = check_variant_classification(variants_merged)

    missing_in_seqr = variants_merged[variants_merged['Chrom'].isnull()]
    matched_in_seqr = variants_merged[~variants_merged['Chrom'].isnull()]

    # Write out all dfs into one Excel spreadsheet
    writer = pd.ExcelWriter(outdir + 'scraping_results.xlsx', engine = 'xlsxwriter')
    variants_merged.to_excel(writer, sheet_name = 'All Scraped Pos Rpt Variants', index = False)
    matched_in_seqr.to_excel(writer, sheet_name = 'Pos Rpt variants in SEQR', index = False)
    missing_in_seqr.to_excel(writer, sheet_name = 'Pos Rpt variants not in SEQR', index = False)
    missing_positive_rpt.to_excel(writer, sheet_name = 'CRIS Rpt variants no Pos Rpt', index = False)
    unscraped.to_excel(writer, sheet_name = 'Cannot scrape Pos Rpt', index = False)
    neg_patients.to_excel(writer, sheet_name = 'All patients w Negative Reports', index = False)
    missing_negative_rpt.to_excel(writer, sheet_name = 'Inconc Neg variants no Neg Rpt', index = False)
    cma.to_excel(writer, sheet_name = 'CMA Pos Rpts', index = False)
    secondary_findings.to_excel(writer, sheet_name = 'Secondary Findings', index = False)

    writer.save()


if __name__ == '__main__':
    main()
