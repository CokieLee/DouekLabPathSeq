import csv
def write_to_sample_sheet(sample_sheet_name, project_type, project_id, flow_cell_ID, lane, sample_ID, species, adaptor, index, phenotype, chain):
    single_sample_row = [project_type, project_id, flow_cell_ID, lane, sample_ID, species, adaptor, index, phenotype, chain]
    with open(sample_sheet_name, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(single_sample_row)
