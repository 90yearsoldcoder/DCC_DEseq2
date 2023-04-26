import csv
import subprocess


def gene_bed_file(input_file, output_file):
    # Open the input TSV file for reading
    with open(input_file, 'r', newline='') as infile:

        # Open the output TSV file for writing
        with open(output_file, 'w', newline='') as outfile:
        
            # Create a CSV reader and writer object
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')
            
            # Skip the first row(header)
            next(reader)
            
            # Loop over each remaining row in the input file
            for row in reader:
                
                # Extract the first three columns
                new_row = row[:3]
                
                # Write the new row to the output file
                writer.writerow(new_row)
                
    print("Extraction complete. Results saved to", output_file)


def extract_last_column(input_file):
    #extract the annotated gene name
    reference = {}
    
    # Open the input TSV file for reading
    with open(input_file, 'r', newline='') as infile:
    
        # Create a CSV reader object
        reader = csv.reader(infile, delimiter='\t')
        
        # Loop over each row in the input file
        for row in reader:
            
            # Extract the last column and add it to the list
            reference[f"{row[0]}:{row[1]}-{row[2]}"] = row[-1]
            
    # Return the list of last column values
    return reference

def concatenate_columns(input_file, output_file,reference):
    with open(input_file, "r") as file_in, open(output_file, "w") as file_out:
        reader = csv.reader(file_in, delimiter="\t")
        writer = csv.writer(file_out, delimiter="\t")

        # Read the header row separately
        header = next(reader)
        # Remove the suffix "_Chimeric.out.junction" from the sample names
        header = [col.replace("_Chimeric.out.junction", "") for col in header]
        header[:3] = ["gene"]
        writer.writerow(header)

        for row in reader:
            if len(row) >= 3:
                # Concatenate first three columns with : and -
                chr_start_end = f"{row[0]}:{row[1]}-{row[2]}(circ_" + reference.get(f"{row[0]}:{row[1]}-{row[2]}", "Not Anot") + ")"
                # Replace the first three columns with the concatenated value
                row[:3] = [chr_start_end]
            writer.writerow(row)


#generate the bed file for further annotation
gene_bed_file("/restricted/projectnb/ncrna/minty/samb_data4/Run/DCC/SP_E4_2_2_ver2/CircRNACount","../Processing/CircRNA.bed")

# Run the bedtool intersect script
subprocess.call(['bash', 'anotate.sh'])

reference = extract_last_column("../Processing/CircRNA_anotated.bed")

#print(reference)

concatenate_columns("/restricted/projectnb/ncrna/minty/samb_data4/Run/DCC/SP_E4_2_2_ver2/CircRNACount","../Processing/CircRNACount_prepared.tsv", reference)
