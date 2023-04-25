import csv

def concatenate_columns(input_file, output_file):
    with open(input_file, "r") as file_in, open(output_file, "w") as file_out:
        reader = csv.reader(file_in, delimiter="\t")
        writer = csv.writer(file_out, delimiter="\t")

        # Read the header row separately
        header = next(reader)
        # Remove the suffix "_Chimeric.out.junction" from the sample names
        header = [col.replace("_Chimeric.out.junction", "") for col in header]
        writer.writerow(header)

        for row in reader:
            if len(row) >= 3:
                # Concatenate first three columns with : and -
                chr_start_end = f"{row[0]}:{row[1]}-{row[2]}"
                # Replace the first three columns with the concatenated value
                row[:3] = [chr_start_end]
            writer.writerow(row)


concatenate_columns("CircRNACount","test.tsv")
