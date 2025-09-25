import math
import gzip
import argparse

## python filter_score_V2.py input.fq.gz output_pass.fq.gz output_fail.fq.gz
parser = argparse.ArgumentParser(description="Filter FASTQ reads")
parser.add_argument("input_path", help="Path to input FASTQ file (e.g., all.fq.gz)")
parser.add_argument("output_pass_path", help="Path to output file for passing reads")
parser.add_argument("output_fail_path", help="Path to output file for failing reads")

# Parse the arguments
args = parser.parse_args()

# Assign variables from arguments
input_path = args.input_path
output_pass_path = args.output_pass_path
output_fail_path = args.output_fail_path

def mean_qscore(qstring: str) -> float:
    """Calculates mean Q-score from a Q-string"""
    # Truncate string is sufficiently long
    qstr = qstring[60:] if len(qstring) > 60 else qstring
    # Add safeguard for empty strings
    if len(qstr) == 0:
        print("The length of the quality string is: ", len(qstring))
        print("Check out the read: ", qstring)
        print("Warning: Quality string is empty after truncation.")
        return 0.0
    # Convert ASCII to Phred quality scores, then to estimated error
    errors = [10**(-(ord(char) - 33)/10) for char in qstr]
    # Calculate mean error
    mean_error = sum(errors) / len(errors)
    # Convert back to q-score
    mean_qscore = -10 * math.log10(mean_error)
    return mean_qscore


with gzip.open(input_path, "rt") as infile, gzip.open(output_pass_path, "wt") as output_pass, gzip.open(output_fail_path, "wt") as output_fail:
    while True:
        header = infile.readline().strip()
        if not header:
            break  # End of file
        sequence = infile.readline().strip()
        plus = infile.readline().strip()
        quality_string = infile.readline().strip()

        all_high_quality = mean_qscore(quality_string) > 10 
        
        if all_high_quality:
            output_pass.write(f"{header}\n{sequence}\n{plus}\n{quality_string}\n")
        if all_high_quality==False:
            output_fail.write(f"{header}\n{sequence}\n{plus}\n{quality_string}\n")


print("done")

