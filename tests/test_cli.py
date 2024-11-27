# Here, we simply verify that the code is working by running a CSV,
# then checking that the 'OptimizedSequence' column matches the 'ExpectedOptimizedSequence' column.

from seqret.cli import optimize_csv
import csv
import os

def test_optimize_csv():
    input_csv = 'test_data/test_medium.csv'
    output_csv = 'test_data/test_medium_out.csv'
    optimize_csv(input_csv, output_csv, None)
    
    with open(output_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            assert row['OptimizedSequence'] == row['ExpectedOptimizedSequence']
    
    os.remove(output_csv)

# Run the test
test_optimize_csv()
