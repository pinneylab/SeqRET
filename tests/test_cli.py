# Here, we simply verify that the code is working by running a CSV,
# then checking that the 'OptimizedSequence' column matches the 'ExpectedOptimizedSequence' column.

import csv
from pathlib import Path

from seqret.cli import optimize_csv

TEST_DATA = Path(__file__).parent / "test_data"


def test_optimize_csv(tmp_path):
    input_csv = TEST_DATA / "test_medium.csv"
    output_csv = tmp_path / "test_medium_out.csv"
    optimize_csv(str(input_csv), str(output_csv), None)

    with open(output_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            assert row['OptimizedSequence'] == row['ExpectedOptimizedSequence']

