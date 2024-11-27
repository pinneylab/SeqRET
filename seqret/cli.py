# cli.py #################################
import argparse
import csv
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description='Optimize DNA sequences for expression.')
    subparsers = parser.add_subparsers(dest='command')

    # Subparser for the default 'optimize' command
    parser_optimize = subparsers.add_parser('optimize', help='Optimize sequences using CLI')
    parser_optimize.add_argument('-s', '--sequence', type=str, help='DNA sequence to optimize')
    #parser_optimize.add_argument('-f', '--file', type=str, help='File containing DNA sequence to optimize')
    parser_optimize.add_argument('-i', '--input_csv', type=str, help='CSV file containing DNA sequences to optimize. Looks for "Sequence" column.')
    parser_optimize.add_argument('-o', '--output_csv', type=str, help='Output CSV file to write optimized sequences. Will include an "OptimizedSequence" column.')
    parser_optimize.add_argument('--filters', type=str, nargs='+', help='List of filters to apply (Not yet implemented.)')

    # Subparser for the 'webapp' command
    parser_webapp = subparsers.add_parser('webapp', help='Launch the web application')
    parser_webapp.add_argument('--host', type=str, default='127.0.0.1', help='Host for the web app')
    parser_webapp.add_argument('--port', type=int, default=8050, help='Port for the web app')

    args = parser.parse_args()
    # Set filters to None (run all of them), because we haven't implemented the filters argument yet
    args.filters = None

    if args.command == 'optimize':
        if args.sequence:
            sequence = args.sequence.upper().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
            optimize_and_print(sequence, args.filters)
        # elif args.file:
        #     with open(args.file, 'r') as f:
        #         sequence = f.read().upper().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
        #     optimize_and_print(sequence, args.filters)
        elif args.input_csv:
            if not args.output_csv:
                print('Error: Output CSV file must be specified when using input CSV.')
                return
            optimize_csv(args.input_csv, args.output_csv, args.filters)
        else:
            parser_optimize.print_help()
    elif args.command == 'webapp':
        from .app import start_app 
        start_app(host=args.host, port=args.port)
    else:
        parser.print_help()

def optimize_and_print(sequence, filters_titles):
    from .filters import get_filters
    from .optimizer import optimize_sequence

    filters_config = get_filters()
    if filters_titles:
        filters_enabled = [filter_conf['title'] in filters_titles for filter_conf in filters_config]
    else:
        filters_enabled = [True] * len(filters_config)
    optimized_sequence = optimize_sequence(sequence, filters_config=filters_config, filters_enabled=filters_enabled)
    print('Original sequence:')
    print(sequence)
    print('Optimized sequence:')
    print(optimized_sequence)

def optimize_csv(input_csv, output_csv, filters_titles):
    from .filters import get_filters
    from .optimizer import optimize_sequence

    filters_config = get_filters()
    if filters_titles:
        filters_enabled = [filter_conf['title'] in filters_titles for filter_conf in filters_config]
    else:
        filters_enabled = [True] * len(filters_config)
    with open(input_csv, 'r') as csvfile_in, open(output_csv, 'w', newline='') as csvfile_out:
        reader = csv.DictReader(csvfile_in)
        fieldnames = reader.fieldnames + ['OptimizedSequence']
        writer = csv.DictWriter(csvfile_out, fieldnames=fieldnames)
        writer.writeheader()
        for row in tqdm(reader):
            sequence = row['Sequence'].upper().replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
            optimized_sequence = optimize_sequence(sequence, filters_config=filters_config, filters_enabled=filters_enabled)
            row['OptimizedSequence'] = optimized_sequence
            writer.writerow(row)

if __name__ == '__main__':
    main()
