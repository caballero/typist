#!/usr/bin/env python3

import csv
import logging
import argparse

def get_markers(file, delim):    
    """
    Reads a gene markers file and returns a dictionary of markers and a list of category names.
    
    The file is expected to be a CSV/TSV file with the first column being the gene
    marker and the first row containing the category names. The other columns are values
    associated with each marker (1 for expressed, 0 for non expressed).
    
    Parameters:
        file (str): The name of the file to be read.
        delim (str): The delimiter used in the file.
    
    Returns:
        tuple: A tuple where the first element is a dictionary of markers and their
        associated values, the maximal scores per category and the second element is a list of
        categories names.
    """
    markers = {}
    max_score = {}
    logging.debug("Reading markers file %s", file)
    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter=delim)
        colnames = [colname.strip() for colname in next(reader)]
        for col in colnames:
            max_score[col] = 0

        for row in reader:
            gen = row[0]
            markers[gen] = {}
            for col in range(1, len(row)):
                markers[gen][colnames[col]] = float(row[col])
                max_score[colnames[col]] += float(row[col])
    colnames.pop(0)

    logging.debug("Found markers: %s", str(len(markers)))
    logging.debug("Found categories: %s", colnames)
    
    return markers, max_score, colnames

def cpm_normalization(expressions):
    """
    Normalizes expression data to counts per million (CPM).
    
    Parameters:
        expressions (dict): A dictionary where keys are gene names and values are expression values.
    
    Returns:
        dict: A dictionary where keys are gene names and values are normalized expression values.
    """
    logging.debug("Normalizing expression data")
    cpm = {}
    for sample in expressions:
        total = sum(expressions[sample].values())
        for gene in expressions[sample]:
            cpm[sample][gene] = expressions[sample][gene] / total * 1000000
    return cpm

def get_expressions(file, delim, min_expr, cpm_norm, average_filter):
    """
    Reads an expression file and returns a dictionary of expressions and a list of sample names.
    
    The file is expected to be a CSV/TSV file with the first column being the gene and the 
    first row containing the sample names. The other columns are valuesassociated with each gene
    (expression values).
    
    Parameters:
        file (str): The name of the file to be read.
        delim (str): The delimiter used in the file.
        min_expr (float): The minimum expression value, below which a value is considered 0.
        no_norm (bool): If True, do not normalize the expression values.
        average_filter (bool): If True, filter the expression values by average expression in each sample.
    
    Returns:
        tuple: A tuple where the first element is a dictionary of expressions and their
        associated values, and the second element is a list of sample names.
    """
    logging.debug("Reading expressions file %s", file)
    expressions = {}

    with open(file, 'r') as f:
        reader = csv.reader(f, delimiter=delim)
        colnames = [colname.strip() for colname in next(reader)]
        for col in range(1, len(colnames)):
            expressions[colnames[col]] = {}

        for row in reader:
            gen = row[0]
            for col in range(1, len(row)):
                expressions[colnames[col]][gen] = float(row[col])

        colnames.pop(0)
    
    logging.debug("Found samples: %s", colnames)
    logging.debug("Found genes per sample: %s", str(len(expressions[colnames[0]])))

    if cpm_norm:
        expressions = cpm_normalization(expressions)
    
    if average_filter:
        logging.debug("Filtering by average expression")
        for sample in expressions:
            avg = sum(expressions[sample].values()) / len(expressions[sample])
            for gene in expressions[sample]:
                if expressions[sample][gene] < avg:
                    expressions[sample][gene] = 0
                else:
                    expressions[sample][gene] = 1
        return expressions, colnames

    else:
        logging.debug("Filtering by minimum expression")
        for sample in expressions:
            for gene in expressions[sample]:
                if expressions[sample][gene] < min_expr:
                    expressions[sample][gene] = 0
                else:
                    expressions[sample][gene] = 1
        return expressions, colnames


def get_predictions(expressions, samples, markers, max_scores, categ):
    """
    Makes predictions based on expression data.

    Parameters:
        expressions (dict): A dictionary of gene expression data, where keys are gene names and values are dictionaries of sample names and expression values.
        samples (list): A list of sample names.
        markers (dict): A dictionary of gene markers, where keys are gene names and values are dictionaries of category names and marker values.
        categ (list): A list of category names.

    Returns:
        dict: A dictionary of predictions, where keys are sample names and values are dictionaries of category names and prediction values.
    """
    logging.debug("Making predictions")
    predicts = {}
    for sample in samples:
        predicts[sample] = {}
        total_markers = 0
        for category in categ:
            predicts[sample][category] = 0
        for gene in expressions[sample]:
            if gene in markers:
                total_markers += 1
                for category in categ:
                    if expressions[sample][gene] >= markers[gene][category]:
                        predicts[sample][category] += markers[gene][category]
        
        for category in categ:
            if total_markers > 0:
                predicts[sample][category] = predicts[sample][category] / max_scores[category]
            else:
                logging.error("No markers found for sample %s", sample)
    return predicts
                    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True, help="input file")
    parser.add_argument("-g", "--genemarkers_file", required=True, help="gene markers file")
    parser.add_argument("-o", "--output_file", required=True, help="output file")
    parser.add_argument("-d", "--delimiter", default="\t", help="field delimiter")
    parser.add_argument("-m", "--min_expression", type=float, default=0.0, help="filter by minimum expression value")
    parser.add_argument("-a", "--average_filter", action="store_true", help="use average expression filter")
    parser.add_argument("-n", "--cpm_normalization", action="store_true", help="use normalization (CPM by default)")
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()

    input_file = args.input_file
    genemarkers_file = args.genemarkers_file
    output_file = args.output_file
    delimiter = args.delimiter
    min_expression = args.min_expression
    cpm_normalization = args.cpm_normalization
    average_filter = args.average_filter
    verbose = args.verbose

    if verbose:
        logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

    markers, max_scores, categ = get_markers(genemarkers_file, delimiter)

    expressions, samples = get_expressions(input_file, delimiter, min_expression, cpm_normalization, average_filter)

    predictions = get_predictions(expressions, samples, markers, max_scores, categ)

    logging.debug("Writing predictions to file %s", output_file)
    with open(output_file, 'w') as out:
        writer = csv.writer(out, delimiter=delimiter)
        writer.writerow(["Sample"] + categ)
        for sample in predictions:
            row = [sample]
            for category in categ:
                row.append(predictions[sample][category])
            writer.writerow(row)
        
    
if __name__ == "__main__":
    main()