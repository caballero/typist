# Typist

This is a simple script to compute a class enrichment for cell types using gene expression information (bulk RNA-seq).

## Requirements
- Only Python 3.X

## Usage

    > python src/typist.py -h
    usage: typist.py [-h] -i INPUT_FILE -g GENEMARKERS_FILE -o OUTPUT_FILE [-d DELIMITER] [-m MIN_EXPRESSION] [-a] [-n] [-v]

    options:
    -h, --help            show this help message and exit
    -i INPUT_FILE, --input_file INPUT_FILE
                            input file
    -g GENEMARKERS_FILE, --genemarkers_file GENEMARKERS_FILE
                            gene markers file
    -o OUTPUT_FILE, --output_file OUTPUT_FILE
                            output file
    -d DELIMITER, --delimiter DELIMITER
                            field delimiter
    -m MIN_EXPRESSION, --min_expression MIN_EXPRESSION
                            filter by minimum expression value
    -a, --average_filter  use average expression filter
    -n, --cpm_normalization
                            use normalization (CPM by default)
    -v, --verbose         increase output verbosity

## Inputs

Two files are required: markers table and expression table.

- Gene makers should be listed with genes as rows, first column should be the gene id (which should match the expression table) and columns indicate different cell-type classes with cell values indicating if the gene is a marker, this can be specified as a binary value (1 = marker, 0 no marker) or weighted values.

Example:

|Gene|ClassA|ClassB|ClassC|
|----|-----:|-----:|-----:|
|GenA|1|0|0|
|GenB|1|0|0|
|GenC|0|1|0|
|GenD|0|1|0|
|GenE|0|0|1|

- Gene expressions should be also a table with gene as rows and samples as columns.

Example:

|Gene|Sample1|Sample2|Sample3|
|----|-----:|-----:|-----:|
|GenA|100|10|0|
|GenB|120|0|6|
|GenC|2|230|2|
|GenD|3|220|8|
|GenE|0|0|1000|

## Output
A simple table with the scores for each sample in each  cell-type category.

## How is the scores computed?
After normalization and filtering, gene markers are searched in samples, if the gene is expressed (based on minimal expression), the score adds the value from the markers table. 
Final score value is generated dividing the total score over the total expected score.


(C) 
Juan Caballero 2025