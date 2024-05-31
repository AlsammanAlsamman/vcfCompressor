# vcfCompressor

`vcfCompressor` is a command-line tool designed for compressing and reconstructing Variant Call Format (VCF) files. It provides efficient handling of large genomic data sets, enabling storage and retrieval in a compact binary format.

## Features

- **Compress VCF files**: Converts VCF files into a binary format along with sample names and marker information files.
- **Reconstruct VCF files**: Restores the original VCF files from the binary format and associated metadata files.
- **Reduced file size**: Achieves significant compression ratios compared to standard compression methods like ZIP.
- **Keeps only GT field**: Retains only the Genotype (GT) field for each sample, reducing the file size.
- **Fast and efficient**: Optimized for large genomic data sets with minimal memory usage and processing time.
- **Command-line interface**: Simple and intuitive commands for compressing and reconstructing VCF files.
- **Cross-platform compatibility**: Works on Windows, macOS, and Linux operating systems.
- **Open-source**: Available under the MIT License for free use and modification.
- **Customizable**: Easily adaptable to different VCF file formats and data structures.
- **Scalable**: Handles thousands of samples and millions of markers with ease.

## Usage

- Please keep in mind that this tool is designed for reducing the size of VCF files and may not be suitable for all use cases.
- It's not a filtering tool, and it focuses on compressing the data while maintaining the essential information for downstream analysis.
- It can be used after filtering and processing the VCF files to reduce storage requirements and improve data transfer speeds.
- If you don't care the DP, GQ, and other fields, this tool can be useful for you. As it only keeps the GT field, it can significantly reduce the file size.
- You can combine this tool with standard compression methods like ZIP for further reduction in file size.
  For example, you can compress the binary file generated by `vcfCompressor` using ZIP to achieve better compression ratios.

### General Usage

```sh
vcfCompressor -c <vcf_filename>
vcfCompressor -r -b <binary_file> -m <marker_file> -s <sample_file> -o <output_vcf>
vcfCompressor -h | --help
```

### Compress a VCF file

```sh
vcfCompressor -c input.vcf -ns 554 -nm 2388
```

### Reconstruct a VCF file

```sh
vcfCompressor -r -b input.bin -m input_markers.txt -s input_samples.txt -o output.vcf -ns 554 -nm 2388
```

### Display Help

```sh
vcfCompressor -h
```

## Compression Ratios

To calculate the compression ratio, you can use the formula:

```math
\text{Compression Ratio} = \frac{\text{Original Size}}{\text{Compressed Size}} 
```

Given the sizes:
- Original size: 445.8 MB
- Zip size: 9.2 MB
- vcfCompressor size: 27.1 MB
- vcfCompressor + Zip size: 547.4 KB

### Zip Compression Ratio
```math
 \text{Zip Compression Ratio} = \frac{445.8 \, \text{MB}}{9.2 \, \text{MB}} = 48.46 
```
### vcfCompressor Compression Ratio
```math
\text{vcfCompressor Compression Ratio} = \frac{445.8 \, \text{MB}}{27.1 \, \text{MB}} = 16.45 
```
### vcfCompressor + Zip Compression Ratio

Since the vcfCompressor + Zip size is in kilobytes, we need to convert it to megabytes for a consistent comparison:

```math
 547.4 \, \text{KB} = 0.5474 \, \text{MB}

 \text{vcfCompressor + Zip Compression Ratio} = \frac{445.8 \, \text{MB}}{0.5474 \, \text{MB}} = 814.54 
```
### Summary

- **Zip Compression Ratio**: 48.46
- **vcfCompressor Compression Ratio**: 16.45
- **vcfCompressor + Zip Compression Ratio**: 814.54


## Building the Program

### Compiler

Ensure you have `gcc` installed.

### Compilation

To compile the program, use:

```sh
make
```

To run the program:

```sh
make run
```

To clean up object files and binaries:

```sh
make clean
```

To test the program:

```sh
make test
```

## License

This project is licensed under the MIT License.

## Citation
Ther is no publication yet. Please cite this repository if you use this tool.

## Authors

- Alsamman M. Alsamman

## Email
- A.Alsamman[at]cgiar.org
- samman.mohammed[at]gmail.com
- smahmoud[at]ageri.sci.eg

## Acknowledgments

- Any acknowledgments you want to include