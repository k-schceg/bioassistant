# BioAssistant

<img align="right" src="bioassistant.jpg" alt="BioAssistant" width="400">

**BioAssistant** is a toolkit for basic operations with DNA and RNA sequences. Your bioassistant will help you with such operations on DNA and RNA as transcription, finding the reversed, complementary and reversed complementary sequence, determining the gc-content and belonging of the sequence to a palindrome. It can also filter fastq by gc-content, sequence length and read quality.
## Content

* [Installation](#installation)
* [Examples](#examples)
* [FAQ](#faq)
* [Discussion and Development](#discussion-and-development)

## Installation

Clone repository
~~~
git clone git@github.com:k-schceg/bioassistant.git
~~~
Move contents to your project folder
~~~
cp -r bioassistant/* /path/to/your/project/
~~~
Import bioassistant module into your script
```python
import bioassistant
```
Enjoy


## Examples

* #### Transcription
```python
>>> bioassistant.run_dna_rna_tools("ATCG", "GATAATGGC", "transcribe")
['AUCG', 'GAUAAUGGC']
```
* #### GC_content
```python
>>> bio.filter_fastq(
      {
          '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGA', 
                         'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGE'),
          '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAG', 
                         'BFFFFFFFB@B@A<@D>BDDAC'),
          '@SRX079803': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTT', 
                         'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFB')
      }, 
      gc_bounds=(100),
      length_bounds=(0, 25),
      quality_threshold=17
  )
{'@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAG', 'BFFFFFFFB@B@A<@D>BDDAC')}
```

## FAQ

**Q** What updates to expect in the new version?

**A** It will be possible to work directly with fasta and fasta files.

**Q** Can I process several sequences at once?

**A** Yes, bioassistant can work with a list of sequences at once

**Q** What scale is used to assess the quality of a read?

**A** BioAssistant uses the phred33 scale.

## Discussion and Development

You can send your feedback directly to the GitHub [issue tracker](https://github.com/k-schceg/bioassistant/issues)

All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome.

