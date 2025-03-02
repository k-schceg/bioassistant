# BioAssistant 3.0

<img align="right" src="bioassistant.jpg" alt="BioAssistant" width="400">

**BioAssistant** is a toolkit for basic operations with DNA and RNA sequences. Your bioassistant will help you with such operations on DNA and RNA as transcription, finding the reversed, complementary and reversed complementary sequence. It can also filter fastq by gc-content, sequence length and read quality. 
Have a fasta file, but can't make a blast? BioAssistant will remove unnecessary gaps in the sequence. Got too much unnecessary information after multiple alignments from blast? No problem! BioAssistant will help you leave only the most important alignments.

Just call and BioAssistant will help you!
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
>>> dna = bioassistant.DNASequence(GATAATGGC)
dna.transcribe
'GAUAAUGGC'
```
* #### Fastq filter
```python
>>> bioassistant.filter_fastq(
      './examples/example_fastq.fastq', 
      './examples/example_fastq_out.fastq'
      gc_bounds=(100),
      length_bounds=(0, 25),
      quality_threshold=17)
```
* #### Parse blast results
```python
>> bioassistant.parse_blast_output(example_blast_results.txt, example_blast_results_out.txt)
```
Examples of input and output data are in the folder "examples"

## FAQ

**Q** What updates have appeared in the new version?

**A** Working with sequences has become even easier! Now sequences are objects of certain classes: DNASequence, RNASequence or AminoAcidSequence

**Q** Can I process several sequences at once?

**A** Yes, bioassistant can work with a list of sequences at once

**Q** What scale is used to assess the quality of a read?

**A** BioAssistant uses the phred33 scale.

## Discussion and Development

You can send your feedback directly to the GitHub [issue tracker](https://github.com/k-schceg/bioassistant/issues)

All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome.
