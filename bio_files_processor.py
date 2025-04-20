import os
from datetime import datetime
import logging


logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=datetime.now().strftime("%Y_%m_%d-%H_%M_%S_%f_bio_files_processor_log.txt"),
    encoding='utf-8',
    level=logging.NOTSET,
    format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %I:%M:%S %p'
)

def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None, mode: str = None
) -> None:
    """Function convert_multiline_fasta_to_oneline

    Args:
        input_fasta: str
            path to input fasta file,
        output_fasta: str
            if None, default file name will be passed in format
            '%Y_%m_%d-%H_%M_%S_%f_convert_multiline_fasta_to_oneline.txt',
        mode: str = None
            if 'overwrite', function will overwrite output file
            if 'append', function will append into output file
            if None, function will raise error if file already exists
    Returns: None
        writes corrected sequences to output fasta file
    """
    if output_fasta is None:
        output_fasta = datetime.now().strftime(
            "%Y_%m_%d-%H_%M_%S_%f_convert_multiline_fasta_to_oneline.txt"
        )
    if os.path.isfile(output_fasta):
        if mode is None:
            logger.error(f"File {output_fasta} already exists! Change file name or 'mode' argument.")
            raise ValueError(
                f"File {output_fasta} already exists! Change file name or 'mode' argument."
            )
        elif mode == "overwrite":
            os.remove(output_fasta)
    with open(input_fasta, "r") as input_file:
        with open(output_fasta, "a") as output_file:
            flag_first_row = True
            for line in input_file:
                if not line.startswith(">"):
                    line = line.strip()
                elif not flag_first_row:
                    line = "\n" + line
                output_file.write(line)
                flag_first_row = False


def parse_blast_output(
    input_blast: str, output_blast: str = None, mode: str = None
) -> None:
    """Function parse_blast_output

    Args:
        input_blast: str
            path to input blast file,
        output_blast: str = None
            path to output blast file
            if None, default file name will be passed in format '%Y_%m_%d-%H_%M_%S_%f_parse_blast_output.txt',
        mode: str = None
            if 'overwrite', function will overwrite output file
            if 'append', function will append into output file
            if None, function will raise error if file already exists

    Returns: None
        writes sorted first row from each description to output blast file
    """
    logger.info('Start blast parsing')
    if not os.path.isfile(input_blast):
        logger.error(f"Input file {input_blast} doesn't exists!")
        raise ValueError(f"Input file {input_blast} doesn't exists!")
    if output_blast is None:
        output_blast = datetime.now().strftime(
            "%Y_%m_%d-%H_%M_%S_%f_parse_blast_output.txt"
        )
        logger.warning(f'Output blast file is not provided => set output_blast = "{output_blast}"')
    if os.path.isfile(output_blast):
        logger.warning(f'Output blast file "{output_blast}" already exists')
        if mode is None:
            logger.error(f"File {output_blast} already exists! Change file name or 'mode' argument.")
            raise ValueError(
                f"File {output_blast} already exists! Change file name or 'mode' argument."
            )
        elif mode == "overwrite":
            logger.info(f'Output file overwrite mode activated')
            os.remove(output_blast)
    first_rows = []
    flag_read_row = False
    logger.info(f'Processing blast...')
    with open(input_blast, "r") as input_file:
        for line in input_file:
            if line.startswith("Description"):
                flag_read_row = True
            elif flag_read_row:
                first_rows.append(line)
                flag_read_row = False
    logger.info(f'Blast file processed')
    first_rows.sort()
    logger.info(f'Saving to output file')
    with open(output_blast, "a") as output_file:
        for line in first_rows:
            output_file.write(line)
    logger.info(f'Success!')
