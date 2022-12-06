import sys
import os
import click
import logging
from typing import *
import subprocess
import pathlib


logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)

@click.command(short_help="Runs the star build a genome index.")
@click.argument("fasta",
                type=click.Path(exists=True,
                                dir_okay=True,
                                file_okay=True,
                                readable=True,
                                resolve_path=True))
@click.argument("gtffile",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.argument("output",
                type=click.Path(exists=False,
                                file_okay=False,
                                dir_okay=True,
                                readable=False,
                                resolve_path=False))
@click.option("--threads", "-t",
              help="Number of threads to be used for genome generation",
              default=10)
@click.option("--length", "-l",
              help="The length of the genomic sequence around the annotated junction to be used in constructing te splice junctions database",
              default=100)

def build_star_index(fasta: str, gtffile: str,
                     output: str, threads: int,
                     length: int) -> None:
    output = pathlib.Path(output)
    if not output.exists():
        os.makedirs(output)
    commands = f'STAR --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {fasta} --sjdbGTFfile {gtffile} --runThreadN {threads} --sjdbOverhang {length}'
    logging.info(commands)
    status = subprocess.run(commands.split(sep=" "), stderr=subprocess.PIPE, stdout=subprocess.PIPE)

if __name__ == '__main__':
    build_star_index()
