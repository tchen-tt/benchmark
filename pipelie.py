import click
import logging
from typing import Any
from collections import OrderedDict
import sys

from index import build_star_index
from alignmentDevo import star_alignment_denovo
from benchmark import benchmark



@click.group()
def cli() -> None:
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    return

cli.add_command(build_star_index)
cli.add_command(star_alignment_denovo)
cli.add_command(benchmark)

if __name__ == '__main__':
    cli()
