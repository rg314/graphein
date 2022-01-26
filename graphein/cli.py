"""Command line interface for Graphein."""
# Graphein
# Author: Arian Jamasb <arian@jamasb.io>
# License: MIT
# Code Repository: https://github.com/a-r-j/graphein
import click

from graphein import __version__
from graphein.utils.utils import parse_config  # , run_from_config


@click.command()
@click.version_option(__version__)
@click.argument(
    "config", type=click.Path(exists=True, file_okay=True, dir_okay=False)
)
def main(path):
    config = parse_config(path)
    # run_from_config(config)


if __name__ == "__main__":
    main()
