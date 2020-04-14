#!/usr/bin/env python
import click
from cream.data_scripts import data


@click.group()
def main():
    """ Command line utility for SCONE Monte Carlo Code.
    """
    pass


# Add nested command groups to main
main.add_command(data)
