"""Main CLI entrypoint for Lotus tooling."""

from __future__ import annotations

import click

from .commands.inject import inject_cmd


@click.group()
def cli() -> None:
    """Lotus command group."""


cli.add_command(inject_cmd, name="inject")


if __name__ == "__main__":
    cli()
