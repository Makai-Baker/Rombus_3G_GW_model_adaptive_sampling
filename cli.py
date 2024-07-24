import sys

import click

import importlib.metadata
import collections

import rombus.plots as plots

from rombus.model import RombusModel
from rombus.samples import Samples
from rombus.rom import ReducedOrderModel
from rombus._core.log import log
from typing import Tuple


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
FLEX_CONTEXT_SETTINGS = dict(
    help_option_names=["-h", "--help"],
    ignore_unknown_options=True,
    allow_extra_args=True,  # needed for the passing of model parameters
)


class _OrderedGroup(click.Group):
    """This class is used to ensure that the ordering of the CLI subcommands are in code-order in the CLI help and documentation."""

    def __init__(self, name=None, commands=None, **attrs):
        super(_OrderedGroup, self).__init__(name, commands, **attrs)
        #: the registered subcommands by their exported names.
        self.commands = commands or collections.OrderedDict()
        """An ordered dict of the registered commands"""

    def list_commands(self, ctx: click.core.Context):
        return self.commands


@click.group(cls=_OrderedGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(importlib.metadata.version("rombus"))
@click.pass_context
def cli(ctx: click.core.Context) -> None:
    """Build and work with Reduced Order Models (ROMs) with Rombus.

    Please see the Rombus documentation at https://rombus.readthedocs.io for information about how to work with Rombus.

    N.B.: The commands listed below are roughly in the order that most users will need to use them.
    """

    # ensure that ctx.obj exists and is a dict (in case `cli()` is called
    # by means other than the `if` block below)
    ctx.ensure_object(dict)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("project_name", type=str)
def quickstart(project_name: str) -> None:
    """Write a project template to kickstart the creation of a new project."""

    with log.context("Creating new project from template", time_elapsed=True):
        RombusModel.write_project_template(project_name)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("model", type=str)
@click.argument("filename_samples", type=click.Path(exists=True))
@click.option(
    "--out",
    "-o",
    default="MODEL_BASENAME.hdf5",
    type=click.Path(exists=False),
    show_default=True,
    help="Output HDF5 filename",
)
@click.option(
    "--do_step",
    type=click.Choice(["RB", "EI"], case_sensitive=False),
    default=None,
    help="Do only one step: RB=reduced basis; EI=empirical interpolant",
)
@click.pass_context
def build(
    ctx: click.core.Context,
    model: str,
    filename_samples: str,
    out: str,
    do_step: click.Choice,
) -> None:
    """Build a Reduced Order Model."""

    with log.context("Building ROM", time_elapsed=True):
        log.comment(f"Model:   {model}")
        log.comment(f"Samples: {filename_samples}")

        # Load model
        model_loaded = RombusModel.load(model)

        # Load samples
        samples = Samples(model_loaded, filename=filename_samples)

        # Build ROM
        assert do_step is None or type(do_step) == str  # keeping mypy happy
        ROM = ReducedOrderModel(model_loaded, samples).build(do_step=do_step)

        # Write ROM
        if out == "MODEL_BASENAME.hdf5":
            filename_out = f"{model_loaded.basename}.hdf5"
        else:
            filename_out = out
        ROM.write(filename_out)


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("filename_ROM", type=click.Path(exists=True))
@click.pass_context
def refine(ctx: click.core.Context, filename_rom: str) -> None:
    """Refine parameter sampling to improve an established reduced order model."""

    with log.context(f"Refining ROM ({filename_rom})", time_elapsed=True):

        # Build model and refine it
        ROM = ReducedOrderModel.from_file(filename_rom).refine()

        # Write results
        filename_out = f"{ROM.basename}_refined.hdf5"
        ROM.write(filename_out)


@cli.command(context_settings=FLEX_CONTEXT_SETTINGS)
@click.argument("filename_ROM", type=str)
@click.pass_context
def plot_bases(
    ctx: click.core.Context,
    filename_rom: str,
) -> None:
    """Plot the bases (and their errors) of a reduced order model."""

    with log.context("Plotting bases", time_elapsed=True):

        # Read ROM
        ROM = ReducedOrderModel.from_file(filename_rom)

        # Generate plots
        plots.bases(ROM)
        plots.bases_errors(ROM)


@cli.command(context_settings=FLEX_CONTEXT_SETTINGS)
@click.argument("filename_ROM", type=str)
@click.argument("parameters", nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def evaluate(
    ctx: click.core.Context, filename_rom: str, parameters: Tuple[str, ...]
) -> None:
    """Evaluate a reduced order model and compare it to truth.

    PARAMETERS is a list of parameter values of the form A=VAL B=VAL ...
    """

    with log.context("Evaluating ROM", time_elapsed=True):

        # Read ROM
        ROM = ReducedOrderModel.from_file(filename_rom)

        # Parse the model parameters, which should have been given as arguments
        model_params = ROM.model.parse_cli_params(parameters)

        # Generate plot
        plots.compare_rom_to_true(ROM, model_params)


@cli.command(context_settings=FLEX_CONTEXT_SETTINGS)
@click.argument("filename_ROM", type=str)
@click.option(
    "-n",
    "--n_samples",
    type=int,
    default=100,
    help="Number of samples to use for timing",
)
@click.pass_context
def timing(ctx: click.core.Context, filename_rom: str, n_samples: int) -> None:
    """Compute timing information for a ROM and it's source model."""

    with log.context(
        f"Creating timing information for ROM {filename_rom}", time_elapsed=True
    ):
        # Read ROM
        ROM = ReducedOrderModel.from_file(filename_rom)

        # Generate the samples to be used
        timing_sample = Samples(ROM.model, n_random=n_samples)

        # Generate timing information for model
        timing_model = ROM.model.timing(timing_sample)

        # Generate timing information for ROM
        timing_ROM = ROM.timing(timing_sample)

    # Report results
    log.comment(
        f"Timing information for ROM:   {timing_ROM:.2e}s for {n_samples} calls ({timing_ROM/n_samples:.2e} per sample)."
    )
    log.comment(
        f"Timing information for model: {timing_model:.2e}s for {n_samples} calls ({timing_model/n_samples:.2e} per sample)."
    )
    if timing_ROM > timing_model:
        log.comment(
            f"ROM is {timing_ROM/timing_model:.2f}X slower than the source model."
        )
    else:
        log.comment(
            f"ROM is {timing_model/timing_ROM:.2f}X faster than the source model."
        )


if __name__ == "__main__":
    cli(obj={})
    sys.exit()
