import numpy as np

from typing import Any, Dict

import pylab as plt  # type: ignore

import rombus.exceptions as exceptions
from rombus.rom import ReducedOrderModel
from rombus._core.log import log

import warnings

warnings.simplefilter("ignore", np.ComplexWarning)


def bases_errors(ROM: ReducedOrderModel) -> None:
    """Generate plot of errors.

    Parameters
    ----------
    err_list : List[float]
        List of errors
    """

    filename_out = f"{ROM.basename}_bases_errors.pdf"

    with log.context("Generating plot of errors"):
        if ROM is None or ROM.reduced_basis is None:
            raise exceptions.RombusPlotError(
                "Basis plots can not be generated for uninitialised ROM."
            )
        err_list = ROM.reduced_basis.error_list
        plt.plot(err_list)
        plt.xlabel("# Basis elements")
        plt.ylabel("Error")
        plt.yscale("log")
        plt.tight_layout()
        plt.savefig(filename_out)
        log.append(f"written to {filename_out}")


def bases(ROM: ReducedOrderModel) -> None:
    """Generate a plot of the Reduced Bases.

    Parameters
    ----------
    rb_matrix : List[np.ndarray]
        List of reduced bases
    """
    filename_out = f"{ROM.basename}_bases.pdf"

    with log.context("Generating plot of bases"):
        if ROM is None or ROM.reduced_basis is None:
            raise exceptions.RombusPlotError(
                "Basis plots can not be generated for uninitialised ROM."
            )
        rb_matrix = ROM.reduced_basis.matrix
        num_elements = len(rb_matrix)
        total_frames = 125
        n_models_per_frame = int(num_elements / total_frames)
        if n_models_per_frame < 1:
            n_models_per_frame = 1
        fig, ax = plt.subplots(total_frames, 1, figsize=(4.5, 2.5 * total_frames))
        for i in range(total_frames):
            start_i = int(i * n_models_per_frame)
            end_i = int(start_i + n_models_per_frame)
            for model_id in range(start_i, end_i):
                if end_i < num_elements:
                    h = rb_matrix[model_id]
                    ax[i].plot(h, color=f"C{model_id}", alpha=0.7)
            ax[i].set_title(f"Basis element {start_i:003}-{end_i:003}")
        plt.tight_layout()
        fig.savefig(filename_out)
        log.append(f"written to {filename_out}")


def compare_rom_to_true(
    ROM: ReducedOrderModel, model_params_in: Dict[str, Any]
) -> None:
    """Generate a plot comparing the ROM to the original source model for a given set of parameters.

    Parameters
    ----------
    ROM : ReducedOrderModel
        The Reduced Order Model to make the comparison for
    model_params_in : Dict[str, Any]
        A dictionary of parameters to use as the input parameters
    """

    filename_out = f"{ROM.basename}_comparison.pdf"
    filename_diff_out = f"{ROM.basename}_ROM_diff"

    with log.context("Generating comparison plot"):
        if ROM.model is None or ROM.empirical_interpolant is None:
            raise exceptions.RomNotInitialised(
                "ROM not initialised when generating comparison plot."
            )
        else:

            # N.B.: mypy struggles with NamedTuples, so typing is turned off for the following
            model_params = ROM.model.sample(model_params_in)
            # model_params = ROM.model.params.params_dtype(**model_params_in)  # type: ignore

            domain = ROM.model.domain

            model_full = ROM.model.compute(model_params, domain)

            if ROM.empirical_interpolant.nodes is None:
                raise exceptions.RomNotInitialised(
                    "ROM's EmpiricalInterpolant is uninitialised when generating comparison plot"
                )
            else:
                model_nodes = ROM.model.compute(
                    model_params, ROM.empirical_interpolant.nodes
                )
                model_rom = ROM.evaluate(model_params)

                plt.xlabel(ROM.model.coordinate.label)
                plt.ylabel(ROM.model.ordinate.label)
                plt.semilogx(domain, model_rom, label="ROM", alpha=0.5, linestyle="--")
                plt.semilogx(domain, model_full, label="Full model", alpha=0.5)
                plt.scatter(ROM.empirical_interpolant.nodes, model_nodes, s=1)
                plt.legend()

                np.save(filename_diff_out, np.subtract(model_rom, model_full))
                plt.savefig(filename_out, bbox_inches="tight")
                log.append(f"written to {filename_out}")
