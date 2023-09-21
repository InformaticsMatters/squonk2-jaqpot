#!/usr/bin/env python

import argparse
import os

from pathlib import Path

from jaqpotpy.models import MolecularModel
from dm_job_utilities.dm_log import DmLog

import rdkit_utils


models = {
    "fUAo2UQO8tTGZFhd5fPB": "Aqueous solubility model",
    "AdIueWr1VDrWC3j90jjX": "hERG model",
    "88NHffXLTX3aBM2vkmOf": "AMES model",
    "7JnhJUBH1wxwB7Vgf8YI": "CYP2C9 inhibition model",
    "4UfyxtoMWFuhN2PBrK42": "CYP3A4 inhibition model",
    "HI7FUfl5phSpxGSYjes1": "CYP2C9 substrate model",
    "dud9GNQZaBZ9grt7VMMA": "CYP2D6 inhibition model",
    "tRgpmWmuBImTw3gC8NXE": "Lipophilicity model",
    "FR60WJT6qZoTE1L7tNzx": "PPBR model",
    "Z7OzhVDtxaTyMLscRJ4v": "HIA model",
    "cp4HGKxIxjAsdiM5T6Oj": "CYP2D6 substrate model",
    "Em70hoXbIqcTvqscjDFu": "Bioavailability model",
    "GfgMcorzIljrChpL824z": "Clearance Microsome model",
    "gpA3uw3FM8TWbAZF2BWP": "LD50 model",
    "Mc8Uc3J1tvEnkkqlsiAE": "CYP3A4 Substrate CarbonMangels  model",
    "2YSbIaKrHFTf1MnyswIk": "CaCO2 Wang model",
    "nL7dOTzvaoGRmp1wiv8m": "DILI model",
    "MLQIb8KSFdLGSIQeUFaV": "VDss Lombardo model",
    "6oO0hMJyz4s0Sz62OLo6": "Clearance Hepatocyte model",
    "y7ymAVUixvLs6tBmwdjZ": "Half Life Obach model",
    "oZZfU6RQgLnmHgk88hnc": "Blood Brain Barrier model",
    "llKNcGM5vuGhf6EGkOpA": "PGP model",
}


model_file = "{}.jmodel"

# testing locally vs. running in container
# model_path = Path(".").joinpath("models").absolute()
model_path = Path(__file__).parent.joinpath("models").absolute()


def run(
    model_id: str,
    input_filename: str,
    output_filename: str,
    delimiter: str = "\t",
    read_header: bool = True,
    id_column=None,
    sdf_read_records: int = 100,
):
    try:
        model = MolecularModel().load(
            str(model_path.joinpath(model_file.format(model_id)))
        )
    except FileNotFoundError:
        DmLog.emit_event(f"Model {model_id} not found!")
        return

    model_title = models[model_id]
    model_name = model_title.replace(" ", "_").replace("__", "_")

    reader = rdkit_utils.create_reader(
        input_filename,
        delimiter=delimiter,
        read_header=read_header,
        id_column=id_column,
        sdf_read_records=sdf_read_records,
    )

    extra_field_names = reader.get_extra_field_names()

    num_outputs = 0
    count = -1
    while True:
        count += 1
        try:
            mol, smi, mol_id, props = reader.read()
        except TypeError as ex:
            DmLog.emit_event(f"{ex}")
            continue
        except StopIteration as ex:
            # end of file
            break

        # actual prediction
        model(mol)
        values = get_calc_values(model)

        num_outputs += 1

        try:
            writer.write(smi, mol, mol_id, props, values)
        except NameError:
            # writer not defined yet, do so now when I have the first
            # results
            calc_prop_names = get_calc_prop_names(model, model_name)
            writer = rdkit_utils.create_writer(
                output_filename,
                extra_field_names=extra_field_names,
                calc_prop_names=calc_prop_names,
                delimiter=delimiter,
            )
            model_type = "classification" if model.probability else "regression"
            DmLog.emit_event(
                "Model title:",
                model_title,
                ", property name",
                model_name,
                ", type:",
                model_type,
            )
            writer.write(smi, mol, mol_id, props, values)

    os.chmod(output_filename, 0o664)

    DmLog.emit_event(num_outputs, "outputs among", count, "molecules")
    DmLog.emit_cost(count)


def get_calc_prop_names(molmod, prefix):
    """
    Get the names of the properties that will be output.
    These will be used as the field names of the values obtained from the get_calc_values method.
    :param molmod: The Jaqpot model to get the data from
    :param prefix: The prefix for the field names
    :return: List of names
    """
    names = [prefix + "_Prediction"]
    try:
        _ = molmod.probability[0][0]
        # has attribute, means classification model
        names.append(prefix + "_Inactive")  # probability of being inactive
        names.append(prefix + "_Active")  # probability of being active
    except IndexError:
        # regression rather than classification
        pass

    try:
        _ = molmod.doa.IN[0]
        names.append(prefix + "_DOA")  # True or False
    except AttributeError:
        pass

    return names


def get_calc_values(molmod):
    """
    Get the values that should be output.
    This depends on the type of the model.
    :param molmod: The Jaqpot model
    :return: List of values
    """
    values = [molmod.prediction[0]]
    try:
        values.append(molmod.probability[0][0])
        values.append(molmod.probability[0][1])
    except IndexError:
        pass

    try:
        values.append(molmod.doa.IN[0])
    except AttributeError:
        pass

    return values


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calqulate aqueous solubility")
    parser.add_argument("-m", "--model", required=True, help="Jaqpot model name")
    parser.add_argument("-i", "--input", required=True, help="Input")
    parser.add_argument("-o", "--output", required=True, help="The output file")
    parser.add_argument("-d", "--delimiter", help="Delimiter when using SMILES")
    parser.add_argument(
        "--id-column",
        help="Column for name field (zero based integer for .smi, text for SDF)",
    )
    parser.add_argument(
        "--read-header",
        action="store_true",
        help="Read a header line with the field names when reading .smi or .txt",
    )
    parser.add_argument(
        "--write-header",
        action="store_true",
        help="Write a header line when writing .smi or .txt",
    )
    parser.add_argument(
        "--sdf-read-records",
        default=100,
        type=int,
        help="Read this many SDF records to determine field names",
    )

    args = parser.parse_args()

    run(
        args.model,
        args.input,
        args.output,
        delimiter=args.delimiter,
        read_header=args.read_header,
        id_column=args.id_column,
        sdf_read_records=args.sdf_read_records,
    )
