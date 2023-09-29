#!/usr/bin/env python

import argparse
import os

from pathlib import Path

from jaqpotpy.models import MolecularModel
from dm_job_utilities.dm_log import DmLog

import rdkit_utils


models_meta = {
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

# memo to self: model execution times (from a non-scientific superficial test):
# oZZfU6RQgLnmHgk88hnc - 11s
# y7ymAVUixvLs6tBmwdjZ -  6s
# the rest: 0.3 - 1.5s

model_file = "{}.jmodel"

# testing locally vs. running in container
# model_path = Path(".").joinpath("models").absolute()
model_path = Path(__file__).parent.joinpath("models").absolute()



def run(
    model_ids: list,
    input_filename: str,
    output_filename: str,
    delimiter: str = "\t",
    read_header: bool = True,
    write_header: bool = True,
    id_column=None,
    sdf_read_records: int = 100,
    reporting_interval: int = 100,
):

    # TODO: when there's more models, reading them in advance may put
    # too much pressure on memory. it's not too bad now, but may need
    # to be evaluated later
    models = {}
    for model_id in set(model_ids):
        try:
            models[model_id] = MolecularModel().load(
                str(model_path.joinpath(model_file.format(model_id)))
            )
        except FileNotFoundError:
            DmLog.emit_event(f"Model {model_id} not found!")
            continue
        
    reader = rdkit_utils.create_reader(
        input_filename,
        delimiter=delimiter,
        read_header=read_header,
        id_column=id_column,
        sdf_read_records=sdf_read_records,
    )

    extra_field_names = reader.get_extra_field_names()
    
    writer = rdkit_utils.create_writer(
        output_filename,
        delimiter=delimiter,
    )

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

        num_outputs += 1
        if (count + 1) % reporting_interval == 0:
            DmLog.emit_event(f'{count + 1} molecules processed')

        values = []
        calc_prop_names = []
                    
        for model_id, model in models.items():
            # actual prediction
            model(mol)
            values.extend(get_calc_values(model))
            calc_prop_names.extend(get_calc_prop_names(model, format_name(models_meta[model_id])))

            # impractical to have them here
            # model_type = "classification" if model.probability else "regression"
            # DmLog.emit_event(f'Running "{models_meta[model_id]}" ({model_type})')


        if count == 1 and write_header:
            headers = rdkit_utils.generate_header_values(extra_field_names, len(props), calc_prop_names)
            writer.write_header(headers)

            
    
        writer.write(
            smiles=smi,
            mol=mol,
            mol_id=mol_id,
            existing_props=props,  # only used in SmilesWriter
            prop_names=calc_prop_names,  # only used in SdfWriter
            new_props=values,
        )

    reader.close()
    writer.close()
    os.chmod(output_filename, 0o664)
    
    DmLog.emit_event(num_outputs, "outputs among", count, "molecules")
    DmLog.emit_cost(count * len(models.keys()))


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

def format_name(name):
    """Return sd-file friendly name"""
    return name.replace(" ", "_").replace("__", "_")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict properties with a Jaqpot model")
    parser.add_argument("models", metavar="Model ID", nargs="+", help="List of Jaqpot model IDs")
    parser.add_argument("-i", "--input", required=True, help="Input")
    parser.add_argument("-o", "--output", default="result.sdf", help="The output file")
    parser.add_argument("-d", "--delimiter", default="\t", help="Delimiter when using SMILES")
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
    parser.add_argument(
        "--reporting-interval",
        default=100,
        type=int,
        help="Log progress messages after N records",
    )    

    args = parser.parse_args()

    run(
        args.models,
        args.input,
        args.output,
        delimiter=args.delimiter,
        read_header=args.read_header,
        id_column=args.id_column,
        sdf_read_records=args.sdf_read_records,
        reporting_interval=args.reporting_interval,
    )
