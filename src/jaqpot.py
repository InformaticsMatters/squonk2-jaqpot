#!/usr/bin/env python

import argparse
import os

from pathlib import Path

from jaqpotpy.models import MolecularModel
from dm_job_utilities.dm_log import DmLog

from urllib.request import urlretrieve
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.parse import urljoin

import rdkit_utils


models_meta = {
    "solubility": "Aqueous solubility model",
    "herg": "hERG model",
    "AMES": "AMES model",
    "CYP2C9_Veith": "CYP2C9 inhibition model",
    "CYP3A4_Veith": "CYP3A4 inhibition model",
    "CYP2C9_Substrate_CarbonMangels": "CYP2C9 substrate model",
    "CYP2D6_Veith": "CYP2D6 inhibition model",
    "lipophilicity": "Lipophilicity model",
    "ppbr_az": "PPBR model",
    "hia_hou": "HIA model",
    "CYP2D6_Substrate_CarbonMangels": "CYP2D6 substrate model",
    "bioavailability_ma": "Bioavailability model",
    "clearance_microsome_az": "Clearance Microsome model",
    "ld50_zhu": "LD50 model",
    "CYP3A4_Substrate_CarbonMangels": "CYP3A4 Substrate CarbonMangels  model",
    "caco2_wang": "CaCO2 Wang model",
    "dili": "DILI model",
    "vdss_lombardo": "VDss Lombardo model",
    "clearance_hepatocyte_az": "Clearance Hepatocyte model",
    "half_life+obach": "Half Life Obach model",
    "BBB": "Blood Brain Barrier model",
    "pgp": "PGP model",
 }




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
    model_base_path: str = "",
):

    # temporary measure
    model_base_path = "https://im-jaqpot-models.s3.eu-central-1.amazonaws.com"
    # if not given, try to extract from env variable
    if not model_base_path:
        try:
            model_base_path = os.environ['BASE_MODEL_URL']
        except KeyError:
            DmLog.emit_event(f"Base model url not set!")
            # assume testing regime and use the default directory
            model_base_path = "/models"

    # TODO: when there's more models, reading them in advance may put
    # too much pressure on memory. it's not too bad now, but may need
    # to be evaluated later
    models = {}
    for model_id in set(model_ids):
        if urlparse(model_base_path).netloc:
            try:
                model_file = urlretrieve(urljoin(model_base_path, f"{model_id}.jmodel"))[0]
            except HTTPError:
                DmLog.emit_event(f"Model {model_id} not available!")
                continue
        else:
            model_file = Path(model_base_path).joinpath(f"{model_id}.jmodel")


        try:
            models[model_id] = MolecularModel().load(model_file)
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
    parser.add_argument(
        "--model-base-path",
        default="",
        type=str,
        help="Model location, URL or path",
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
        model_base_path=args.model_base_path,
    )
