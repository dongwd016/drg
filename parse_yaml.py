"""Contains main driver function for pyMARS program."""
import os
import sys
import logging
from argparse import ArgumentParser
from typing import List, Dict, NamedTuple

import yaml
import cantera as ct

# local imports
from sampling import sample_metrics, parse_ignition_inputs, parse_psr_inputs, parse_flame_inputs, parse_znd_inputs
from sampling import InputIgnition, InputPSR, InputLaminarFlame, InputZND

import sampling
import soln2cti
from drgep import run_drgep
from drg import run_drg
from tools import convert

#: Supported reduction methods
METHODS = ['DRG', 'DRGEP', 'PFA']

class ReductionInputs(NamedTuple):
    """Collects inputs for overall reduction process.
    """
    model: str
    error: float
    ignition_conditions: List[InputIgnition]
    psr_conditions: List[InputPSR]
    flame_conditions: List[InputLaminarFlame]
    znd_conditions: List[InputZND]
    method: str
    target_species: List[str]
    safe_species: List[str] = []
    sensitivity_analysis: bool = False
    upper_threshold: float = 0.1
    sensitivity_type: str = 'greedy'
    phase_name: str = ''


def parse_inputs(input_dict):
    """Parses and checks dictionary of inputs for consistency and correctness.

    Parameters
    ----------
    input_dict : dict
        Inputs for reduction

    Returns
    -------
    ReductionInputs
        Object with checked inputs
    
    """
    model = input_dict.get('model', '')
    assert model, 'Input file requires specifying "model".'
    
    error = input_dict.get('error', 0.0)
    assert error, 'Input file requires an error limit specified by "error".'

    method = input_dict.get('method', '')
    sensitivity_analysis = input_dict.get('sensitivity-analysis', False)
    
    assert method or sensitivity_analysis, (
        'Input file requires either "method" or "sensitivity-analysis" to be given.'
        )
    
    if method:
        assert method in METHODS, 'Reduction method must be one of ' + ', '.join(METHODS)

        target_species = input_dict.get('targets', [])
        assert target_species, (
            'At least one "target" species must be specified for graph-based reduction methods.'
            )
    
    upper_threshold = input_dict.get('upper-threshold', 0.1)
    sensitivity_type = input_dict.get('sensitivity-type', 'initial')

    safe_species = input_dict.get('retained-species', [])

    phase_name = input_dict.get('phase-name', '')
    
    # check that the specified model actually contains the specified phase
    try:
        gas = ct.Solution(model, phase_name)
    except ValueError:
        raise ValueError(model + ' does not contain phase ' + phase_name)

    # check that species are present in model
    for sp in target_species:
        assert sp in gas.species_names, f'Specified target species {sp} not in model'
    
    for sp in safe_species:
        assert sp in gas.species_names, f'Specified retained species {sp} not in model'
    
    ignition_conditions = input_dict.get('autoignition-conditions', {})
    psr_conditions = input_dict.get('psr-conditions', {})
    flame_conditions = input_dict.get('laminar-flame-conditions', {})
    znd_conditions = input_dict.get('znd-conditions', {})
    
    # check validity of input file
    ignition_inputs = parse_ignition_inputs(model, ignition_conditions, phase_name)
    psr_inputs = parse_psr_inputs(model, psr_conditions, phase_name)
    flame_inputs = parse_flame_inputs(model, flame_conditions, phase_name)
    znd_inputs = parse_znd_inputs(model, znd_conditions, phase_name)

    return ReductionInputs(
        model=model, error=error, 
        ignition_conditions=ignition_inputs, 
        psr_conditions=psr_inputs, flame_conditions=flame_inputs, znd_conditions=znd_inputs,
        method=method, target_species=target_species, safe_species=safe_species,
        sensitivity_analysis=sensitivity_analysis, upper_threshold=upper_threshold,
        sensitivity_type=sensitivity_type, phase_name=phase_name
        )