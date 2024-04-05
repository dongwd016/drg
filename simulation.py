"""Autoignition simulation module

.. moduleauthor:: Kyle Niemeyer <kyle.niemeyer@gmail.com>
"""

# Standard libraries
import os
import logging
import json

# Related modules
import numpy as np
import tables
import cantera as ct
from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.znd import zndsolve
from sdtoolbox.utilities import CJspeed_plot, znd_plot, znd_fileout

ct.suppress_thermo_warnings()


def json_convert(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if hasattr(obj, "tolist"):
        return {"$nparray": obj.tolist()}


def json_deconvert(obj):
    if len(obj) == 1:
        key, value = next(iter(obj.items()))
        if key == "$nparray":
            return np.array(value)
    return obj


# class Simulation_ign(object):
#     """Class for ignition delay simulations

#     Parameters
#     ----------
#     idx : int
#         Identifer index for case
#     properties : InputIgnition
#         Object with initial conditions for simulation
#     model : str
#         Filename for Cantera-format model to be used
#     phase_name : str, optional
#         Optional name for phase to load from CTI file (e.g., 'gas').
#     path : str, optional
#         Path for location of output files

#     """
#     def __init__(self, idx, properties, model, phase_name='', path=''):
#         self.idx = idx
#         self.properties = properties
#         self.model = model
#         self.phase_name = phase_name
#         self.path = path

#     def setup_case(self):
#         """Initialize simulation case.
#         """
#         self.gas = ct.Solution(self.model, self.phase_name)

#         # Default maximum number of steps
#         self.max_steps = 10000
#         if self.properties.max_steps:
#             self.max_steps = self.properties.max_steps

#         # By default, simulations will run to steady state, with the maximum number of steps
#         # given by ``self.max_steps``. Alternatively, an end time (in seconds) can be
#         # given in cases where something specific is needed (e.g., longer than normal)
#         self.time_end = 0.0
#         if self.properties.end_time:
#             self.time_end = self.properties.end_time

#         self.gas.TP = (
#             self.properties.temperature, self.properties.pressure * ct.one_atm
#             )
#         # set initial composition using either equivalence ratio or general reactant composition
#         if self.properties.equivalence_ratio:
#             self.gas.set_equivalence_ratio(
#                 self.properties.equivalence_ratio,
#                 self.properties.fuel,
#                 self.properties.oxidizer
#                 )
#         else:
#             if self.properties.composition_type == 'mole':
#                 self.gas.TPX = (
#                     self.properties.temperature, self.properties.pressure * ct.one_atm,
#                     self.properties.reactants
#                     )
#             else:
#                 self.gas.TPY = (
#                     self.properties.temperature, self.properties.pressure * ct.one_atm,
#                     self.properties.reactants
#                     )

#         if self.properties.kind == 'constant pressure':
#             self.reac = ct.IdealGasConstPressureReactor(self.gas)
#         else:
#             self.reac = ct.IdealGasReactor(self.gas)

#         # Create ``ReactorNet`` newtork
#         self.sim = ct.ReactorNet([self.reac])

#         # Set file for later data file
#         self.save_file = os.path.join(self.path, str(self.idx) + '.h5')
#         self.sample_points = []

#         self.ignition_delay = 0.0

#     def run_case(self, stop_at_ignition=False, restart=False):
#         """Run simulation case set up ``setup_case``.

#         If no end time is specified for the integration, the function integrates
#         to steady state (or a maximum of 10,000 steps, by default). This is done
#         by checking whether the system state changes below a certain threshold,
#         with the residual computed using feature checking. This is blatantly stolen
#         from Cantera's :meth:`cantera.ReactorNet.advance_to_steady_state` method.

#         Parameters
#         ----------
#         stop_at_ignition : bool
#             If ``True``, stop integration at ignition point, don't save data.
#         restart : bool
#             If ``True``, skip if results file exists.

#         Returns
#         -------
#         self.ignition_delay : float
#             Computed ignition delay in seconds

#         """

#         if restart and os.path.isfile(self.save_file):
#             print('Skipped existing case ', self.idx)
#             return

#         # Save simulation results in hdf5 table format.
#         table_def = {'time': tables.Float64Col(pos=0),
#                      'temperature': tables.Float64Col(pos=1),
#                      'pressure': tables.Float64Col(pos=2),
#                      'mass_fractions': tables.Float64Col(
#                           shape=(self.reac.thermo.n_species), pos=3
#                           ),
#                      }

#         with tables.open_file(self.save_file, mode='w',
#                               title=str(self.idx)
#                               ) as h5file:

#             table = h5file.create_table(where=h5file.root,
#                                         name='simulation',
#                                         description=table_def
#                                         )
#             # Row instance to save timestep information to
#             timestep = table.row
#             # Save initial conditions
#             timestep['time'] = self.sim.time
#             timestep['temperature'] = self.reac.T
#             timestep['pressure'] = self.reac.thermo.P
#             timestep['mass_fractions'] = self.reac.Y
#             # Add ``timestep`` to table
#             timestep.append()

#             ignition_flag = False

#             # Main time integration loop
#             if self.time_end:
#                 # if end time specified, continue integration until reaching that time
#                 while self.sim.time < self.time_end:
#                     self.sim.step()

#                     # Save new timestep information
#                     timestep['time'] = self.sim.time
#                     timestep['temperature'] = self.reac.T
#                     timestep['pressure'] = self.reac.thermo.P
#                     timestep['mass_fractions'] = self.reac.Y

#                     if self.reac.T >= self.properties.temperature + 400.0 and not ignition_flag:
#                         self.ignition_delay = self.sim.time
#                         ignition_flag = True

#                         if stop_at_ignition:
#                             break

#                     # Add ``timestep`` to table
#                     timestep.append()

#             else:
#                 # otherwise, integrate until steady state, or maximum number of steps reached
#                 self.sim.reinitialize()
#                 max_state_values = self.sim.get_state()
#                 residual_threshold = 10. * self.sim.rtol
#                 absolute_tolerance = self.sim.atol

#                 for step in range(self.max_steps):
#                     previous_state = self.sim.get_state()

#                     self.sim.step()

#                     # Save new timestep information
#                     timestep['time'] = self.sim.time
#                     timestep['temperature'] = self.reac.T
#                     timestep['pressure'] = self.reac.thermo.P
#                     timestep['mass_fractions'] = self.reac.Y

#                     if self.reac.T >= self.properties.temperature + 400.0 and not ignition_flag:
#                         self.ignition_delay = self.sim.time
#                         ignition_flag = True

#                         if stop_at_ignition:
#                             break

#                     # Add ``timestep`` to table
#                     timestep.append()

#                     state = self.sim.get_state()
#                     max_state_values = np.maximum(max_state_values, state)
#                     residual = np.linalg.norm(
#                         (state - previous_state) / (max_state_values + absolute_tolerance)
#                         ) / np.sqrt(self.sim.n_vars)

#                     if residual < residual_threshold:
#                         break

#                 if step == self.max_steps - 1:
#                     logging.error(
#                         'Maximum number of steps reached before '
#                         f'convergence for ignition case {self.idx}'
#                         )
#                     raise RuntimeError(
#                         'Maximum number of steps reached before '
#                         f'convergence for ignition case {self.idx}'
#                     )

#             # Write ``table`` to disk
#             table.flush()

#             if not ignition_flag:
#                 logging.error(f'No ignition detected for ignition case {self.idx}')
#                 raise RuntimeError(f'No ignition detected for ignition case {self.idx}')

#         return self.ignition_delay

#     def calculate_ignition(self):
#         """Run simulation case set up ``setup_case``, just for ignition delay.
#         """
#         # Main time integration loop
#         if self.time_end:
#             # if end time specified, continue integration until reaching that time
#             while self.sim.time < self.time_end:
#                 self.sim.step()
#                 if self.reac.T >= self.properties.temperature + 400.0:
#                     self.ignition_delay = self.sim.time
#                     break
#             if not self.ignition_delay:
#                 logging.warning(
#                     f'No ignition detected before end time for ignition case {self.idx}'
#                     )
#         else:
#             # otherwise, integrate until steady state, or maximum number of steps reached
#             for step in range(self.max_steps):
#                 self.sim.step()
#                 if self.reac.T >= self.properties.temperature + 400.0:
#                     self.ignition_delay = self.sim.time
#                     break
#             if step == self.max_steps - 1:
#                 logging.warning(
#                     'Maximum number of steps reached before '
#                     f'convergence for ignition case {self.idx}'
#                     )

#         return self.ignition_delay

#     def process_results(self, skip_data=False):
#         """Process integration results to sample data

#         Parameters
#         ----------
#         skip_data : bool
#             Flag to skip sampling thermochemical data

#         Returns
#         -------
#         tuple of float, numpy.ndarray or float
#             Ignition delay, or ignition delay and sampled data

#         """
#         delta = 0.05
#         deltas = np.arange(delta, 1 + delta, delta)

#         # Load saved integration results
#         self.save_file = os.path.join(self.path, str(self.idx) + '.h5')
#         with tables.open_file(self.save_file, 'r') as h5file:
#             # Load Table with Group name simulation
#             table = h5file.root.simulation

#             times = table.col('time')
#             temperatures = table.col('temperature')
#             pressures = table.col('pressure')
#             mass_fractions = table.col('mass_fractions')

#         temperature_initial = temperatures[0]
#         temperature_max = temperatures[len(temperatures)-1]
#         temperature_diff = temperature_max - temperature_initial

#         sampled_data = np.zeros((len(deltas), 2 + mass_fractions.shape[1]))

#         # need to add processing to get the 20 data points here
#         self.ignition_delay = 0.0
#         ignition_flag = False
#         idx = 0
#         for time, temp, pres, mass in zip(
#             times, temperatures, pressures, mass_fractions
#             ):
#             if temp >= temperature_initial + 400.0 and not ignition_flag:
#                     self.ignition_delay = time
#                     ignition_flag = True
#                     if skip_data:
#                         return self.ignition_delay

#             if temp >= temperature_initial + (deltas[idx] * temperature_diff):
#                 sampled_data[idx, 0:2] = [temp, pres]
#                 sampled_data[idx, 2:] = mass

#                 idx += 1
#                 if idx == 20:
#                     self.sampled_data = sampled_data
#                     return self.ignition_delay, sampled_data

#     def clean(self):
#         """Delete HDF5 file with full integration data.
#         """
#         try:
#             os.remove(self.save_file)
#         except OSError:
#             pass


class Simulation_fls(object):
    """Class for flame speed simulations

    Parameters
    ----------
    idx : int
        Identifer index for case
    properties : InputLaminarFlame
        Object with initial conditions for simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas').
    path : str, optional
        Path for location of output files

    """

    def __init__(self, idx, properties, model, phase_name="", path=""):
        self.idx = idx
        self.properties = properties
        self.model = model
        self.phase_name = phase_name
        self.path = path
        self.save_file = os.path.join(self.path, str(self.idx) + ".json")

    def setup_case(self):
        """Initialize simulation case."""
        self.gas = ct.Solution(self.model, self.phase_name)
        self.gas.TP = (self.properties.temperature, self.properties.pressure * ct.one_atm)
        # set initial composition using either equivalence ratio or general reactant composition
        if self.properties.equivalence_ratio:
            self.gas.set_equivalence_ratio(self.properties.equivalence_ratio, self.properties.fuel,
                                           self.properties.oxidizer)
        else:
            if self.properties.composition_type == "mole":
                self.gas.TPX = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, self.properties.reactants)
            else:
                self.gas.TPY = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, self.properties.reactants)

        # Create simulation object
        self.sim = ct.FreeFlame(self.gas, width=2e-3)

        self.sample_points = []
        self.flame_speed = 0.0

    def run_case(self, stop_at_ignition=False, restart=False):
        """Run simulation case set up ``setup_case``.

        If no end time is specified for the integration, the function integrates
        to steady state (or a maximum of 10,000 steps, by default). This is done
        by checking whether the system state changes below a certain threshold,
        with the residual computed using feature checking. This is blatantly stolen
        from Cantera's :meth:`cantera.ReactorNet.advance_to_steady_state` method.

        Parameters
        ----------
        stop_at_ignition : bool
            If ``True``, stop integration at ignition point, don't save data.
        restart : bool
            If ``True``, skip if results file exists.

        Returns
        -------
        self.ignition_delay : float
            Computed ignition delay in seconds

        """

        if restart and os.path.isfile(self.save_file):
            print("Skipped existing case ", self.idx)
            return

        self.sim.set_refine_criteria(ratio=10, slope=0.6, curve=0.6)
        self.sim.solve(loglevel=0, auto=True)

        self.sim.set_refine_criteria(ratio=7, slope=0.3, curve=0.3)
        self.sim.solve(loglevel=0, auto=True)

        self.sim.set_refine_criteria(ratio=5, slope=0.1, curve=0.1)
        self.sim.solve(loglevel=0, auto=True)

        # self.sim.set_refine_criteria(ratio=2, slope=0.05, curve=0.05)
        # self.sim.transport_model = 'Multi'
        # self.sim.soret_enabled = True
        # self.sim.solve(loglevel=1, auto=True)

        self.flame_speed = self.sim.velocity[0]
        solution_dict = dict()
        solution_dict["x_locs"] = self.sim.flame.grid
        solution_dict["temperature"] = self.sim.T
        solution_dict["pressure"] = np.ones(len(self.sim.flame.grid)) * self.sim.P
        solution_dict["mass_fractions"] = self.sim.Y.T
        solution_dict["flame_speed"] = self.sim.velocity[0]
        json.dump(solution_dict, open(self.save_file, "w"), default=json_convert)

        return self.flame_speed

    def calculate_flamespeed(self):
        """Run simulation case set up ``setup_case``, just for ignition delay."""
        # Main time integration loop
        self.sim.set_refine_criteria(ratio=10, slope=0.6, curve=0.6)
        self.sim.solve(loglevel=0, auto=True)

        self.sim.set_refine_criteria(ratio=7, slope=0.3, curve=0.3)
        self.sim.solve(loglevel=0, auto=True)

        self.sim.set_refine_criteria(ratio=5, slope=0.1, curve=0.1)
        self.sim.solve(loglevel=0, auto=True)

        # self.sim.set_refine_criteria(ratio=2, slope=0.05, curve=0.05)
        # self.sim.transport_model = 'Multi'
        # self.sim.soret_enabled = True
        # self.sim.solve(loglevel=0, auto=True)

        self.flame_speed = self.sim.velocity[0]

        return self.flame_speed

    def process_results(self, skip_data=False):
        """Process integration results to sample data

        Parameters
        ----------
        skip_data : bool
            Flag to skip sampling thermochemical data

        Returns
        -------
        tuple of float, numpy.ndarray or float
            Ignition delay, or ignition delay and sampled data

        """
        delta = 0.05
        deltas = np.arange(delta, 1 + delta, delta)

        # Load saved integration results
        saved_dict = json.load(open(self.save_file, "r"), object_hook=json_deconvert)

        x_locs = saved_dict["x_locs"]
        temperatures = saved_dict["temperature"]
        pressures = saved_dict["pressure"]
        mass_fractions = saved_dict["mass_fractions"]
        flame_speed = saved_dict["flame_speed"]

        temperature_initial = temperatures[0]
        temperature_max = temperatures[len(temperatures) - 1]
        temperature_diff = temperature_max - temperature_initial

        sampled_data = np.zeros((len(deltas), 2 + mass_fractions.shape[1]))

        idx = 0
        for x_loc, temp, pres, mass in zip(x_locs, temperatures, pressures, mass_fractions):
            if temp >= temperature_initial + (deltas[idx] * temperature_diff):
                sampled_data[idx, 0:2] = [temp, pres]
                sampled_data[idx, 2:] = mass
                idx += 1
                if idx == 20:
                    self.sampled_data = sampled_data
                    return flame_speed, sampled_data

    def clean(self):
        """Delete HDF5 file with full integration data."""
        try:
            os.remove(os.path.join(self.path, str(self.idx) + "_fls.json"))
        except OSError:
            pass


class Simulation_znd(object):
    """Class for flame speed simulations

    Parameters
    ----------
    idx : int
        Identifer index for case
    properties : InputLaminarFlame
        Object with initial conditions for simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas').
    path : str, optional
        Path for location of output files

    """

    def __init__(self, idx, properties, model, phase_name="", path=""):
        self.idx = idx
        self.properties = properties
        self.model = model
        self.phase_name = phase_name
        self.path = path
        self.save_file = os.path.join(self.path, str(self.idx) + "_znd.json")

    def setup_case(self):
        """Initialize simulation case."""
        self.gas_pre = ct.Solution(self.model, self.phase_name)
        self.gas_pre.TP = (self.properties.temperature, self.properties.pressure * ct.one_atm)
        # set initial composition using either equivalence ratio or general reactant composition
        if self.properties.equivalence_ratio:
            self.gas_pre.set_equivalence_ratio(self.properties.equivalence_ratio, self.properties.fuel,
                                               self.properties.oxidizer)
        else:
            if self.properties.composition_type == "mole":
                self.gas_pre.TPX = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, self.properties.reactants)
            else:
                self.gas_pre.TPY = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, self.properties.reactants)

        self.cj_speed, _, _ = CJspeed(self.properties.pressure * ct.one_atm, self.properties.temperature,
                                      self.gas_pre.X, self.model, fullOutput=True)

        self.gas_post = PostShock_fr(self.cj_speed, self.properties.pressure * ct.one_atm, self.properties.temperature,
                                     self.gas_pre.X, self.model)

        self.sample_points = []
        self.induction_length = 0.0

    def run_case(self, stop_at_ignition=False, restart=False):
        """Run simulation case set up ``setup_case``.

        If no end time is specified for the integration, the function integrates
        to steady state (or a maximum of 10,000 steps, by default). This is done
        by checking whether the system state changes below a certain threshold,
        with the residual computed using feature checking. This is blatantly stolen
        from Cantera's :meth:`cantera.ReactorNet.advance_to_steady_state` method.

        Parameters
        ----------
        stop_at_ignition : bool
            If ``True``, stop integration at ignition point, don't save data.
        restart : bool
            If ``True``, skip if results file exists.

        Returns
        -------
        self.ignition_delay : float
            Computed ignition delay in seconds

        """

        if restart and os.path.isfile(self.save_file):
            print("Skipped existing case ", self.idx)
            return

        znd_out = zndsolve(self.gas_post, self.gas_pre, self.cj_speed, t_end=1e-5, advanced_output=True)

        self.induction_length = znd_out["ind_len_ZND"]

        solution_dict = dict()
        solution_dict["temperature"] = znd_out["T"]
        solution_dict["pressure"] = znd_out["P"]
        solution_dict["mass_fractions"] = znd_out["species"]
        solution_dict["induction_length"] = znd_out["ind_len_ZND"]
        json.dump(solution_dict, open(self.save_file, "w"), default=json_convert)

        return self.induction_length

    def calculate_inductionlength(self):
        """Run simulation case set up ``setup_case``, just for ignition delay."""

        znd_out = zndsolve(self.gas_post, self.gas_pre, self.cj_speed, t_end=1e-5, advanced_output=True)
        self.induction_length = znd_out["ind_len_ZND"]
        return self.induction_length

    def process_results(self, skip_data=False):
        """Process integration results to sample data

        Parameters
        ----------
        skip_data : bool
            Flag to skip sampling thermochemical data

        Returns
        -------
        tuple of float, numpy.ndarray or float
            Ignition delay, or ignition delay and sampled data

        """
        delta = 0.05
        deltas = np.arange(delta, 1 + delta, delta)

        # Load saved integration results
        saved_dict = json.load(open(self.save_file, "r"), object_hook=json_deconvert)

        temperatures = saved_dict["temperature"]
        pressures = saved_dict["pressure"]
        mass_fractions = saved_dict["mass_fractions"]
        induction_length = saved_dict["induction_length"]

        temperature_initial = temperatures[0]
        temperature_max = temperatures[len(temperatures) - 1]
        temperature_diff = temperature_max - temperature_initial

        sampled_data = np.zeros((len(deltas), 2 + mass_fractions.shape[1]))

        idx = 0
        for temp, pres, mass in zip(temperatures, pressures, mass_fractions):
            if temp >= temperature_initial + (deltas[idx] * temperature_diff):
                sampled_data[idx, 0:2] = [temp, pres]
                sampled_data[idx, 2:] = mass
                idx += 1
                if idx == 20:
                    # self.sampled_data = sampled_data
                    break
                    # return induction_length, sampled_data
        return induction_length, sampled_data

    def clean(self):
        """Delete HDF5 file with full integration data."""
        try:
            os.remove(os.path.join(self.path, str(self.idx) + "_znd.json"))
        except OSError:
            pass


class Simulation_ign(object):
    """Class for ignition delay simulations

    Parameters
    ----------
    idx : int
        Identifer index for case
    properties : InputIgnition
        Object with initial conditions for simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas').
    path : str, optional
        Path for location of output files

    """

    def __init__(self, idx, properties, model, phase_name="", path=""):
        self.idx = idx
        self.properties = properties
        self.model = model
        self.phase_name = phase_name
        self.path = path
        self.save_file = os.path.join(self.path, str(self.idx) + ".h5")

    def setup_case(self):
        """Initialize simulation case."""
        self.gas = ct.Solution(self.model, self.phase_name)

        # Default maximum number of steps
        self.max_steps = 10000
        if self.properties.max_steps:
            self.max_steps = self.properties.max_steps

        # By default, simulations will run to steady state, with the maximum number of steps
        # given by ``self.max_steps``. Alternatively, an end time (in seconds) can be
        # given in cases where something specific is needed (e.g., longer than normal)
        self.time_end = 0.0
        if self.properties.end_time:
            self.time_end = self.properties.end_time

        self.gas.TP = (self.properties.temperature, self.properties.pressure * ct.one_atm)
        # set initial composition using either equivalence ratio or general reactant composition
        if self.properties.equivalence_ratio:
            self.gas.set_equivalence_ratio(self.properties.equivalence_ratio, self.properties.fuel,
                                           self.properties.oxidizer)
        else:
            if self.properties.composition_type == "mole":
                self.gas.TPX = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, self.properties.reactants)
            else:
                self.gas.TPY = (
                    self.properties.temperature, self.properties.pressure * ct.one_atm, self.properties.reactants)

        if self.properties.kind == "constant pressure":
            self.reac = ct.IdealGasConstPressureReactor(self.gas)
        else:
            self.reac = ct.IdealGasReactor(self.gas)

        # Create ``ReactorNet`` newtork
        self.sim = ct.ReactorNet([self.reac])
        self.sample_points = []
        self.ignition_delay = 0.0

    def run_case(self, stop_at_ignition=False, restart=False):
        """Run simulation case set up ``setup_case``.

        If no end time is specified for the integration, the function integrates
        to steady state (or a maximum of 10,000 steps, by default). This is done
        by checking whether the system state changes below a certain threshold,
        with the residual computed using feature checking. This is blatantly stolen
        from Cantera's :meth:`cantera.ReactorNet.advance_to_steady_state` method.

        Parameters
        ----------
        stop_at_ignition : bool
            If ``True``, stop integration at ignition point, don't save data.
        restart : bool
            If ``True``, skip if results file exists.

        Returns
        -------
        self.ignition_delay : float
            Computed ignition delay in seconds

        """

        if restart and os.path.isfile(self.save_file):
            print("Skipped existing case ", self.idx)
            return

        # Save simulation results in hdf5 table format.
        table_def = {
            "time": tables.Float64Col(pos=0),
            "temperature": tables.Float64Col(pos=1),
            "pressure": tables.Float64Col(pos=2),
            "mass_fractions": tables.Float64Col(shape=(self.reac.thermo.n_species), pos=3),
        }

        with tables.open_file(self.save_file, mode="w", title=str(self.idx)) as h5file:

            table = h5file.create_table(where=h5file.root, name="simulation", description=table_def)
            # Row instance to save timestep information to
            timestep = table.row
            # Save initial conditions
            timestep["time"] = self.sim.time
            timestep["temperature"] = self.reac.T
            timestep["pressure"] = self.reac.thermo.P
            timestep["mass_fractions"] = self.reac.Y
            # Add ``timestep`` to table
            timestep.append()

            ignition_flag = False

            # Main time integration loop
            if self.time_end:
                # if end time specified, continue integration until reaching that time
                while self.sim.time < self.time_end:
                    self.sim.step()

                    # Save new timestep information
                    timestep["time"] = self.sim.time
                    timestep["temperature"] = self.reac.T
                    timestep["pressure"] = self.reac.thermo.P
                    timestep["mass_fractions"] = self.reac.Y

                    if self.reac.T >= self.properties.temperature + 400.0 and not ignition_flag:
                        self.ignition_delay = self.sim.time
                        ignition_flag = True

                        if stop_at_ignition:
                            break

                    # Add ``timestep`` to table
                    timestep.append()

            else:
                # otherwise, integrate until steady state, or maximum number of steps reached
                self.sim.reinitialize()
                max_state_values = self.sim.get_state()
                residual_threshold = 10.0 * self.sim.rtol
                absolute_tolerance = self.sim.atol

                for step in range(self.max_steps):
                    previous_state = self.sim.get_state()

                    self.sim.step()

                    # Save new timestep information
                    timestep["time"] = self.sim.time
                    timestep["temperature"] = self.reac.T
                    timestep["pressure"] = self.reac.thermo.P
                    timestep["mass_fractions"] = self.reac.Y

                    if self.reac.T >= self.properties.temperature + 400.0 and not ignition_flag:
                        self.ignition_delay = self.sim.time
                        ignition_flag = True

                        if stop_at_ignition:
                            break

                    # Add ``timestep`` to table
                    timestep.append()

                    state = self.sim.get_state()
                    max_state_values = np.maximum(max_state_values, state)
                    residual = np.linalg.norm(
                        (state - previous_state) / (max_state_values + absolute_tolerance)) / np.sqrt(self.sim.n_vars)

                    if residual < residual_threshold:
                        break

                if step == self.max_steps - 1:
                    logging.error("Maximum number of steps reached before " f"convergence for ignition case {self.idx}")
                    raise RuntimeError(
                        "Maximum number of steps reached before " f"convergence for ignition case {self.idx}")

            # Write ``table`` to disk
            table.flush()

            if not ignition_flag:
                logging.error(f"No ignition detected for ignition case {self.idx}")
                raise RuntimeError(f"No ignition detected for ignition case {self.idx}")

        return self.ignition_delay

    def calculate_ignition(self):
        """Run simulation case set up ``setup_case``, just for ignition delay."""
        # Main time integration loop
        if self.time_end:
            # if end time specified, continue integration until reaching that time
            while self.sim.time < self.time_end:
                self.sim.step()
                if self.reac.T >= self.properties.temperature + 400.0:
                    self.ignition_delay = self.sim.time
                    break
            if not self.ignition_delay:
                logging.warning(f"No ignition detected before end time for ignition case {self.idx}")
        else:
            # otherwise, integrate until steady state, or maximum number of steps reached
            for step in range(self.max_steps):
                self.sim.step()
                if self.reac.T >= self.properties.temperature + 400.0:
                    self.ignition_delay = self.sim.time
                    break
            if step == self.max_steps - 1:
                logging.warning("Maximum number of steps reached before " f"convergence for ignition case {self.idx}")

        return self.ignition_delay

    def process_results(self, skip_data=False):
        """Process integration results to sample data

        Parameters
        ----------
        skip_data : bool
            Flag to skip sampling thermochemical data

        Returns
        -------
        tuple of float, numpy.ndarray or float
            Ignition delay, or ignition delay and sampled data

        """
        delta = 0.05
        # deltas = np.arange(delta, 1 + delta, delta)
        deltas = np.arange(0, 1, delta)

        # Load saved integration results
        self.save_file = os.path.join(self.path, str(self.idx) + ".h5")
        with tables.open_file(self.save_file, "r") as h5file:
            # Load Table with Group name simulation
            table = h5file.root.simulation

            times = table.col("time")
            temperatures = table.col("temperature")
            pressures = table.col("pressure")
            mass_fractions = table.col("mass_fractions")

        # temperature_initial = temperatures[0]
        # temperature_max = temperatures[len(temperatures) - 1]
        # temperature_diff = temperature_max - temperature_initial
        #
        # sampled_data = np.zeros((len(deltas), 2 + mass_fractions.shape[1]))
        #
        # # need to add processing to get the 20 data points here
        # self.ignition_delay = 0.0
        # ignition_flag = False
        # idx = 0
        # for time, temp, pres, mass in zip(times, temperatures, pressures, mass_fractions):
        #     if temp >= temperature_initial + 400.0 and not ignition_flag:
        #         self.ignition_delay = time
        #         ignition_flag = True
        #         if skip_data:
        #             return self.ignition_delay
        #
        #     if temp >= temperature_initial + (deltas[idx] * temperature_diff):
        #         sampled_data[idx, 0:2] = [temp, pres]
        #         sampled_data[idx, 2:] = mass
        #
        #         idx += 1
        #         if idx == 20:
        #             self.sampled_data = sampled_data
        #             return self.ignition_delay, sampled_data

        # change from sampling thermal runaway to radical explosion (modified by Kevin D. on 2024/4/4)
        temperature_initial = temperatures[0]
        sampled_data = np.zeros((len(deltas), 2 + mass_fractions.shape[1]))

        # need to add processing to get the 20 data points here
        self.ignition_delay = 0.0
        ignition_flag = False
        for time, temp, pres, mass in zip(times, temperatures, pressures, mass_fractions):
            if temp >= temperature_initial + 400.0 and not ignition_flag:
                self.ignition_delay = time
                ignition_flag = True
                if skip_data:
                    return self.ignition_delay

        idx = 0
        for time, temp, pres, mass in zip(times, temperatures, pressures, mass_fractions):
            if time >= deltas[idx] * self.ignition_delay:
                sampled_data[idx, 0:2] = [temp, pres]
                sampled_data[idx, 2:] = mass

                idx += 1
                if idx == 20:
                    self.sampled_data = sampled_data
                    return self.ignition_delay, sampled_data

    def clean(self):
        """Delete HDF5 file with full integration data."""
        try:
            os.remove(self.save_file)
        except OSError:
            pass
