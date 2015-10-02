""" PDE solver. """
import os
import pyurdme

class PDESolver(pyurdme.URDMESolver):
    """ PDE solver class. """
    NAME = 'pde'

    def __init__(self, model, report_level=0, error_tolarance=1e-03):
        self.model = model
        self.is_compiled = False
        self.report_level = report_level
        self.model_name = self.model.name
    
    def run(self, number_of_trajectories=1, seed=None, input_file=None, loaddata=False):
        """ Run one simulation of the model.

        number_of_trajectories: How many trajectories should be run.
        seed: the random number seed (incremented by one for multiple runs).
        input_file: the filename of the solver input data file .
        loaddata: boolean, should the result object load the data into memory on creation.

        Returns:
            URDMEResult object.
                or, if number_of_trajectories > 1
            a list of URDMEResult objects
        """
        raise Exception('TODO')
    

    def __getstate__(self):
        """ Save the state of the solver, saves all instance variables
            and reads all the files necessary to compile the solver off
            of the file system and stores it in a separate state variable.
            If the solver model files is specified, it saves that too.
            This is used by Pickle.
        """
        raise Exception('TODO')

    def __setstate__(self, state):
        """ Set all instance variables for the object.  This is used by Pickle.
        """
        raise Exception('TODO')

    def __del__(self):
        """ Deconstructor.  Removes the compiled solver."""
        raise Exception('TODO')

    def serialize(self, filename=None, report_level=0, sopts=None):
        """ Not needed. """
        pass

    def compile(self):
        """ Compile the model. Is this Needed?"""
        pass

    def create_propensity_file(self, file_name=None):
        """ Is this Needed? """
        raise Exception('TODO')

