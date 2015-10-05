""" PDE solver. """
import os
import pyurdme
import sys
import numpy
from scipy.integrate import ode

class PDESolver(pyurdme.URDMESolver):
    """ PDE solver class. """
    NAME = 'pde'

    def __init__(self, model, report_level=0, error_tolarance=1e-06):
        self.model = model
        self.is_compiled = False
        self.report_level = report_level
        self.model_name = self.model.name
        self.error_tolarance = error_tolarance
        
        self.derivative_fn = None
        self.solver_initial_value = None
    
    def run(self, number_of_trajectories=1, seed=None):
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
        if not self.is_compiled:
            self.compile()
                
        solver = ode(self.derivative_fn).set_integrator('vode', method='bdf', order=15, atol=self.error_tolarance)
        solver.set_initial_value(self.solver_initial_value, 0)
        
        result = PDEResult(self.model)
        ndofs = result.U.shape[1]
        for tndx, t in enumerate(self.model.tspan):
            while solver.successful() and solver.t < t:
                solver.integrate(t)
            result.U[tndx:tndx+1,:] = solver.y.reshape(1,ndofs)

        if number_of_trajectories > 1:
            return [result]*number_of_trajectories
        else:
            return result
    
    def compile(self):
        """ Compile the model. Create the derivative function for the ODE solver. """
        self._create_python_derivative_function()
        self._create_solver_initial_value()

    def _create_python_derivative_function(self):
        """ Create the python function for the RHS of the ODE system. """
        sd = self.model.get_solver_datastructure()
        D = sd['D']
        N = sd['N']
        (xlen, ylen) = D.shape
        species_offset={}
        for snum, sname in enumerate(self.model.listOfSpecies.keys()):
            species_offset[sname]=snum
                
        func_str = "def pyurdme_pde_derivative_function(t, x):\n"
        func_str+= "  return [\n"
        for vndx in range(self.model.u0.shape[1]):
            for sndx in range(self.model.u0.shape[0]):
                ndx = self.model.u0.shape[0] * vndx + sndx
                #sys.stderr.write('{0} '.format(ndx))
                suddomain_num = sd['sd'][vndx]
                species_name = self.model.listOfSpecies.keys()[sndx]
                #sys.stderr.write('vndx={0} sndx={1} species_name={2}\n'.format(vndx,sndx,species_name))
                
                terms = []
                rowsum=0
                for yndx in range(ylen):
                    val = sd['D'][yndx,ndx]
                    rowsum+=val
                    if ndx == yndx:
                        continue
                    if val != 0:
                        terms.append("{0}*x[{1}]".format(val, yndx))
                terms.append("{0}*x[{1}]".format(sd['D'][ndx,ndx]+rowsum, ndx))
                for rndx,R in enumerate(self.model.listOfReactions):
                    #sys.stderr.write('rxn {0} '.format(R))
                    # Check if reaction is valid in subdomain
                    if not (self.model.listOfReactions[R].restrict_to == None or
                        (isinstance(self.model.listOfReactions[R].restrict_to, list) and
                            (len(self.model.listOfReactions[R].restrict_to) == 0 or suddomain_num in self.model.listOfReactions[R].restrict_to)) or
                        self.model.listOfReactions[R].restrict_to == suddomain_num):
                        #sys.stderr.write('not the right subdomain   R.restrict_to={0} suddomain_num={1}\n'.format(self.model.listOfReactions[R].restrict_to, suddomain_num))
                        continue
                    # Check if reaction involves this species (reactant or product)
                    if species_name not in self.model.listOfReactions[R].reactants and species_name not in self.model.listOfReactions[R].products:
                        #sys.stderr.write('this species not involved in this reaction    species_name={0} reactants={1} products={2}\n'.format(species_name, self.model.listOfReactions[R].reactants, self.model.listOfReactions[R].products))
                        continue
                    
                    if self.model.listOfReactions[R].massaction:
                        #sys.stderr.write('is massaction ')
                        #sys.stderr.write('marate={0}  '.format(self.model.listOfReactions[R].marate.value))
                        reactnt_terms = [str(self.model.listOfReactions[R].marate.value)]
                        for rname in self.model.listOfReactions[R].reactants.keys():
                            reactnt_ndx = self.model.u0.shape[0] * vndx + species_offset[rname]
                            reactnt_cnt = self.model.listOfReactions[R].reactants[rname]
                            for _ in range(reactnt_cnt):
                                reactnt_terms.append("x[{0}]".format(reactnt_ndx))
                        #sys.stderr.write("term={0}*{1}\n".format(N[sndx,rndx], "*".join(reactnt_terms)))
                        terms.append("{0}*{1}".format(N[sndx,rndx], "*".join(reactnt_terms)))
                    else:
                        #sys.stderr.write('is custom: {0} {1}\n'.format(N[sndx,rndx], self.model.listOfReactions[R].propensity_function))
                        #TODO replace DataFunctions definitions with values
                        terms.append("{0}*({1})".format(N[sndx,rndx], self.model.listOfReactions[R].propensity_function))
                func_str+= "    {0} ,\n".format(" + ".join(terms))
        func_str+= "    ]\n"
        #print func_str
        
        self.derivative_fn = self._exec_python_derivative_function(func_str)
        

    def _exec_python_derivative_function(self, func_str):
        exec func_str
        return pyurdme_pde_derivative_function


    
    def _create_solver_initial_value(self):
        """ Create initial value for the ODE system. """
        N = self.model.u0.shape[0]*self.model.u0.shape[1]
        self.solver_initial_value = numpy.zeros((N,1),dtype=float)
        for vndx in range(self.model.u0.shape[1]):
            for sndx in range(self.model.u0.shape[0]):
                self.solver_initial_value[ self.model.u0.shape[0] * vndx + sndx ] = self.model.u0[sndx, vndx]


    ##################################################################
    def __getstate__(self):
        """ Save the state of the solver, saves all instance variables
            and reads all the files necessary to compile the solver off
            of the file system and stores it in a separate state variable.
            If the solver model files is specified, it saves that too.
            This is used by Pickle.
        """
        return self.__dict__

    def __setstate__(self, state):
        """ Set all instance variables for the object.  This is used by Pickle.
        """
        self.__dict__ = state
        return self
    
    def __del__(self):
        """ Deconstructor. """
        pass

    def serialize(self, filename=None, report_level=0, sopts=None):
        """ Not needed. """
        pass

    def create_propensity_file(self, file_name=None):
        """ Not needed. """
        pass


###########################################################################################

class PDEResult(pyurdme.URDMEResult):
    """ Result object for a PDE simulation, extends the URDMEResult object. """

    def __init__(self, model=None):
        if model is None:
            raise URDMEError("Model must not be None")
        self.model = model
        N = model.u0.shape[0]*model.u0.shape[1]
        self.U = numpy.zeros((len(model.tspan),N),dtype=float)
        self.tspan = model.tspan
        self.stdout = None #kept for consistency
        self.stderr = None


    def __getstate__(self):
        """ Used by pickle to get state when pickling. We need to read the contents of the
        output file since we can't pickle file objects. """
        state = self.__dict__
        state["v2d"] = self.get_v2d()
        state["d2v"] = self.get_d2v()
        return state

    def __setstate__(self, state):
        """ Used by pickle to set state when unpickling. """
        self.__dict__ = state
        return self

    def read_solution(self):
        pass
    
    def get_timespan(self):
        return self.tspan

    def get_species(self, species, timepoints="all", concentration=True):
        """ Returns a slice (view) of the output matrix U that contains one species for the timepoints
            specified by the time index array. The default is to return all timepoints.

            Data is loaded by slicing directly in the hdf5 dataset, i.e. it the entire
            content of the file is not loaded in memory and the U matrix
            is never added to the object.

            if concentration is False (default), the integer, raw, trajectory data is returned,
            if set to True, the concentration (=copy_number/volume) is returned.

        """

        if isinstance(species, pyurdme.Species):
            spec_name = species.name
        else:
            spec_name = species

        U = self.U
        
        species_map = self.model.get_species_map()
        num_species = self.model.get_num_species()

        spec_indx = species_map[spec_name]
        
        num_species = self.model.u0.shape[0]
        Ncells = self.model.u0.shape[1]

        if timepoints  ==  "all":
            Uslice= U[:,(spec_indx*Ncells):(spec_indx*Ncells+Ncells)]
        else:
            Uslice = U[timepoints,(spec_indx*Ncells):(spec_indx*Ncells+Ncells)]

        if not concentration:
            raise Exception("PDE solver can only report concentration")

        # Reorder the dof from dof ordering to voxel ordering
        Uslice = self._reorder_dof_to_voxel(Uslice, num_species=1)

        # Make sure we return 1D slices as flat arrays
        dims = numpy.shape(Uslice)
        if dims[0] == 1:
            Uslice = Uslice.flatten()

        return Uslice


    def __del__(self):
        """ Deconstructor. """
        pass

    def __setupitems__(self, k):
        pass
    def _initialize_sol(self):
        raise Exception("Not implemented for PDE solver")
