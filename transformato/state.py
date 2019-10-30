import transformato
import os
import shutil
from .utils import get_toppar_dir
import logging

logger = logging.getLogger(__name__)


class IntermediateStateFactory(object):

    def __init__(self, system:transformato.system, mutation_list:list, configuration:dict):
        """
        Generate the intermediate directories with for the provided systems with the provided mutations.
        Parameters
        ----------
        system : transformato.system
            definition of the two states for a given system (either waterbox and vacuum or waterbox and complex)
        mutation_list : list
            list of mutations defined by the transformato.ProposeMutationRoute object
        configuration : dict
            configuration dictionary
        """

        self.system = system
        self.mutation_list = mutation_list
        self.path = f"{configuration['analysis_dir_base']}/{self.system.name}"
        self._init_base_dir()
        self.configuration = configuration

    def generate_specific_intermediate_state(self, mutation, state:int):

        output_file_base = self._init_intermediate_state_dir(state)
        logger.info('Writing to {}'.format(output_file_base))
        logger.info('#########################################')
        for psf, offset, env in zip([self.system.complex_psf, self.system.waterbox_psf], [self.system.complex_offset, self.system.waterbox_offset], ['complex', 'waterbox']):
            mutation.mutate(psf, offset, state)
            self._write_psf(psf, output_file_base, env)
        self._write_rtf_file(psf, output_file_base, self.system.tlc)
        self._write_prm_file(psf, output_file_base, self.system.tlc)
        self._write_toppar_str(output_file_base, self.system.tlc)
        self._copy_files(output_file_base)
        return output_file_base


    
    def generate_intermediate_states(self, strategy='seperate'):
        """
        Generate the intermediate states as defined the the mutation list.
        """

        intst_nr = 1
        if strategy == 'seperate':
            # no mixing of the different mutation states - first electrostatics is turend off,
            # then VdW and the the bonded terms are transformed 
            
            for m in self.mutation_list:
                for current_step in range(0, m.nr_of_steps):
                    logger.info('Current step: {}'.format(current_step))
                    output_file_base = self._init_intermediate_state_dir(intst_nr)
                    logger.info('#########################################')
                    logger.info('#########################################')
                    for psf, offset, env in zip([self.system.complex_psf, self.system.waterbox_psf], [self.system.complex_offset, self.system.waterbox_offset], ['complex', 'waterbox']):
                        m.mutate(psf, offset, current_step)
                        self._write_psf(psf, output_file_base, env)
                    self._write_rtf_file(psf, output_file_base, self.system.tlc)
                    self._write_prm_file(psf, output_file_base, self.system.tlc)
                    self._write_toppar_str(output_file_base, self.system.tlc)
                    self._copy_files(output_file_base)
                    intst_nr += 1


    def _copy_files(self, intermediate_state_file_path):
        """
        Copy the files from the original CHARMM-GUI output folder in the intermediate directories.
        """
        # copy crd files
        basedir = self.system.charmm_gui_base
        for env in ['waterbox', 'complex']:
            crd_file_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['crd_file_name']}.crd"
            crd_file_target = f"{intermediate_state_file_path}/lig_in_{env}.crd"
            shutil.copyfile(crd_file_source , crd_file_target)


        # copy rst files
        for env in ['waterbox', 'complex']:
            rst_file_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['rst_file_name']}.rst"
            rst_file_target = f"{intermediate_state_file_path}/lig_in_{env}.rst"
            shutil.copyfile(rst_file_source , rst_file_target)


        # copy ligand rtf file
        ligand_rtf = f"{basedir}/complex/{self.system.tlc.lower()}/{self.system.tlc.lower()}_g.rtf"
        toppar_target = f"{intermediate_state_file_path}/{self.system.tlc.lower()}_g.rtf" 
        shutil.copyfile(ligand_rtf, toppar_target)

        # copy ligand prm file
        ligand_prm = f"{basedir}/complex/{self.system.tlc.lower()}/{self.system.tlc.lower()}.prm"
        toppar_target = f"{intermediate_state_file_path}/{self.system.tlc.lower()}.prm" 
        shutil.copyfile(ligand_prm, toppar_target)



        # copy diverse set of helper functions
        omm_barostat_source = f"{basedir}/complex/openmm/omm_barostat.py"
        omm_barostat_target = f"{intermediate_state_file_path}/omm_barostat.py"
        shutil.copyfile(omm_barostat_source, omm_barostat_target)

        omm_readinputs_source = f"{basedir}/complex/openmm//omm_readinputs.py"
        omm_readinputs_target = f"{intermediate_state_file_path}/omm_readinputs.py"
        shutil.copyfile(omm_readinputs_source, omm_readinputs_target)

        omm_readparams_source = f"{basedir}/complex/openmm/omm_readparams.py"
        omm_readparams_target = f"{intermediate_state_file_path}/omm_readparams.py"
        shutil.copyfile(omm_readparams_source, omm_readparams_target)

        omm_restraints_source = f"{basedir}/complex/openmm/omm_restraints.py"
        omm_restraints_target = f"{intermediate_state_file_path}/omm_restraints.py"
        shutil.copyfile(omm_restraints_source, omm_restraints_target)

        omm_rewrap_source = f"{basedir}/complex/openmm/omm_rewrap.py"
        omm_rewrap_target = f"{intermediate_state_file_path}/omm_rewrap.py"
        shutil.copyfile(omm_rewrap_source, omm_rewrap_target)

        omm_vfswitch_source = f"{basedir}/complex/openmm/omm_vfswitch.py"
        omm_vfswitch_target = f"{intermediate_state_file_path}/omm_vfswitch.py"
        shutil.copyfile(omm_vfswitch_source, omm_vfswitch_target)


        # parse omm simulation paramter
        for env in ['waterbox', 'complex']:
            omm_simulation_parameter_source = f"{basedir}/{env}/openmm/{self.configuration['system'][self.system.structure][env]['simulation_parameter']}" 
            omm_simulation_parameter_target = f"{intermediate_state_file_path}/{self.configuration['system'][self.system.structure][env]['intermediate-filename']}"
            input_simulation_parameter = open(omm_simulation_parameter_source, 'r')
            output_simulation_parameter = open(omm_simulation_parameter_target + '.inp', 'w+')
        
            for l in input_simulation_parameter.readlines():
                if l.strip():
                    t1, t2 = l.split('=')
                    t1 = t1.strip()
                    t2, comment = t2.split('#')
                    t2 = t2.strip()
                    comment = comment.strip()
                    if t1 == 'nstep':
                        t2 = self.configuration['simulation']['nsteps']
                    output_simulation_parameter.write(f"{t1:<25} = {t2:<25} # {comment:<30}\n")
                else:
                    output_simulation_parameter.write('\n')
            input_simulation_parameter.close()
            output_simulation_parameter.close()


        # copy omm simulation script
        omm_simulation_script_source = f"{basedir}/complex/openmm/openmm_run.py"
        omm_simulation_script_target = f"{intermediate_state_file_path}/openmm_run.py"
        shutil.copyfile(omm_simulation_script_source, omm_simulation_script_target)

        # adding serializer functions
        f = open(omm_simulation_script_target, 'a')
        f.write(
'''
# mw: adding xml serializer to the simulation script
file_name = str(args.psffile).replace('.psf', '')
print(file_name)
serialized_integrator = XmlSerializer.serialize(integrator)
outfile = open(file_name + '_integrator.xml','w')
outfile.write(serialized_integrator)
outfile.close()
serialized_system = XmlSerializer.serialize(system)
outfile = open(file_name + '_system.xml','w')
outfile.write(serialized_system)
outfile.close()
'''
        )
        f.close()

        # copy toppar folder
        toppar_dir = get_toppar_dir()
        toppar_source = f"{toppar_dir}"
        toppar_target = f"{intermediate_state_file_path}/toppar" 
        shutil.copytree(toppar_source, toppar_target)

        omm_simulation_submit_script_source = f"{self.configuration['bin_dir']}/simulation.sh"
        omm_simulation_submit_script_target = f"{intermediate_state_file_path}/simulation.sh"
        shutil.copyfile(omm_simulation_submit_script_source, omm_simulation_submit_script_target)  
    
    
    
    def _write_rtf_file(self, psf, output_file_base, tlc): # NOTE: thisneeds some refactoring!
        """
        Generates the dummy atom parameter rtf.
        """

        header_rtf = '''* Dummy atom parameters 
* generated by transformato
*
36  1
'''
        rtf_file_handler = open(output_file_base +'/dummy_atom_definitions.rtf', 'w')
        rtf_file_handler.write(header_rtf)
        for atom in psf.view[f":{tlc}"].atoms:            
            if hasattr(atom, 'real_type'):
                print('- Setting dummy parameters ...')
                print('  + Atom-Name: ', atom.name)
                print('  + Atom-Type: ', atom.real_type)
                print('  + Atom Dummy Type: ', atom.type)

                rtf_file_handler.write('{:7} {:6} {:6} {:6}\n'.format('MASS', '-1', atom.type, atom.mass))
            
        rtf_file_handler.close()    




    def _write_prm_file(self, psf, output_file_base, tlc):
    
        header_prm = '''* Parameters generated by analogy by
* CHARMM General Force Field (CGenFF) program version 1.0.0
*
! Automatically obtained dummy parameters 
! from transformato
'''

        prm_file_handler = open(output_file_base + '/dummy_parameters.prm', 'w')
        prm_file_handler.write(header_prm)
        prm_file_handler.write('\nATOMS\n')
        # for debugging porpose:
        atom_set = set()
        bond_set = set()
        angle_set = set()
        dihedral_set = set()
        improper_set = set()
        nb_set = set()


        view = psf.view[f":{tlc}"]
        # writing atom parameters
        for atom in view.atoms:
            if hasattr(atom, 'real_type'):
                if atom.type in atom_set:
                    continue
                atom_set.add(atom.type)
                print('- Setting dummy parameters ...')
                print('  + Atom-Name: ', atom.name)
                print('  + Atom-Type: ', atom.real_type)
                print('  + Atom Dummy Type: ', atom.type)
                prm_file_handler.write('{:7} {:6} {:6} {:9.5f}\n'.format('MASS', '-1', atom.type, atom.mass))
            
        prm_file_handler.write('\n\n')

        ##############################################################################
        # write bond parameters - again there are two ways to use this:
        # - keeping bonded terms between real/dummy and dummy atoms intact
        # - changing bonded parameters between real atoms - this again needs dummy atoms

        prm_file_handler.write('BONDS\n')
        for bond in view.bonds:
            atom1, atom2 = bond.atom1, bond.atom2
            if hasattr(atom1, 'real_type') or hasattr(atom2, 'real_type'):
                s = frozenset([atom1.type, atom2.type])
            else:
                continue
            if s in bond_set:
                continue
            bond_set.add(s)
            #logger.info(' >> Setting dummy bond parameters for: {} - {}'.format(str(atom1.type),str(atom2.type)))
            prm_file_handler.write('{:7} {:7} {:9.5f} {:9.5f} \n'.format(str(atom1.type), str(atom2.type), bond.type.k ,bond.type.req))


        #################################################################
        prm_file_handler.write('\n\n')
        prm_file_handler.write('ANGLES\n')
        for angle in view.angles:
            atom1, atom2, atom3 = angle.atom1, angle.atom2, angle.atom3
            if hasattr(atom1, 'real_type') or hasattr(atom2, 'real_type') or hasattr(atom3, 'real_type'):
                s = frozenset([atom1.type, atom2.type, atom3.type])
            else:
                continue
            if s in angle_set:
                continue
            angle_set.add(s)

            #print(' >> Setting dummy angle parameters for: {}-{}-{}'.format(str(atom1.type),str(atom2.type),str(atom3.type)))
            prm_file_handler.write('{:7} {:7} {:7} {:9.5f} {:9.5f} \n'.format(str(atom1.type), str(atom2.type), str(atom3.type), angle.type.k , angle.type.theteq))


        #################################################################
        prm_file_handler.write('\n\n')
        prm_file_handler.write('DIHEDRALS\n')
        for dihedral in view.dihedrals:
            atom1, atom2, atom3, atom4 = dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
            if hasattr(atom1, 'real_type') or hasattr(atom2, 'real_type') or hasattr(atom3, 'real_type') or hasattr(atom4, 'real_type'):
                s = frozenset([atom1.type, atom2.type, atom3.type, atom4.type])
            else:
                continue
            if s in dihedral_set:
                continue
            dihedral_set.add(s)

            #print(' >> Setting dummy dihedral parameters for: {}-{}-{}-{}'.format(str(atom1.type),str(atom2.type),str(atom3.type),str(atom4.type)))
            for i in range(len(dihedral.type)):
                prm_file_handler.write('{:7} {:7} {:7} {:7} {:6.5f} {:9.5f} {:9.5f} \n'.format(str(atom1.type), str(atom2.type), str(atom3.type), str(atom4.type), dihedral.type[i].phi_k ,dihedral.type[i].per, dihedral.type[i].phase))

        #################################################################
        # get all unique improper and parameters
        prm_file_handler.write('\n\n')
        prm_file_handler.write('IMPROPERS\n')
        for impr in view.impropers:
            atom1, atom2, atom3, atom4 = impr.atom1, impr.atom2, impr.atom3, impr.atom4
            if hasattr(atom1, 'real_type') or hasattr(atom2, 'real_type') or hasattr(atom3, 'real_type') or hasattr(atom4, 'real_type'):
                s = frozenset([atom1.type, atom2.type, atom3.type, atom4.type])
            else:
                continue

            if s in improper_set:
                continue

            improper_set.add(s)
            
            #print('>> Setting dummy improper parameters for: {}-{}-{}-{}'.format(str(atom1.type),str(atom2.type),str(atom3.type),str(atom4.type)))
            # carefull with this solution - > central atom has to be set in the beginning
            prm_file_handler.write('{:7} {:7} {:7} {:7} {:9.5f} {:9.5f} \n'.format(str(atom1.type), str(atom2.type), str(atom3.type), str(atom4.type), impr.type.psi_k , impr.type.psi_eq))

        #################################################################
        prm_file_handler.write('\n\n')
        prm_file_handler.write('''NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5''')
        prm_file_handler.write('\n\n')

        for atom in view.atoms:
            if not hasattr(atom, 'real_type'):
                continue
            if atom.type in nb_set:
                continue
            nb_set.add(atom.type)
            prm_file_handler.write('{:7} {:6} {:9.5f} {:9.5f}\n'.format(atom.type, 0.0, atom.epsilon, atom.rmin))

        prm_file_handler.write('\n')
        prm_file_handler.write('END')
        prm_file_handler.close()



    def _init_base_dir(self):
        """
        Generates the base directory which all intermediate states are located.
        """
       
        if os.path.isdir(self.path):
            shutil.rmtree(self.path)
            os.makedirs(self.path)
        else:
            os.makedirs(self.path)
    
    def _write_toppar_str(self, output_file_base, tlc):

        toppar_format = """
toppar/top_all36_prot.rtf
toppar/par_all36m_prot.prm
toppar/top_all36_na.rtf
toppar/par_all36_na.prm
toppar/top_all36_carb.rtf
toppar/par_all36_carb.prm
toppar/top_all36_lipid.rtf
toppar/par_all36_lipid.prm
toppar/top_all36_cgenff.rtf
toppar/par_all36_cgenff.prm
toppar/toppar_water_ions.str
toppar/toppar_dum_noble_gases.str
toppar/toppar_all36_prot_d_aminoacids.str
toppar/toppar_all36_prot_fluoro_alkanes.str
toppar/toppar_all36_prot_heme.str
toppar/toppar_all36_prot_na_combined.str
toppar/toppar_all36_prot_retinol.str
toppar/toppar_all36_na_nad_ppi.str
toppar/toppar_all36_na_rna_modified.str
toppar/toppar_all36_lipid_bacterial.str
toppar/toppar_all36_lipid_cardiolipin.str
toppar/toppar_all36_lipid_cholesterol.str
toppar/toppar_all36_lipid_inositol.str
toppar/toppar_all36_lipid_lps.str
toppar/toppar_all36_lipid_miscellaneous.str
toppar/toppar_all36_lipid_model.str
toppar/toppar_all36_lipid_prot.str
toppar/toppar_all36_lipid_pyrophosphate.str
toppar/toppar_all36_lipid_sphingo.str
toppar/toppar_all36_lipid_yeast.str
toppar/toppar_all36_lipid_hmmm.str
toppar/toppar_all36_lipid_detergent.str
toppar/toppar_all36_lipid_ether.str
toppar/toppar_all36_carb_glycolipid.str
toppar/toppar_all36_carb_glycopeptide.str
toppar/toppar_all36_carb_imlab.str
toppar/toppar_all36_label_spin.str
toppar/toppar_all36_label_fluorophore.str
{}_g.rtf
{}.prm
dummy_atom_definitions.rtf
dummy_parameters.prm
""".format(tlc.lower(), tlc.lower())
        
        f = open(output_file_base + '/toppar.str', 'w+')
        f.write(toppar_format)
        f.close()


    def _write_psf(self, psf, output_file_base:str, env:str):
        """
        Writes the new psf.
        """
           
        psf.write_psf(f"{output_file_base}/lig_in_{env}.psf")
        psf.write_pdb(f"{output_file_base}/lig_in_{env}.pdb")
                    


    def _init_intermediate_state_dir(self, nr:int):
        """
        Generates the intermediate state directory.
        """
        output_file_base = f"{self.path}/intst{nr}/" 

        logger.info(' - Created directory: - {}'.format(os.path.abspath(output_file_base)))
        os.makedirs(output_file_base)
        logger.info(' - Writing in - {}'.format(os.path.abspath(output_file_base)))
        return output_file_base