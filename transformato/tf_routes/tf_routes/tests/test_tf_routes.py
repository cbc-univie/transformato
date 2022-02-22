"""
Unit and regression test for the tf_routes package.
"""

# Import package, test suite, and other packages as needed
import random
import sys
from rdkit import Chem
from copy import deepcopy
import networkx as nx

from tf_routes import preprocessing
from tf_routes import routes



def test_tf_routes_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "tf_routes" in sys.modules

def test_routes_compare():

    #pdbbind_smiles_selection = ['CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(C)=O', 'CCCCC(NC(=O)CCC(=O)O)P(=O)(O)Oc1ccccc1', 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O', 'CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccc(OP(=O)(O)O)cc1)NC(C)=O', 'CC(=O)NC(Cc1ccc(OP(=O)(O)O)cc1)C(=O)NC(CCC(=O)O)C(=O)N(C)CCCC1CCCC1', 'CCCCC1CCCN(C(=O)C(CCC(=O)O)NC(=O)C(Cc2ccc(OP(=O)(O)O)cc2)NC(C)=O)C1', 'CC(C)CC(NC(=O)C(O)Cc1ccc(O)cc1)C(=O)N1C(C(=O)NC(CO)CCCNC(=N)N)CC2CCC(O)CC21', 'CCC(C)C(NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(=O)C(CC(=O)O)NC(=O)C[NH3+])C(=O)N1CC=CC1C(=O)NCC(=O)NC(CCC(=O)O)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(C)C(=O)O', 'CC(=O)NC1C(NC(=N)N)C=C(C(=O)O)OC1C(O)C(O)CO', 'COC1=C2CC(C)CC(OC)C(O)C(C)/C=C(\\C)C(OC(N)=O)C(OC)/C=C\\C=C(/C)C(=O)NC(=CC1=O)C2=O', 'CC(=O)Nc1ccc(N2C(=O)C3C4CCC(NC(=O)OCC(=O)O)(CC4)C3C2=O)cc1', 'O=c1[nH]cnc2c1ncn2C1OC(CO)C(O)C1O', 'CCCN(CCc1ccccc1)C(=O)C1OC(C(=O)O)=CC([NH3+])C1NC(C)=O', 'Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(N)(=O)O)C(O)C2O)c(=O)[nH]1', 'CN(C)c1cccc2c(S(=O)(=O)NC(CCCNC(=N)N)C(=O)N3CCCCC3CCC(=O)c3nccs3)cccc12', 'N=C(N)NCCCC(NC(=O)C1CC[NH+]2CCC([NH3+])(Cc3ccccc3)C(=O)N12)C(O)C(=O)NCCc1ccccc1', 'N=C(N)c1ccc(CC2CCCCC(Cc3ccc(C(=N)N)cc3)C2=O)cc1', 'CC(=O)Nc1cc(S(=O)(=O)O)cc2cc(S(=O)(=O)O)cc(O)c12', 'CC(=O)NC(C(=O)NC(C(=O)NC(C)C(=O)NC(CO)C(=O)NC(CO)C(N)=O)C(C)C)C(C)O', 'O=S(=O)(O)CC[NH+]1CCOCC1']
    pdbbind_smiles_selection = ['CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(C)=O']
    smiless_selection = ["CC1=CC2=CC=CC=C2N1", "CC(C)(C)C", "CC1=CC=CC=C1", "CC", "CC1=CC=CO1"]

    all_smiles = smiless_selection + pdbbind_smiles_selection

    smiless = random.choices(all_smiles, k=5)

    molarray1 = []
    for smile in smiless:
        mol = Chem.MolFromSmiles(smile)
        molarray1.append(mol)
    molarray2 = molarray1

    for mol1 in molarray1:
        for mol2 in molarray2:

            mol1 = preprocessing.generate_apply_dicts(mol1)
            mol2 = preprocessing.generate_apply_dicts(mol2)

            mol1coreindex, mol2coreindex, hit_ats1, hit_ats2 = preprocessing.get_common_core(mol1, mol2)

            graphmol1 = preprocessing._mol_to_nx_full_weight(mol1)
            graphmol2 = preprocessing._mol_to_nx_full_weight(mol2)

            subg1, G_dummy1 = preprocessing._find_connected_dummy_regions_mol(mol1, mol1coreindex, graphmol1)
            subg2, G_dummy2 = preprocessing._find_connected_dummy_regions_mol(mol2, mol2coreindex, graphmol2)

            terminaldummy1, terminalreal1 = preprocessing._find_terminal_atom(mol1coreindex,  mol1)
            terminaldummy2, terminalreal2 = preprocessing._find_terminal_atom(mol2coreindex,  mol2)

            matchterminal1 = preprocessing._match_terminal_real_and_dummy_atoms(mol1, terminalreal1, terminaldummy1)
            matchterminal2 = preprocessing._match_terminal_real_and_dummy_atoms(mol2, terminalreal2, terminaldummy2)

            order1 = routes._calculate_order_of_LJ_mutations(
            subg1, matchterminal1, graphmol1
        )
            order2 = routes._calculate_order_of_LJ_mutations(
            subg2, matchterminal2, graphmol2
        )
            

            #new mutation order function - code with iteration
            order1new = routes._calculate_order_of_LJ_mutations_new(
            subg1, matchterminal1, graphmol1
        )
            order2new = routes._calculate_order_of_LJ_mutations_new(
            subg2, matchterminal2, graphmol2
        )
            
            #new mutation order function - code with iteration and different steps according to state
            order1newiter = routes._calculate_order_of_LJ_mutations_new_iter(
            subg1, matchterminal1, graphmol1
        )
            order2newiter = routes._calculate_order_of_LJ_mutations_new_iter(
            subg2, matchterminal2, graphmol2
        )

        #new mutation order function - code with iteration and different steps according to state
            order1newiter_change = routes._calculate_order_of_LJ_mutations_new_iter_change(
            subg1, matchterminal1, graphmol1
        )
            order2newiter_change = routes._calculate_order_of_LJ_mutations_new_iter_change(
            subg2, matchterminal2, graphmol2
        )

    

            sortorder1 = deepcopy(order1)
            sortorder1new = deepcopy(order1new)
            sortorder1newiter = deepcopy(order1newiter)
            sortorder1newiter_change = deepcopy(order1newiter_change)

            for i in sortorder1:
                (i.sort())
            for i in sortorder1new:
                (i.sort())
            for i in sortorder1newiter:
                (i.sort())
            for i in sortorder1newiter_change:
                (i.sort())

            #all mutation algorithms have to return an array of same size and with same entries
            assert (sortorder1 == sortorder1new == sortorder1newiter == sortorder1newiter_change)
        

            sortorder2 = deepcopy(order2)
            sortorder2new = deepcopy(order2new)
            sortorder2newiter = deepcopy(order2newiter)
            sortorder2newiter_change = deepcopy(order2newiter_change)

            for i in sortorder2:
                (i.sort())
            for i in sortorder2new:
                (i.sort())
            for i in sortorder2newiter:
                (i.sort())
            for i in sortorder2newiter_change:
                (i.sort())

            #all mutation algorithms have to return an array of same size and with same entries
            assert (sortorder2 == sortorder2new == sortorder2newiter == sortorder2newiter_change)


            subgnodes1 = []
            for i in subg1:
                for j in i:
                    subgnodes1.append(j)
            subgnodes1 = set(subgnodes1)

            subgnodes2 = []
            for i in subg2:
                for j in i:
                    subgnodes2.append(j)
            subgnodes2 = set(subgnodes2)
            order1_flat = [el for region in order1 for el in region]
            order1_flat = set(order1_flat)

            #the number of different entries returned by the mutation algorithms has to be identical to the number of atoms in the dummy regions 
            assert (len(subgnodes1) == len(order1_flat))

            
            order2_flat = [el for region in order2 for el in region]
            order2_flat = set(order2_flat)
                                
            #the number of different entries returned by the mutation algorithms has to be identical to the number of atoms in the dummy regions
            assert (len(subgnodes2) == len(order2_flat))



def test_cycle_checks():
   
    cholesterol = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"
    cortisol = "CC12CCC(=O)C=C1CCC3C2C(CC4(C3CCC4(C(=O)CO)O)C)O"  
    toluene = "CC1=CC=CC=C1"
    ethane = "CC"  

    mol1 = Chem.MolFromSmiles(cholesterol)
    mol2 = Chem.MolFromSmiles(cortisol)
    mol3 = Chem.MolFromSmiles(toluene)
    mol4 = Chem.MolFromSmiles(ethane)


    mol1 = preprocessing.generate_apply_dicts(mol1)
    mol2 = preprocessing.generate_apply_dicts(mol2)
    mol3 = preprocessing.generate_apply_dicts(mol3)
    mol4 = preprocessing.generate_apply_dicts(mol4)

    graphmol1 = preprocessing._mol_to_nx_full_weight(mol1)
    graphmol2 = preprocessing._mol_to_nx_full_weight(mol2)
    graphmol3 = preprocessing._mol_to_nx_full_weight(mol3)
    graphmol4 = preprocessing._mol_to_nx_full_weight(mol4)

    
    cdict_cholesterol, degree_cholesterol = (routes.cycle_checks(graphmol1))
    cdict_cortisol, degree_cortisol = routes.cycle_checks(graphmol2)
    cdict_toluene, degree_toluene = routes.cycle_checks(graphmol3)
    cdict_ethane, degree_ethane = routes.cycle_checks(graphmol4)
    
    from collections import Counter

    countercholesterol = Counter(cdict_cholesterol.values()).most_common()
    countercholesterol = dict(countercholesterol)

    countercortisol = Counter(cdict_cortisol.values()).most_common()
    countercortisol = dict(countercortisol)

    countertoluene = Counter(cdict_toluene.values()).most_common()
    countertoluene = dict(countertoluene)

    counterethane = Counter(cdict_ethane.values()).most_common()
    counterethane = dict(counterethane)

    #check if the number of atoms with a specific cycle participation is correct
    for key, value in countercholesterol.items():
        if (key == 1):
            assert (value == 11)
        if (key == 2):
            assert (value == 6)
        assert(key < 3)
    
    for key, value in countercortisol.items():
        if (key == 1):
            assert (value == 11)
        if (key == 2):
            assert (value == 6)
        assert(key < 3)

    keys = [0, 1, 2, 3]

    for key, value in countertoluene.items():
        if (key == 1):
            assert (value == 6)
        if (key == 0):
            assert (value == 1)
        assert (key < 2)
    

    for key, value in counterethane.items():
        assert (key < 1)


def test_cycle_checks_nx():
   
    cholesterol = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"
    cortisol = "CC12CCC(=O)C=C1CCC3C2C(CC4(C3CCC4(C(=O)CO)O)C)O"  
    toluene = "CC1=CC=CC=C1"
    ethane = "CC"  

    mol1 = Chem.MolFromSmiles(cholesterol)
    mol2 = Chem.MolFromSmiles(cortisol)
    mol3 = Chem.MolFromSmiles(toluene)
    mol4 = Chem.MolFromSmiles(ethane)


    mol1 = preprocessing.generate_apply_dicts(mol1)
    mol2 = preprocessing.generate_apply_dicts(mol2)
    mol3 = preprocessing.generate_apply_dicts(mol3)
    mol4 = preprocessing.generate_apply_dicts(mol4)

    graphmol1 = preprocessing._mol_to_nx_full_weight(mol1)
    graphmol2 = preprocessing._mol_to_nx_full_weight(mol2)
    graphmol3 = preprocessing._mol_to_nx_full_weight(mol3)
    graphmol4 = preprocessing._mol_to_nx_full_weight(mol4)

    
    cdict_cholesterol, degree_cholesterol = routes.cycle_checks(graphmol1)
    cdict_cortisol, degree_cortisol = routes.cycle_checks(graphmol2)
    cdict_toluene, degree_toluene = routes.cycle_checks(graphmol3)
    cdict_ethane, degree_ethane = routes.cycle_checks(graphmol4)
  
    from collections import Counter

    countercholesterol = Counter(cdict_cholesterol.values()).most_common()
    countercholesterol = dict(countercholesterol)

    countercortisol = Counter(cdict_cortisol.values()).most_common()
    countercortisol = dict(countercortisol)

    countertoluene = Counter(cdict_toluene.values()).most_common()
    countertoluene = dict(countertoluene)

    counterethane = Counter(cdict_ethane.values()).most_common()
    counterethane = dict(counterethane)

    #check if the number of atoms with a specific cycle participation is correct

    for key, value in countercholesterol.items():
        if (key == 1):
            assert (value == 11)
        if (key == 2):
            assert (value == 6)
        assert(key < 3)
    
    for key, value in countercortisol.items():
        if (key == 1):
            assert (value == 11)
        if (key == 2):
            assert (value == 6)
        assert(key < 3)

    keys = [0, 1, 2, 3]

    for key, value in countertoluene.items():
        if (key == 1):
            assert (value == 6)
        if (key == 0):
            assert (value == 1)
        assert (key < 2)
    
    

    for key, value in counterethane.items():
        assert (key < 1)


    graph_cholesterol = routes.cycle_checks_nx(graphmol1)
    graph_cortisol = routes.cycle_checks_nx(graphmol2)
    graph_toluene = routes.cycle_checks_nx(graphmol3)
    graph_ethane = routes.cycle_checks_nx(graphmol4)

    #check if the update of weights according to cycle participation is correct
    
    maxweight = 50
    present_weights = set()
    for i in graph_cholesterol.edges.data():      
        present_weights.add(i[2]['weight'])  
    present_weights = sorted(present_weights)
    if len(present_weights) > 1:
        diff_weights = present_weights[-1] - present_weights[-2]
    else:
        diff_weights = 0

    
    for i in graph_cholesterol.edges.data():
        key1 = cdict_cholesterol[i[0]]    
        key2 = cdict_cholesterol[i[1]]    
        assert (i[2]['weight'] == maxweight - diff_weights *  (key1 + key2))

    present_weights = set()
    for i in graph_cortisol.edges.data():     
        present_weights.add(i[2]['weight'])  
    present_weights = sorted(present_weights)
    if len(present_weights) > 1:
        diff_weights = present_weights[-1] - present_weights[-2]
    else:
        diff_weights = 0

    for i in graph_cortisol.edges.data():
        key1 = cdict_cortisol[i[0]]    
        key2 = cdict_cortisol[i[1]]    
        assert (i[2]['weight'] == maxweight - diff_weights *  (key1 + key2))


    present_weights = set()
    for i in graph_toluene.edges.data():    
        present_weights.add(i[2]['weight'])  
    present_weights = sorted(present_weights)
    if len(present_weights) > 1:
        diff_weights = present_weights[-1] - present_weights[-2]
    else:
        diff_weights = 0

    for i in graph_toluene.edges.data():
        key1 = cdict_toluene[i[0]]    
        key2 = cdict_toluene[i[1]]    
        assert (i[2]['weight'] == maxweight - diff_weights *  (key1 + key2))


    present_weights = set()
    for i in graph_ethane.edges.data():
        present_weights.add(i[2]['weight'])  
    present_weights = sorted(present_weights)
    if len(present_weights) > 1:
        diff_weights = present_weights[-1] - present_weights[-2]
    else:
        diff_weights = 0

    for i in graph_ethane.edges.data():
        key1 = cdict_ethane[i[0]]    
        key2 = cdict_ethane[i[1]]    
        assert (i[2]['weight'] == maxweight - diff_weights *  (key1 + key2))



def test_mutations_new_result():
    pdbbind_smiles_selection = ['CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(C)=O', 'CCCCC(NC(=O)CCC(=O)O)P(=O)(O)Oc1ccccc1', 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O', 'CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccc(OP(=O)(O)O)cc1)NC(C)=O', 'CC(=O)NC(Cc1ccc(OP(=O)(O)O)cc1)C(=O)NC(CCC(=O)O)C(=O)N(C)CCCC1CCCC1', 'CCCCC1CCCN(C(=O)C(CCC(=O)O)NC(=O)C(Cc2ccc(OP(=O)(O)O)cc2)NC(C)=O)C1', 'CC(C)CC(NC(=O)C(O)Cc1ccc(O)cc1)C(=O)N1C(C(=O)NC(CO)CCCNC(=N)N)CC2CCC(O)CC21', 'CCC(C)C(NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(=O)C(CC(=O)O)NC(=O)C[NH3+])C(=O)N1CC=CC1C(=O)NCC(=O)NC(CCC(=O)O)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(C)C(=O)O', 'CC(=O)NC1C(NC(=N)N)C=C(C(=O)O)OC1C(O)C(O)CO', 'COC1=C2CC(C)CC(OC)C(O)C(C)/C=C(\\C)C(OC(N)=O)C(OC)/C=C\\C=C(/C)C(=O)NC(=CC1=O)C2=O', 'CC(=O)Nc1ccc(N2C(=O)C3C4CCC(NC(=O)OCC(=O)O)(CC4)C3C2=O)cc1', 'O=c1[nH]cnc2c1ncn2C1OC(CO)C(O)C1O', 'CCCN(CCc1ccccc1)C(=O)C1OC(C(=O)O)=CC([NH3+])C1NC(C)=O', 'Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(N)(=O)O)C(O)C2O)c(=O)[nH]1', 'CN(C)c1cccc2c(S(=O)(=O)NC(CCCNC(=N)N)C(=O)N3CCCCC3CCC(=O)c3nccs3)cccc12', 'N=C(N)NCCCC(NC(=O)C1CC[NH+]2CCC([NH3+])(Cc3ccccc3)C(=O)N12)C(O)C(=O)NCCc1ccccc1', 'N=C(N)c1ccc(CC2CCCCC(Cc3ccc(C(=N)N)cc3)C2=O)cc1', 'CC(=O)Nc1cc(S(=O)(=O)O)cc2cc(S(=O)(=O)O)cc(O)c12', 'CC(=O)NC(C(=O)NC(C(=O)NC(C)C(=O)NC(CO)C(=O)NC(CO)C(N)=O)C(C)C)C(C)O', 'O=S(=O)(O)CC[NH+]1CCOCC1']

    smiless_selection = ["CC1=CC2=CC=CC=C2N1", "CC(C)(C)C", "CC1=CC=CC=C1", "CC", "CC1=CC=CO1"]

    all_smiles = smiless_selection + pdbbind_smiles_selection

    smiless = random.sample(all_smiles, k=2)

    molarray1 = []
    for smile in smiless:
        mol = Chem.MolFromSmiles(smile)
        molarray1.append(mol)
    molarray2 = molarray1

    for mol1 in molarray1:
        for mol2 in molarray2:

            mol1 = preprocessing.generate_apply_dicts(mol1)
            mol2 = preprocessing.generate_apply_dicts(mol2)

            mol1coreindex, mol2coreindex, hit_ats1, hit_ats2 = preprocessing.get_common_core(mol1, mol2)

            graphmol1 = preprocessing._mol_to_nx_full_weight(mol1)
            graphmol2 = preprocessing._mol_to_nx_full_weight(mol2)

            subg1, G_dummy1 = preprocessing._find_connected_dummy_regions_mol(mol1, mol1coreindex, graphmol1)
            subg2, G_dummy2 = preprocessing._find_connected_dummy_regions_mol(mol2, mol2coreindex, graphmol2)

            terminaldummy1, terminalreal1 = preprocessing._find_terminal_atom(mol1coreindex,  mol1)
            terminaldummy2, terminalreal2 = preprocessing._find_terminal_atom(mol2coreindex,  mol2)

            matchterminal1 = preprocessing._match_terminal_real_and_dummy_atoms(mol1, terminalreal1, terminaldummy1)
            matchterminal2 = preprocessing._match_terminal_real_and_dummy_atoms(mol2, terminalreal2, terminaldummy2)

       
            #only necessary for 'illegal' cases, i.e. when a dummy region is multiple times connected to the common core
            matchterminal1 = preprocessing.reduce_terminal(matchterminal1, subg1, G_dummy1)
            matchterminal2 = preprocessing.reduce_terminal(matchterminal2, subg2, G_dummy2)

       

            #new mutation order function - code with iteration
            order1new = routes._calculate_order_of_LJ_mutations_new(
            subg1, matchterminal1, graphmol1
        )
            order2new = routes._calculate_order_of_LJ_mutations_new(
            subg2, matchterminal2, graphmol2
        )

   
            removed_nodes = []
            #remove dummy nodes sequentially
            for nextnodes in order1new:
                for count, nextnode in enumerate(nextnodes):
                  
                    removed_nodes.append(nextnode)
                    graphmol1.remove_node(nextnode)
                   
                    #check that no disconnected components emerge
                    assert(nx.number_connected_components(graphmol1) == 1)

                    #check if finally removed node is root node
                    if (count == len(nextnodes)-1):                       
                        assert ((nextnode in terminaldummy1) == True)

            sublists = []
            for i in subg1:
                sublist = list(i)
                sublists.append(sublist)
            glist = [item for sublist in sublists for item in sublist]

            #finally all nodes from the dummy atom list have to be removed
            assert (set(glist) == set(removed_nodes))


            #same for second molecule
            removed_nodes = []
            #remove dummy nodes sequentially
            for nextnodes in order2new:
                for count, nextnode in enumerate(nextnodes):
                  
                    removed_nodes.append(nextnode)
                    graphmol2.remove_node(nextnode)
                   
                    #check that no disconnected components emerge
                    assert(nx.number_connected_components(graphmol2) == 1)

                    #check if last removed node is root node
                    if (count == len(nextnodes)-1):                       
                        assert ((nextnode in terminaldummy2) == True)

            sublists = []
            for i in subg2:
                sublist = list(i)
                sublists.append(sublist)
            glist = [item for sublist in sublists for item in sublist]

            #finally all nodes from the dummy atom list have to be removed
            assert (set(glist) == set(removed_nodes))


            #after removal off all dummy atoms, i.e. all atoms of the mutation route, the graph representation of both molecules has to be the same (because only the common core is left)
            assert(nx.is_isomorphic(graphmol1, graphmol2))


def test_mutations_new_iter_result():
    pdbbind_smiles_selection = ['CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(C)=O', 'CCCCC(NC(=O)CCC(=O)O)P(=O)(O)Oc1ccccc1', 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O', 'CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccc(OP(=O)(O)O)cc1)NC(C)=O', 'CC(=O)NC(Cc1ccc(OP(=O)(O)O)cc1)C(=O)NC(CCC(=O)O)C(=O)N(C)CCCC1CCCC1', 'CCCCC1CCCN(C(=O)C(CCC(=O)O)NC(=O)C(Cc2ccc(OP(=O)(O)O)cc2)NC(C)=O)C1', 'CC(C)CC(NC(=O)C(O)Cc1ccc(O)cc1)C(=O)N1C(C(=O)NC(CO)CCCNC(=N)N)CC2CCC(O)CC21', 'CCC(C)C(NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(=O)C(CC(=O)O)NC(=O)C[NH3+])C(=O)N1CC=CC1C(=O)NCC(=O)NC(CCC(=O)O)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(C)C(=O)O', 'CC(=O)NC1C(NC(=N)N)C=C(C(=O)O)OC1C(O)C(O)CO', 'COC1=C2CC(C)CC(OC)C(O)C(C)/C=C(\\C)C(OC(N)=O)C(OC)/C=C\\C=C(/C)C(=O)NC(=CC1=O)C2=O', 'CC(=O)Nc1ccc(N2C(=O)C3C4CCC(NC(=O)OCC(=O)O)(CC4)C3C2=O)cc1', 'O=c1[nH]cnc2c1ncn2C1OC(CO)C(O)C1O', 'CCCN(CCc1ccccc1)C(=O)C1OC(C(=O)O)=CC([NH3+])C1NC(C)=O', 'Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(N)(=O)O)C(O)C2O)c(=O)[nH]1', 'CN(C)c1cccc2c(S(=O)(=O)NC(CCCNC(=N)N)C(=O)N3CCCCC3CCC(=O)c3nccs3)cccc12', 'N=C(N)NCCCC(NC(=O)C1CC[NH+]2CCC([NH3+])(Cc3ccccc3)C(=O)N12)C(O)C(=O)NCCc1ccccc1', 'N=C(N)c1ccc(CC2CCCCC(Cc3ccc(C(=N)N)cc3)C2=O)cc1', 'CC(=O)Nc1cc(S(=O)(=O)O)cc2cc(S(=O)(=O)O)cc(O)c12', 'CC(=O)NC(C(=O)NC(C(=O)NC(C)C(=O)NC(CO)C(=O)NC(CO)C(N)=O)C(C)C)C(C)O', 'O=S(=O)(O)CC[NH+]1CCOCC1']

    smiless_selection = ["CC1=CC2=CC=CC=C2N1", "CC(C)(C)C", "CC1=CC=CC=C1", "CC", "CC1=CC=CO1"]

    all_smiles = smiless_selection + pdbbind_smiles_selection

    smiless = random.sample(all_smiles, k=2)

    molarray1 = []
    for smile in smiless:
        mol = Chem.MolFromSmiles(smile)
        molarray1.append(mol)
    molarray2 = molarray1

    for mol1 in molarray1:
        for mol2 in molarray2:

            mol1 = preprocessing.generate_apply_dicts(mol1)
            mol2 = preprocessing.generate_apply_dicts(mol2)

            mol1coreindex, mol2coreindex, hit_ats1, hit_ats2 = preprocessing.get_common_core(mol1, mol2)

            graphmol1 = preprocessing._mol_to_nx_full_weight(mol1)
            graphmol2 = preprocessing._mol_to_nx_full_weight(mol2)

            subg1, G_dummy1 = preprocessing._find_connected_dummy_regions_mol(mol1, mol1coreindex, graphmol1)
            subg2, G_dummy2 = preprocessing._find_connected_dummy_regions_mol(mol2, mol2coreindex, graphmol2)

            terminaldummy1, terminalreal1 = preprocessing._find_terminal_atom(mol1coreindex,  mol1)
            terminaldummy2, terminalreal2 = preprocessing._find_terminal_atom(mol2coreindex,  mol2)

            matchterminal1 = preprocessing._match_terminal_real_and_dummy_atoms(mol1, terminalreal1, terminaldummy1)
            matchterminal2 = preprocessing._match_terminal_real_and_dummy_atoms(mol2, terminalreal2, terminaldummy2)

       
            #only necessary for 'illegal' cases, i.e. when a dummy region is multiple times connected to the common core; do not use for well-defined common cores
            matchterminal1 = preprocessing.reduce_terminal(matchterminal1, subg1, G_dummy1)
            matchterminal2 = preprocessing.reduce_terminal(matchterminal2, subg2, G_dummy2)

       

            #new mutation order function - code with iteration
            order1new = routes._calculate_order_of_LJ_mutations_new_iter(
            subg1, matchterminal1, graphmol1
        )
            order2new = routes._calculate_order_of_LJ_mutations_new_iter(
            subg2, matchterminal2, graphmol2
        )

   
            removed_nodes = []
            #remove dummy nodes sequentially
            for nextnodes in order1new:
                for count, nextnode in enumerate(nextnodes):
                  
                    removed_nodes.append(nextnode)
                    graphmol1.remove_node(nextnode)
                   
                    #check that no disconnected components emerge
                    assert(nx.number_connected_components(graphmol1) == 1)

                    #check if finally removed node is root node
                    if (count == len(nextnodes)-1):                       
                        assert ((nextnode in terminaldummy1) == True)

            sublists = []
            for i in subg1:
                sublist = list(i)
                sublists.append(sublist)
            glist = [item for sublist in sublists for item in sublist]

            #finally all nodes from the dummy atom list have to be removed
            assert (set(glist) == set(removed_nodes))


            #same for second molecule
            removed_nodes = []
            #remove dummy nodes sequentially
            for nextnodes in order2new:
                for count, nextnode in enumerate(nextnodes):
                  
                    removed_nodes.append(nextnode)
                    graphmol2.remove_node(nextnode)
                   
                    #check that no disconnected components emerge
                    assert(nx.number_connected_components(graphmol2) == 1)

                    #check if last removed node is root node
                    if (count == len(nextnodes)-1):                       
                        assert ((nextnode in terminaldummy2) == True)

            sublists = []
            for i in subg2:
                sublist = list(i)
                sublists.append(sublist)
            glist = [item for sublist in sublists for item in sublist]

            #finally all nodes from the dummy atom list have to be removed
            assert (set(glist) == set(removed_nodes))


            #after removal off all dummy atoms, i.e. all atoms of the mutation route, the graph representation of both molecules has to be the same (because only the common core is left)
            assert(nx.is_isomorphic(graphmol1, graphmol2))


def test_mutations_new_iter_change_result():
    pdbbind_smiles_selection = ['CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(C)=O', 'CCCCC(NC(=O)CCC(=O)O)P(=O)(O)Oc1ccccc1', 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O', 'CCCCCN(CCCCC)C(=O)C(CCC(=O)O)NC(=O)C(Cc1ccc(OP(=O)(O)O)cc1)NC(C)=O', 'CC(=O)NC(Cc1ccc(OP(=O)(O)O)cc1)C(=O)NC(CCC(=O)O)C(=O)N(C)CCCC1CCCC1', 'CCCCC1CCCN(C(=O)C(CCC(=O)O)NC(=O)C(Cc2ccc(OP(=O)(O)O)cc2)NC(C)=O)C1', 'CC(C)CC(NC(=O)C(O)Cc1ccc(O)cc1)C(=O)N1C(C(=O)NC(CO)CCCNC(=N)N)CC2CCC(O)CC21', 'CCC(C)C(NC(=O)C(CCC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(Cc1ccccc1)NC(=O)C(CC(=O)O)NC(=O)C[NH3+])C(=O)N1CC=CC1C(=O)NCC(=O)NC(CCC(=O)O)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(C)C(=O)O', 'CC(=O)NC1C(NC(=N)N)C=C(C(=O)O)OC1C(O)C(O)CO', 'COC1=C2CC(C)CC(OC)C(O)C(C)/C=C(\\C)C(OC(N)=O)C(OC)/C=C\\C=C(/C)C(=O)NC(=CC1=O)C2=O', 'CC(=O)Nc1ccc(N2C(=O)C3C4CCC(NC(=O)OCC(=O)O)(CC4)C3C2=O)cc1', 'O=c1[nH]cnc2c1ncn2C1OC(CO)C(O)C1O', 'CCCN(CCc1ccccc1)C(=O)C1OC(C(=O)O)=CC([NH3+])C1NC(C)=O', 'Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(N)(=O)O)C(O)C2O)c(=O)[nH]1', 'CN(C)c1cccc2c(S(=O)(=O)NC(CCCNC(=N)N)C(=O)N3CCCCC3CCC(=O)c3nccs3)cccc12', 'N=C(N)NCCCC(NC(=O)C1CC[NH+]2CCC([NH3+])(Cc3ccccc3)C(=O)N12)C(O)C(=O)NCCc1ccccc1', 'N=C(N)c1ccc(CC2CCCCC(Cc3ccc(C(=N)N)cc3)C2=O)cc1', 'CC(=O)Nc1cc(S(=O)(=O)O)cc2cc(S(=O)(=O)O)cc(O)c12', 'CC(=O)NC(C(=O)NC(C(=O)NC(C)C(=O)NC(CO)C(=O)NC(CO)C(N)=O)C(C)C)C(C)O', 'O=S(=O)(O)CC[NH+]1CCOCC1']

    smiless_selection = ["CC1=CC2=CC=CC=C2N1", "CC(C)(C)C", "CC1=CC=CC=C1", "CC", "CC1=CC=CO1"]

    all_smiles = smiless_selection + pdbbind_smiles_selection

    smiless = random.sample(all_smiles, k=2)

    molarray1 = []
    for smile in smiless:
        mol = Chem.MolFromSmiles(smile)
        molarray1.append(mol)
    molarray2 = molarray1

    for mol1 in molarray1:
        for mol2 in molarray2:

            mol1 = preprocessing.generate_apply_dicts(mol1)
            mol2 = preprocessing.generate_apply_dicts(mol2)

            mol1coreindex, mol2coreindex, hit_ats1, hit_ats2 = preprocessing.get_common_core(mol1, mol2)

            graphmol1 = preprocessing._mol_to_nx_full_weight(mol1)
            graphmol2 = preprocessing._mol_to_nx_full_weight(mol2)

            subg1, G_dummy1 = preprocessing._find_connected_dummy_regions_mol(mol1, mol1coreindex, graphmol1)
            subg2, G_dummy2 = preprocessing._find_connected_dummy_regions_mol(mol2, mol2coreindex, graphmol2)

            terminaldummy1, terminalreal1 = preprocessing._find_terminal_atom(mol1coreindex,  mol1)
            terminaldummy2, terminalreal2 = preprocessing._find_terminal_atom(mol2coreindex,  mol2)

            matchterminal1 = preprocessing._match_terminal_real_and_dummy_atoms(mol1, terminalreal1, terminaldummy1)
            matchterminal2 = preprocessing._match_terminal_real_and_dummy_atoms(mol2, terminalreal2, terminaldummy2)

       
            #only necessary for 'illegal' cases, i.e. when a dummy region is multiple times connected to the common core; do not use for well-defined common cores
            matchterminal1 = preprocessing.reduce_terminal(matchterminal1, subg1, G_dummy1)
            matchterminal2 = preprocessing.reduce_terminal(matchterminal2, subg2, G_dummy2)

       

            #new mutation order function - code with iteration
            order1new = routes._calculate_order_of_LJ_mutations_new_iter_change(
            subg1, matchterminal1, graphmol1
        )
            order2new = routes._calculate_order_of_LJ_mutations_new_iter_change(
            subg2, matchterminal2, graphmol2
        )

   
            removed_nodes = []
            #remove dummy nodes sequentially
            for nextnodes in order1new:
                for count, nextnode in enumerate(nextnodes):
                  
                    removed_nodes.append(nextnode)
                    graphmol1.remove_node(nextnode)
                   
                    #check that no disconnected components emerge
                    assert(nx.number_connected_components(graphmol1) == 1)

                    #check if finally removed node is root node
                    if (count == len(nextnodes)-1):                       
                        assert ((nextnode in terminaldummy1) == True)

            sublists = []
            for i in subg1:
                sublist = list(i)
                sublists.append(sublist)
            glist = [item for sublist in sublists for item in sublist]

            #finally all nodes from the dummy atom list have to be removed
            assert (set(glist) == set(removed_nodes))


            #same for second molecule
            removed_nodes = []
            #remove dummy nodes sequentially
            for nextnodes in order2new:
                for count, nextnode in enumerate(nextnodes):
                  
                    removed_nodes.append(nextnode)
                    graphmol2.remove_node(nextnode)
                   
                    #check that no disconnected components emerge
                    assert(nx.number_connected_components(graphmol2) == 1)

                    #check if last removed node is root node
                    if (count == len(nextnodes)-1):                       
                        assert ((nextnode in terminaldummy2) == True)

            sublists = []
            for i in subg2:
                sublist = list(i)
                sublists.append(sublist)
            glist = [item for sublist in sublists for item in sublist]

            #finally all nodes from the dummy atom list have to be removed
            assert (set(glist) == set(removed_nodes))


            #after removal off all dummy atoms, i.e. all atoms of the mutation route, the graph representation of both molecules has to be the same (because only the common core is left)
            assert(nx.is_isomorphic(graphmol1, graphmol2))
