import numpy as np
import scipy as sp 
import networkx as nx
import pandas as pd
import operator
import OpenPNM
import os

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def extract_section_as_array(file_name, keyword):
    with open(file_name, 'r') as f:
        line = f.readline()
        while line != '':
            if line.strip() == keyword:
                key = f.readline().split();
                num = int(f.readline().split()[0])

                for i in range(0, num):
                    line = f.readline().split()
                    value = [float(x) if is_float(x) else x for x in line]
                    
                    if i==0:
                        data = np.zeros((num, len(value)))

                    data[i, :] = value
    
            line = f.readline()
    
    try:
        return data
    except UnboundLocalError:
        print('Specified keyword not found in file.')
        raise

def add_nodes_to_nx_graph(G, file_name, keyword, remove_dict_key=[]):
    with open(file_name, 'r') as f:
        line = f.readline()
        while line != '':
            if line.strip() == keyword:
                key = f.readline().split();
                num = int(f.readline().split()[0])

                for i in range(0, num):
                    line = f.readline().split()
                    value = [float(x) if is_float(x) else x for x in line]
                    d = dict(zip(key, value))
                    d['coords'] = str(d['x']) + ' ' + str(d['y']) + ' ' + str(d['z'])
                    for k in key:
                        if k in remove_dict_key:
                            del d[k]

                    G.add_node(i, **d)

            line = f.readline()

    return G

def add_edges_to_nx_graph(G, file_name, keyword, remove_dict_key=[]):
    with open(file_name, 'r') as f:
        line = f.readline()
        while line != '':
            if line.strip() == keyword:
                key = f.readline().split()[2:];
                num = int(f.readline().split()[0])

                for i in range(0, num):
                    line = f.readline().split()
                    value = [float(x) if is_float(x) else x for x in line[2:]]
                    id_from = int(line[0]) - 1
                    id_to = int(line[1]) - 1
                    d = dict(zip(key, value))
                    for k in key:
                        if k in remove_dict_key:
                            del d[k]

                    G.add_edge(id_from, id_to, **d)

            line = f.readline()

    return G

def clean_network_and_get_lcc(G):
    # Remove dead ends from the network (unless the node is an inlet/outlet)
    deg = dict(G.degree())
    to_remove = [n for n in deg if ((deg[n] == 1) and (G.node[n]['inlets'] == 0) and G.node[n]['outlets'] == 0)]
    G.remove_nodes_from(to_remove)

    # Get largest connected component
    LG = sorted(nx.connected_component_subgraphs(G, copy = True), key = len, reverse = True)[0]

    # Restore sequential indexing
    mapping = dict(zip(LG.nodes(), range(max(LG.nodes()))))
    LG = nx.relabel_nodes(LG,mapping)
    
    return LG

def write_nx_csv(G, csv_path):
    
    pore_arrays = []
    pore_names  = []
    for k in list(G.nodes(data = True))[0][1]:
        if k not in ['x', 'y', 'z']:
            pore_arrays.append(pd.Series(nx.get_node_attributes(G, k)).sort_index().reset_index(drop = True))
            pore_names.append('pore.' + k)

    throat_arrays = []
    throat_names  = []
    for k in list(G.edges(data = True))[0][2]:
        throat_arrays.append(pd.Series(nx.get_edge_attributes(G, k)).sort_index().reset_index(drop = True))
        throat_names.append('throat.' + k)
        
    # edges = sorted(G.edges(), key = operator.itemgetter(0,1))
    edges = G.edges()
    df = pd.concat([pd.Series([str(x[0]) + ' ' + str(x[1]) for x in edges])] + pore_arrays + throat_arrays, 
                   ignore_index=True, axis=1)
    
    df.columns = ['throat.conns'] + pore_names + throat_names
    df.to_csv(csv_path, index=False)
    
    return

def update_network_data_file(file_name, networkType, node_index, pressure, edge_index, flowrate, permeability_inlet, permeability_outlet):
    os.remove(file_name)
    
    with open('temp.out', 'w+') as f_out:
        # Write out the node pressures        
        for i in range(0, len(pressure)):
            f_out.write(str(pressure[i]) + '\n')

    os.rename('temp.out', file_name)

def hagen_poiseuille_conductance(network, phase, physics, **kwargs):
    C = network['throat.equivalent_conductance']

    return C[network.throats(physics.name)]

def corrected_boundary_pressure(pore_coord, inletPressure, outletPressure, inletLocation, outletLocation, flowDirection):
    if flowDirection=='x':
        pore_position = pore_coord[:, 0]
    elif flowDirection=='y':
        pore_position = pore_coord[:, 1]
    elif flowDirection=='z':
        pore_position = pore_coord[:, 2]

    grad                = (inletPressure - outletPressure)/(inletLocation - outletLocation)
    corrected_pressure  = inletPressure + grad*(pore_position -  inletLocation)

    return corrected_pressure

def stokes_flow(file_name, csv_path, networkType, inletPressure, outletPressure, inletLocation, outletLocation, flowDirection, mu=1.e-3, temp=298.0):
    
    # Ensure that inlet pressure is higher than outlet pressure. Flow is
    # assumed to go from the inlet to the outlet in all cases.
    assert(inletPressure>outletPressure)

    # Enforce that the 'networkType' input variable can take either the values
    # of 'delaunay' or 'modified'
    assert(networkType=='delaunay' or networkType=='modified')

    # Enforce that the 'flowDirection' input variable can take either the
    # values of 'x', 'y' or 'z'.
    assert(flowDirection=='x' or flowDirection=='y' or flowDirection=='z')

    # **********************************************************************#
    
    # REDUNDANT: this is no longer used. Previous version read the output.txt
    # file and added nodes and edges to a networkx graph, which was
    # subsequently converted to a csv file format and read in by OpenPNM.

    # G = nx.Graph()
    # G = add_nodes_to_nx_graph(G, file_name, networkType+'_node', ['diameter'])
    # G = add_edges_to_nx_graph(G, file_name, networkType+'_edge', ['id_particle1', 'id_particle2', 'id_particle3', 'throat_area'])
    # G = clean_network_and_get_lcc(G)

    # csv_path = 'output.csv'
    # csv_path = 'temp.csv'
    # write_nx_csv(G, csv_path)

    # **********************************************************************#

    # Read network directly from the csv file output from
    # tessellation2network.m, and initialise the network, geometry, phase and
    # physics dictionary. Note that the physics dictionary points to the pre-
    # calculated conductance values in the csv file.
    gen_pn  = OpenPNM.Network.GenericNetwork(name = 'network') 
    pn      = OpenPNM.Utilities.IO.CSV.load(filename = csv_path, network = gen_pn)
    geom    = OpenPNM.Geometry.GenericGeometry(network = pn, pores = pn.Ps, throats = pn.Ts)
    water   = OpenPNM.Phases.GenericPhase(network=pn)
    water['pore.temperature'] = temp
    water['pore.viscosity'] = mu
    phys    = OpenPNM.Physics.GenericPhysics(network = pn, phase = water, geometry = geom)
    phys.models.add(propname = 'throat.hydraulic_conductance', model = hagen_poiseuille_conductance)

    # Initialise the Stokes flow algprithm. Add the dirchelet boundary
    # conditions at the inlet and outlet nodes.
    flow = OpenPNM.Algorithms.StokesFlow(network=pn, phase=water)
    corrected_inletPressure = corrected_boundary_pressure(pn['pore.coords'][pn.pores('inlets'), :], 
        inletPressure, outletPressure, inletLocation, outletLocation, flowDirection)
    corrected_outletPressure = corrected_boundary_pressure(pn['pore.coords'][pn.pores('outlets'), :], 
        inletPressure, outletPressure, inletLocation, outletLocation, flowDirection)
    flow.set_boundary_conditions(pores=pn.pores('inlets'),  bctype='Dirichlet', bcvalue=corrected_inletPressure)
    flow.set_boundary_conditions(pores=pn.pores('outlets'), bctype='Dirichlet', bcvalue=corrected_outletPressure)

    # Run Stokes flow algorithm and get results.
    flow.run()
    flow.return_results()

    # os.remove('temp.csv')

    # Calculate permeability from Darcy's Law. Note that the method 'flow.rate'
    # sums all the flows associated with throats containing pores listed in
    # the specified label. These can be positive or negative depending on the
    # direction of the pressure gradient, and hence, care should be taken when
    # labelling throats. Permeability is calculated from flow rate data on the
    # 'inlets' and 'outlets' boundaries. These should be roughly the same.
    Q_in  = sp.absolute(flow.rate(pores=pn.pores('inlets')))
    Q_out = sp.absolute(flow.rate(pores=pn.pores('outlets')))
    
    Lx    = sp.amax(pn['pore.coords'][:, 0]) - sp.amin(pn['pore.coords'][:, 0])
    Ly    = sp.amax(pn['pore.coords'][:, 1]) - sp.amin(pn['pore.coords'][:, 1])
    Lz    = sp.amax(pn['pore.coords'][:, 2]) - sp.amin(pn['pore.coords'][:, 2])
    if flowDirection=='x':
        L = Lx
        A = Ly*Lz
    elif flowDirection=='y':
        L = Ly
        A = Lx*Lz
    elif flowDirection=='z':
        L = Lz
        A = Lx*Ly

    k_in  = Q_in*mu*L/((inletPressure - outletPressure)*A)
    k_out = Q_out*mu*L/((inletPressure - outletPressure)*A)

    # Extract pressures at each node
    ps      = pn.pores()
    idx_ps  = pn['pore.id']
    P       = flow['pore.pressure']

    # Extract flow rate for each throat
    ts      = pn['throat.conns']
    idx_ts  = idx_ps[ts]
    g       = flow['throat.conductance']
    Q       = g*(P[ts][:,0] - P[ts][:,1])

    # Write to output file
    update_network_data_file(file_name, networkType, idx_ps, P, idx_ts, Q, k_in, k_out)

    # Clear the workspace. This is important when looping through several assemblies.
    OpenPNM.Base.Workspace().clear()
