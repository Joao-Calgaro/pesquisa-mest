#!/usr/bin/env python3

# import packages
import sys
import argparse
import io
import pandas as pd
import numpy as np
from sklearn import metrics
import phy
from ete3 import Tree

sys.setrecursionlimit(100000)

# import file with a table of parameter values and the maximum time of simulation (large number e.g. 500)
parser = argparse.ArgumentParser(description='Generates a tree')
parser.add_argument('inputFile', help='an input file with parameters for the tree')
parser.add_argument('maxTime', type=float, help='an input float for the maximal simulation time')

args = parser.parse_args()

with open(args.inputFile, 'r') as des:
    des_data = des.read()
des.close()

design = pd.read_table(io.StringIO(des_data), index_col='index')
# design = design.loc[:,['R_nought_null', 'transmission_rate', 'removal_rate', 'sample_proba', 'pseudo_sampling_rate']]

#design = design.loc[:, ['R_nought', 'tr_11', 'tr_22', 'tr_12', 'tr_21', 'removal_rate',
#                        'sampling_proba', 'R_nought_1', 'R_nought_2', 'R_nought_verif', 'tree_size', '  nsmission',
#                        'fraction_1', 'infectious_period']]

design = design.loc[:, ['R_nought_1', 'R_nought_2', 'tr_11', 'tr_22', 'tr_12', 'tr_21', 'removal_rate',
                        'sampling_proba', 'infectious_time', 'tree_size', 'proportion_t22_t11',
                        'fraction_1']]


nb_samples = len(design)

# should not be constraining (eg Infinity), depends on your experiment
maxTime = args.maxTime

# initiate attributes of nodes
# reasons why a branch stops (transmission, removed/sampled/unsampled tips before the end of simulations)
STOP_REASON = 'stop_reason'
STOP_UNKNOWN = 0
STOP_TRANSMISSION = 1
STOP_REMOVAL_WOS = 2
STOP_SAMPLING = 3
STOP_TIME = 4
# É para ter apenas um stop_mutation
STOP_MUTATION = 5



# infectious types: their meaning (eg normal/superspreading) depends on the values of parameters
I_1 = 1
I_2 = 2
I_T = 'i_t'

HISTORY = 'history'
SAMPLING = 'sampling'
TRANSMISSION = 'transmission'
DIST_TO_START = 'DIST_TO_START'
PROCESSED = 'processed'



def simulate_bdss_tree_gillespie(tr_r11, tr_r12, tr_r21, tr_r22, removal_r, sampling_p, max_s, max_t,
                                 fraction_1):
    # vou precisar simular R0, gama, a razão do tr22/tr11, tr12, 
    
    """
    Simulates the tree evolution with heterogeneous hosts (of type t1 and t2) based on the given transmission rates,
     removal rate, sampling probabilities and number of tips
    :param tr_r11: float of transmission rate from t1 to t1 type spreader
    :param tr_r12: float of transmission rate from t1 to t2 type spreader
    :param tr_r21: float of transmission rate from t2 to t1 type spreader
    :param tr_r22: float of transmission rate from t2 to t2 type spreader
    :param removal_r: float of removal rate of both t1 and t2 type spreaders
    :param sampling_p: float, between 0 and 1, probability for removed t1 tip/spreader to be immediately sampled
    :param max_s: int, maximum number of sampled leaves in a tree
    :param max_t: float, maximum time from root simulation
    :param fraction_1: float, between 0 and 1, fraction of type 1 spreaders
    :return: the simulated tree (ete3.Tree).
    """
    
    right_size = 0
    trial = 0

    def initialize(init_type):
        if init_type == 1:
            metrics['number_inf1_leaves'] = 1
            leaves_dict['inf1_leaves'] = root.get_leaves()

            metrics['total_inf1'] = 1
        else:
            metrics['number_inf2_leaves'] = 1
            leaves_dict['inf2_leaves'] = root.get_leaves()

            metrics['total_inf2'] = 1
        return None

    def update_rates(rates, metrics_dc):
        """
        updates rates dictionary
        :param rates: dict, all rate values at previous step
        :param metrics_dc: dict, counts of different individuals, list(s) of different branch types
        :return:
        """
        
        #taxa inicial sem mutação
        rates['transmission_rate_1_1_i'] = tr_r11 * metrics_dc['number_inf1_leaves']

        # taxa de surgir uma nova variante:
        # quando surgir uma nova variante, zerar esse parametro para evitar novas mutações
        rates['transmission_rate_1_2_i'] = tr_r12 * metrics_dc['number_inf1_leaves']
        if FLAG_MUTATION == 1:
            rates['transmission_rate_1_2_i'] = 0

        # do 2 pra 1 é zero, ocorreu a mutação
        # rates['transmission_rate_2_1_i'] = tr_r21 * metrics_dc['number_inf2_leaves']
        rates['transmission_rate_2_1_i'] = 0 
        
        #taxa da mutação:
        rates['transmission_rate_2_2_i'] = tr_r22 * metrics_dc['number_inf2_leaves']


        rates['removal_rate_t1_i'] = removal_r * metrics_dc['number_inf1_leaves']
        rates['removal_rate_t2_i'] = removal_r * metrics_dc['number_inf2_leaves']

        return None

    def transmission(t_donor, t_recipient):
        # t_donor is the type of the donor (1 or 2), t_recipient is the type of the recipient (1 or 2)
        # Quando ocorrer a mutação, t_donor = 1 e t_recipient = 2. Depois disso, não pode mais ocorrer a mutação
        global FLAG_MUTATION
        
        if t_donor == 1:
            # which leaf will be affected by the event?
            nb_which_leaf = np.random.randint(0, metrics['number_inf1_leaves'])
            which_leaf = leaves_dict['inf1_leaves'][nb_which_leaf]
            del leaves_dict['inf1_leaves'][nb_which_leaf]
        else:
            # which leaf will be affected by the event?
            nb_which_leaf = np.random.randint(0, metrics['number_inf2_leaves'])
            which_leaf = leaves_dict['inf2_leaves'][nb_which_leaf]
            del leaves_dict['inf2_leaves'][nb_which_leaf]

        # which_leaf becomes an internal node
        # Anotar em qual leaf surge a variante
        # Depois que ocorre a mutação, o parametro de transmissão de t1 para t2 é zerado
        # então a função transmission(1, 2) nunca é chamada dnv (amém)
        #cuidado que o nó que recebe o STOP_REASON=5 possui i_t=1, porque ele é o doador, a mutação ocorre depois da transmissão, então o nó da mutação é do tipo 1, e os descendentes são do tipo 2
        
        if t_donor == 1 and t_recipient == 2:
            which_leaf.dist = abs(time - which_leaf.DIST_TO_START)
            which_leaf.add_feature(DIST_TO_START, time)
            which_leaf.add_feature(STOP_REASON, STOP_MUTATION)
            
            

            # laderalizar a árvore, encontrar a posição das folhas mutantes, 
            # e a partir do formato vetorizado, encontrar a posição do nó da mutação
            # as folhas mutantes tem t=2, preciso achá-las.                                     OK
            # Talvez encontrar todos os nós que possuem stop_reason = 3 (SAMPLING), e t=2       OK
            # o entre os colchetes do CBLV, é o método padrao para passar os atributos dos nós   
        # Aqui entra (1, 1) e (2, 2)
        else:
            which_leaf.dist = abs(time - which_leaf.DIST_TO_START)
            which_leaf.add_feature(DIST_TO_START, time)
            which_leaf.add_feature(STOP_REASON, STOP_TRANSMISSION)
        
        
        
        
        # which_leaf gives birth to recipent and donor
        recipient, donor = which_leaf.add_child(dist=0), which_leaf.add_child(dist=0)
        recipient.add_feature(DIST_TO_START, which_leaf.DIST_TO_START)
        donor.add_feature(DIST_TO_START, which_leaf.DIST_TO_START)

        

        #I_1 e I_2 é o tipo do infectado
        # add recipient to its lists, add it its attributes:
        if t_recipient == 2:
            recipient.add_feature(I_T, I_2)
            leaves_dict['inf2_leaves'].append(recipient)
            metrics['total_inf2'] += 1
            metrics['number_inf2_leaves'] += 1
            FLAG_MUTATION = 1
        else:
            recipient.add_feature(I_T, I_1)
            leaves_dict['inf1_leaves'].append(recipient)
            metrics['total_inf1'] += 1
            metrics['number_inf1_leaves'] += 1

        # add donor to its lists, add it its attributes:
        if t_donor == 1:
            donor.add_feature(I_T, I_1)
            leaves_dict['inf1_leaves'].append(donor)
        else:
            donor.add_feature(I_T, I_2)
            leaves_dict['inf2_leaves'].append(donor)

        metrics['total_branches'] += 1

        
        



        return None

   

    def removal(removed_type):
        """
        updates the tree, the metrics and leaves_dict following a removal event
        :param removed_type: int, either 1: a branch of type 1 undergoes removal; 2: a branch of type 2 undergoes
        removal
        :return:
        """
        # which leaf is removed?: "which_leaf"
        if removed_type == 1:
            nb_which_leaf = np.random.randint(0, metrics['number_inf1_leaves'])
            which_leaf = leaves_dict['inf1_leaves'][nb_which_leaf]
            del leaves_dict['inf1_leaves'][nb_which_leaf]
            metrics['number_inf1_leaves'] -= 1
        else:
            nb_which_leaf = np.random.randint(0, metrics['number_inf2_leaves'])
            which_leaf = leaves_dict['inf2_leaves'][nb_which_leaf]
            del leaves_dict['inf2_leaves'][nb_which_leaf]
            metrics['number_inf2_leaves'] -= 1

        # which_leaf becomes a tip
        which_leaf.dist = abs(time - which_leaf.DIST_TO_START)
        which_leaf.add_feature(DIST_TO_START, time)
        which_leaf.add_feature(PROCESSED, True)

        # was which_leaf sampled?
        if np.random.rand() < sampling_p:
            metrics['number_sampled'] += 1
            which_leaf.add_feature(STOP_REASON, STOP_SAMPLING)
            metrics['total_removed'] += 1
            if removed_type == 1:
                metrics['sampled_inf1'] += 1
            else:
                metrics['sampled_inf2'] += 1
        else:
            which_leaf.add_feature(STOP_REASON, STOP_REMOVAL_WOS)
            metrics['total_removed'] += 1
        return None

    # up to 100 times retrial of simulation until reaching correct size
    while right_size == 0 and trial < 100:

        root = Tree(dist=0)
        root.add_feature(DIST_TO_START, 0)

        # initiate the time of simulation
        time = 0

        # INITIATE: metrics counting leaves and branches of different types, leaves_dict storing all leaves alive, and
        # rates_i with all rates at given time, for Gillespie algorithm
        metrics = {'total_branches': 1, 'total_removed': 0, 'number_sampled': 0, 'total_inf1': 0, 'total_inf2': 0,
                   'sampled_inf1': 0, 'sampled_inf2': 0, 'number_inf1_leaves': 0, 'number_inf2_leaves': 0}
        leaves_dict = {'inf1_leaves': [], 'inf2_leaves': []}

        rates_i = {'removal_rate_t1_i': 0, 'removal_rate_t2_i': 0, 'transmission_rate_1_1_i': 0,
                   'transmission_rate_1_2_i': 0, 'transmission_rate_2_1_i': 0, 'transmission_rate_2_2_i': 0}

        # INITIATE: of which type is the first branch?
        # first individual will be of type 1 with probability frac_1 (frequence of type 1 at equilibrium)
        
        # Fraction 1 será 1 na inicialização 
        # então não precisa da verificação do fraction, da para inicializar direto
        
        initialize(1)

        # simulate while [1] the epidemics do not go extinct, [2] given number of patients were not sampled,
        # [3] maximum time of simulation was not reached
        while (metrics['number_inf1_leaves'] + metrics['number_inf2_leaves']) > 0 \
                and (metrics['number_sampled'] < max_s) and (time < max_t):
            # first we need to re-calculate the rates and take its sum
            update_rates(rates_i, metrics)
            sum_rates_i = sum(rates_i.values())

            # when does the next event take place?
            time_to_next = np.random.exponential(1 / sum_rates_i, 1)[0]
            time = time + time_to_next

            # which event will happen
            random_event = np.random.uniform(0, 1, 1) * sum_rates_i

            if random_event < rates_i['transmission_rate_1_1_i']:
                # there will be a transmission event from t1 to t1 type spreader
                transmission(1, 1)

            elif random_event < (rates_i['transmission_rate_1_1_i'] + rates_i['transmission_rate_1_2_i']):
                # there will be a transmission event from t1 to t2 type spreader
                transmission(1, 2)
                                    # mutação ocorreu, dentro da função transmission, o FLAG_MUTATION é atualizado para 1, e a taxa de transmissão de t1 para t2 é zerada
                                    # agora tenho que salvar em qual nó surgiu a mutação


            elif random_event < (rates_i['transmission_rate_1_1_i'] + rates_i['transmission_rate_1_2_i'] +
                                 rates_i['transmission_rate_2_1_i']):
                # there will be a transmission event from t2 to t1 type spreader
                transmission(2, 1)

            elif random_event < (rates_i['transmission_rate_1_1_i'] + rates_i['transmission_rate_1_2_i'] +
                                 rates_i['transmission_rate_2_1_i'] + rates_i['transmission_rate_2_2_i']):
                transmission(2, 2)

            elif random_event < \
                    (sum_rates_i - rates_i['removal_rate_t2_i']):
                # there will be a removal event of t1 type spreader
                removal(removed_type=1)

            else:
                removal(removed_type=2)
                # there will be a removal event of a t2 type spreader

        # tag non-removed tips at the end of simulation
        for leaflet in root.get_leaves():
            if getattr(leaflet, STOP_REASON, False) != 2 and getattr(leaflet, STOP_REASON, False) != 3:
                leaflet.dist = abs(time - leaflet.DIST_TO_START)
                leaflet.add_feature(DIST_TO_START, time)
                leaflet.add_feature(STOP_REASON, STOP_TIME)

        if metrics['number_sampled'] == max_s:
            right_size = 1
        else:
            trial += 1

    # statistics on the simulation
    del metrics['number_inf1_leaves']
    del metrics['number_inf2_leaves']
    vector_count = list(metrics.values())
    vector_count.extend([time, trial])

    return root, vector_count


def _merge_node_with_its_child(nd, child=None, state_feature=STOP_REASON):
    if not child:
        child = nd.get_children()[0]
    nd_hist = getattr(nd, HISTORY, [(getattr(nd, state_feature, ''), 0)])
    nd_hist += [('!', nd.dist - sum(it[1] for it in nd_hist))] \
               + getattr(child, HISTORY, [(getattr(child, state_feature, ''), 0)])
    child.add_features(**{HISTORY: nd_hist})
    child.dist += nd.dist
    if nd.is_root():
        child.up = None
    else:
        parent = nd.up
        parent.remove_child(nd)
        parent.add_child(child)
    return child


def remove_certain_leaves(tre, to_remove=lambda node: False, state_feature=STOP_REASON):
    """
    Removes all the branches leading to naive leaves from the given tree.
    :param tre: the tree of interest (ete3 Tree)
    [(state_1, 0), (state_2, time_of_transition_from_state_1_to_2), ...]. Branch removals will be added as '!'.
    :param to_remove: a method to check if a leaf should be removed.
    :param state_feature: the node feature to store the state
    :return: the tree with naive branches removed (ete3 Tree) or None is all the leaves were naive in the initial tree.
    """

    for nod in tre.traverse("postorder"):
        # If this node has only one child branch
        # it means that the other child branch used to lead to a naive leaf and was removed.
        # We can merge this node with its child
        # (the child was already processed and either is a leaf or has 2 children).
        if len(nod.get_children()) == 1:
            merged_node = _merge_node_with_its_child(nod, state_feature=state_feature)
            if merged_node.is_root():
                tre = merged_node
        elif nod.is_leaf() and to_remove(nod):
            if nod.is_root():
                return None
            nod.up.remove_child(nod)
    return tre


# PREPARE EXPORT

col = ['tree']
forest_export = pd.DataFrame(index=design.index, columns=col)
full_forest_export = pd.DataFrame(index=design.index, columns=col)
tree_graph_export = pd.DataFrame(index=design.index, columns=col)
cblv_export = pd.DataFrame(index=design.index, columns=col)
tree_arrays = []

col2 = ['total_leaves', 'removed_leaves', 'sampled_leaves', 'total_inf1_leaves', 'total_inf2_leaves',
        'sampled_inf1_leaves', 'sampled_inf2_leaves', 'time_of_simulation', 'nb_trials']
stats_export = pd.DataFrame(index=design.index, columns=col2)

# SIMULATE the trees
FLAG_MUTATION = 0  
for experiment_id in range(nb_samples):
    params = design.iloc[experiment_id, ]

    # simulation
    tr, vector_counter = simulate_bdss_tree_gillespie(tr_r11=params[2], tr_r12=params[4], tr_r21=params[5],
                                                      tr_r22=params[3], removal_r=params[6], sampling_p=params[7],
                                                      max_s=params[9], max_t=maxTime, fraction_1=params[11])
    # for display purposes
    i = 0
    for node in tr.traverse("levelorder"):
        node.name = "n" + str(i)
        i += 1

    # salvar o gráfico laderalizado antes de poda
    
    

    tr.write(outfile="testando_formato_before.nw", format=3, features=['DIST_TO_START', 'stop_reason', 'i_t'], format_root_node=True)
    tr.ladderize(direction=1)
    

    with open("tree_before_removal.txt", "w") as file:
        file.write(tr.get_ascii(show_internal=True, attributes=["name"]))


    if tr is not None:
        full_forest_export.iloc[experiment_id][0] = tr.write(features=['DIST_TO_START', 'stop_reason', 'i_t'],
                                                        format_root_node=True, format=3)
    else:
        full_forest_export.iloc[experiment_id][0] = "NA"

    
    tr = remove_certain_leaves(tr, to_remove=lambda node: getattr(node, STOP_REASON) != STOP_SAMPLING)

    # salvar o gráfico laderalizado depois da poda
    tr.ladderize(direction=1)
    tr.write(outfile="testando_formato_after.nw", format=3, features=['DIST_TO_START', 'stop_reason', 'i_t'], format_root_node=True)


    with open("tree_after_removal.txt", "w") as file:
        file.write(tr.get_ascii(show_internal=True, attributes=["name", "DIST_TO_START", "i_t"]))

    if tr is not None:
        forest_export.iloc[experiment_id][0] = tr.write(features=['DIST_TO_START', 'stop_reason', 'i_t'],
                                                        format_root_node=True, format=3)
    else:
        forest_export.iloc[experiment_id][0] = "NA"

    stats_export.iloc[experiment_id] = vector_counter


    tuple_tree, _ = phy.encode_into_most_recent(tr, 1)
    tree_array = tuple_tree.values.flatten()
    tree_arrays.append(tree_array)

    cblv_export.iloc[experiment_id][0] = tuple_tree

    FLAG_MUTATION = 0  
    
# EXPORT
# subpopulations to export as csv

#stats_export.to_csv(path_or_buf="subpopulations.txt", sep='\t', index=True, header=True)

tree_matrix = np.array(tree_arrays)
rescale_factor = np.array(_)

save_data = {
    'tree_data': tree_matrix,
    'rescale_factor': rescale_factor
}

np.savez_compressed('mutation_deep_learning.npz', **save_data)

# For the pipe : export to stdout
# printando a árvore  
#sys.stdout.write(forest_export.to_csv(sep='\t', index=True, header=True))


cblv_export.to_csv("cblv_export.csv", index=True, header=True)

#CBLV
#tuple_tree, _ = phy.encode_into_most_recent(tr, 1)
#print(tuple_tree)