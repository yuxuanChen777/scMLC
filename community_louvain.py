# -*- coding: utf-8 -*-

from __future__ import print_function

import array

import numbers
import warnings

import networkx as nx
import numpy as np
import math

from community_status import Status

__PASS_MAX = -1
__MIN = 0.0000001


def check_random_state(seed):
    #print('seed :',seed)
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError("%r cannot be used to seed a numpy.random.RandomState"
                     " instance" % seed)


def partition_at_level(dendrogram, level):
    partition = dendrogram[0].copy()
    for index in range(1, level + 1):
        for node, community in partition.items():
            partition[node] = dendrogram[index][community]
    return partition


def modularity(partition, multi_graph, weight='weight'):
    for i in multi_graph.keys():
        if multi_graph[i].is_directed():
            raise TypeError("Bad graph type, use only non directed graph")

            
    sum_res=0
    
    #print('len(multi_graph) :',len(multi_graph))
    
    for i in multi_graph.keys():
        graph=multi_graph[i]
        inc = dict([])
        deg = dict([])
        links = graph.size(weight=weight)
        if links == 0:
            raise ValueError("A graph without link has an undefined modularity")

        for node in graph:
#             print('node :',node)
            com = partition[node]#当前node所在的社区编号
#             print('当前node所在的社区编号 ：',com)
            deg[com] = deg.get(com, 0.) + graph.degree(node, weight=weight)#此node所在社区的度数
            for neighbor, datas in graph[node].items():#寻找当前node的邻居
                edge_weight = datas.get(weight, 1)
                if partition[neighbor] == com:#若是邻居所在的社区和当前node所在的社区相同
#                     if neighbor == node:#邻居和当前node是同一个点
#                         inc[com] = inc.get(com, 0.) + float(edge_weight)#此社区的内部度数/权值加上边权,这个情况不会发生
#                     else:
                     inc[com] = inc.get(com, 0.) + float(edge_weight)#此社区的内部度数/权值加上边权除以2

        res = 0.
        for com in set(partition.values()):
            print('com :',com)
            print('inc[%d] :%d'%(com,inc.get(com,0)))
            print('deg[%d] :%d'%(com,deg.get(com,0)))
             
            res += (inc.get(com, 0.) / (2.*links)) - \
                   (deg.get(com, 0.) / (2. * links)) ** 2
        sum_res +=res
    return sum_res


def best_partition(multi_graph,cell_id_w_dict,
                   partition=None,
                   weight='weight',
                   resolution=1.,
                   randomize=None,
                   random_state=None):
    print('resolution :',resolution)
    dendo = generate_dendrogram(multi_graph,cell_id_w_dict,
                                partition,
                                weight,
                                resolution,
                                randomize,
                                random_state)
    return partition_at_level(dendo, len(dendo) - 1)


def generate_dendrogram(multi_graph,cell_id_w_dict,
                        part_init=None,
                        weight='weight',
                        resolution=1.,
                        randomize=None,
                        random_state=None):
#     print('resolution :',resolution)
    for i in multi_graph.keys():
        G=multi_graph[i]
        if G.is_directed():
            raise TypeError("Bad graph type, use only non directed graph")

    # Properly handle random state, eventually remove old `randomize` parameter
    # NOTE: when `randomize` is removed, delete code up to random_state = ...
    if randomize is not None:
        warnings.warn("The `randomize` parameter will be deprecated in future "
                      "versions. Use `random_state` instead.", DeprecationWarning)
        # If shouldn't randomize, we set a fixed seed to get determinisitc results
        if randomize is False:
            random_state = 0

    # We don't know what to do if both `randomize` and `random_state` are defined
    if randomize and random_state is not None:
        raise ValueError("`randomize` and `random_state` cannot be used at the "
                         "same time")

    random_state = check_random_state(random_state)
    # special case, when there is no link
    # the best partition is everyone in its community
    flag=0
    for i in multi_graph.keys(): #所有网络无边才这么分配
        G = multi_graph[i]
        if G.number_of_edges() != 0:
            flag=1
            break
        elif G.number_of_edges()== 0:
            continue
    if flag==0:
        print('图中的点之间都没有边，每个点都分别属于一个单独的社区')
        part = dict([])
        for i, node in enumerate(multi_graph[0].nodes()):
            part[node] = i
            return [part]

    S=list()
    for i in multi_graph.keys(): #对各个层的status初始化
        name='status'+str(i)
        locals()['v'+str(i)]= Status()
        #print('name :',name)
        S.append(locals()['v'+str(i)])
        current_graph = multi_graph[i]
        locals()['v'+str(i)].init(current_graph, weight,part_init)
    
    if len(multi_graph) > 1:#将不同层的相同node分到相同的社区编号
        for i in multi_graph.keys():
            if i==0:
                continue
            else:
                for node in multi_graph[i].nodes():
                    S[i].node2com[node]=S[0].node2com[node]
    status_list = list()
    best_status_list = list()
    Max_mod = 0 
##################################################################################################################
    print('_one_level_first')
    __one_level(multi_graph, cell_id_w_dict ,S, weight, resolution, random_state)
    new_mod = __modularity(S, resolution) 
    for i in multi_graph.keys():
        partition = __renumber(S[i].node2com)
    status_list.append(partition)
    mod = new_mod
    
    for i in multi_graph.keys():
        if i==0:
            multi_graph[i],cell_id_w_dict = induced_graph(partition, multi_graph[i],cell_id_w_dict, weight)
        else :
            multi_graph[i] = induced_graph01(partition, multi_graph[i], weight)  
    sorted(cell_id_w_dict.items(), key=lambda item:item[1], reverse=False)
    for i in multi_graph.keys(): #对各个层的status初始化
        current_graph = multi_graph[i]
        S[i].init(current_graph, weight, part_init)
    while True:
        print('_one_level_second')
        __one_level(multi_graph, cell_id_w_dict, S, weight, resolution, random_state)
        new_mod = __modularity(S, resolution)
        if new_mod - mod < __MIN:
            break
        for i in multi_graph.keys():
            partition = __renumber(S[i].node2com)
        status_list.append(partition)
        mod = new_mod
        for i in multi_graph.keys():
            if i==0:
                multi_graph[i],cell_id_w_dict = induced_graph(partition, multi_graph[i],cell_id_w_dict, weight)
            else :
                multi_graph[i] = induced_graph01(partition, multi_graph[i], weight)  
        sorted(cell_id_w_dict.items(), key=lambda item:item[1], reverse=False)
        for i in multi_graph.keys(): #对各个层的status初始化
            current_graph = multi_graph[i]
            S[i].init(current_graph, weight, part_init)
##################################################################################################################
    return status_list[:]


def induced_graph(partition, graph, cell_id_w_dict, weight="weight"):
    cell_id_w_dict_induced= dict()
    ret = nx.Graph()
    ret.add_nodes_from(partition.values())
    #print('len(cell_id_w_dict) :',len(cell_id_w_dict))
    for node1, node2, datas in graph.edges(data=True):
        edge_weight = datas.get(weight, 1)
        com1 = partition[node1]
        com2 = partition[node2]
        if com1 not in cell_id_w_dict_induced:
            cell_id_w_dict_induced[com1]=0
        if com2 not in cell_id_w_dict_induced:
            cell_id_w_dict_induced[com2]=0
        cell_id_w_dict_induced[com1]=cell_id_w_dict_induced[com1]+cell_id_w_dict[node1]
        
        cell_id_w_dict_induced[com2]=cell_id_w_dict_induced[com2]+cell_id_w_dict[node2]
        
        
        w_prec = ret.get_edge_data(com1, com2, {weight: 0}).get(weight, 1)
        ret.add_edge(com1, com2, **{weight: w_prec + edge_weight})

    return ret,cell_id_w_dict_induced

############################################################################################
def induced_graph01(partition, graph, weight="weight"):
    ret = nx.Graph()
    ret.add_nodes_from(partition.values())
   
    for node1, node2, datas in graph.edges(data=True):
        edge_weight = datas.get(weight, 1)
        com1 = partition[node1]
        com2 = partition[node2]
        w_prec = ret.get_edge_data(com1, com2, {weight: 0}).get(weight, 1)
        ret.add_edge(com1, com2, **{weight: w_prec + edge_weight})

    return ret
############################################################################################


def __renumber(dictionary):
    values = set(dictionary.values())
#     print('旧社区编号values :',values)
    target = set(range(len(values)))
#     print('新社区编号target :',target)
    if values == target:
        ret = dictionary.copy()
    else:
        # add the values that won't be renumbered
        renumbering = dict(zip(target.intersection(values),
                               target.intersection(values)))
        renumbering.update(dict(zip(values.difference(target),
                                    target.difference(values))))

        ret = {k: renumbering[v] for k, v in dictionary.items()}

    return ret



def __one_level(multi_graph, cell_id_w_dict,S, weight_key, resolution, random_state):
    """Compute one level of communities
    """
    #print('开始执行_one_level')
    modified = True
    nb_pass_done = 0
    cur_mod = __modularity(S, resolution)
    new_mod = cur_mod
    graph=multi_graph[0]#以融合网络为主视图
    status=S[0]
    while modified and nb_pass_done != __PASS_MAX:
        cur_mod = new_mod
        modified = False
        nb_pass_done += 1
        
        for node in cell_id_w_dict:  #开始选择节点
#             print('权值：',cell_id_w_dict[node])
            remove_cost=0
#             print('当前节点node :',node)
            com_node = status.node2com[node]#当前点所属的社区
            degc_totw=list()
            neigh_communities=list()
            for i in multi_graph.keys():
                degc_totw .append( S[i].gdegrees.get(node, 0.) / (S[i].total_weight * 2.))  # 该点的度(总数加起来是重复计算的并且等于边数的两倍)占总度数(总权值)的比例  print('进入_neighcom')
           #存储的是字典
            for i in multi_graph.keys():
                neigh_communities.append( __neighcom(node, multi_graph[i], S[i], weight_key))
            neigh_com = __neighcom_2(node, multi_graph, S, weight_key)
            for i in multi_graph.keys():
                com=S[i].node2com[node]
#                 print('不同层同一个节点的社区号：',com)
                remove_cost += - resolution * neigh_communities[i].get(com,0) + (S[i].degrees.get(com, 0.) - S[i].gdegrees.get(node, 0.)) * degc_totw[i]    #某节点离开属于他的社区，某节点同属于一个社区的邻居节点的边权都去掉和该节点离开所在社区的度数remove_cost之和（remove cost越小，证明该节点不该移动，越大，证明该节点应该移动）
                __remove(node,S[i].node2com[node],neigh_communities[i].get(S[i].node2com[node], 0.), S[i])

            best_com = com_node
            best_increase = 0
            best_w = 1000000000000000
            
            for com, dnc in __randomize(neigh_com.items(), random_state):
#                 print('邻居com :',com,'社区内和该node当邻居边权之和dnc :',dnc)
                incr = remove_cost + resolution * dnc
                for i in  multi_graph.keys():
                    incr += - S[i].degrees.get(com, 0.) * degc_totw[i]

           
                if incr >= best_increase : #改进随机性
                    if incr > best_increase:
                        best_increase = incr
                        best_com = com
#                         print('com :',com)
                        best_w = cell_id_w_dict[com]
                    elif  incr == best_increase and cell_id_w_dict[com] > best_w:
#                         print('ok')
                        best_increase = incr
                        best_com = com
                        best_w = cell_id_w_dict[com]
            for i in multi_graph.keys():
                __insert(node,best_com,neigh_communities[i].get(best_com, 0.), S[i])
            if best_com != com_node:
                modified = True#继续修改
        new_mod = __modularity(S, resolution)
        if new_mod - cur_mod < __MIN:
            break


def __neighcom(node, graph, status, weight_key):
#     print('node :',node)
    weights = {}#字典,每次轮到一个新的node,weights都是重新计算的
    #print('_neighcom')
#     print('type(weights) :',type(weights))
    for neighbor, datas in graph[node].items():
#         print('neighbor :',neighbor)
        if neighbor != node:
            edge_weight = datas.get(weight_key, 1)
#             print('邻居和该点的边权edge_weight :',edge_weight)
            neighborcom = status.node2com[neighbor]#邻居所在的社区
#             print('邻居所在的社区编号 ：',neighborcom)
            weights[neighborcom] = weights.get(neighborcom, 0) + edge_weight 
#             print('weights[%d] :%d'%(neighborcom,weights[neighborcom]))
    return weights

def __neighcom_2(node, multi_graph, S, weight_key):
#     print('node :',node)
    weights = {}#字典,每次轮到一个新的node,weights都是重新计算的
#     print('_neighcom')
#     print('type(weights) :',type(weights))
    for i in multi_graph.keys():
        graph=multi_graph[i]
        for neighbor, datas in graph[node].items():
#             print('neighbor :',neighbor)
            if neighbor != node:
                edge_weight = datas.get(weight_key, 1)
#                 print('邻居和该点的边权edge_weight :',edge_weight)
                neighborcom = S[i].node2com[neighbor]#邻居所在的社区
#                 print('邻居所在的社区编号 ：',neighborcom)
                if neighborcom in weights.keys():
                    weights[neighborcom] += edge_weight 
                else:
                    weights[neighborcom] = weights.get(neighborcom, 0) + edge_weight 
    return weights


def __remove(node, com, weight, status):
    """ Remove node from community com and modify status"""
    status.degrees[com] = (status.degrees.get(com, 0.)
                           - status.gdegrees.get(node, 0.))#社区的度（内外都有）改变
#     print('remove status.degrees[%d] :%d'%(com,status.degrees[com]))
    status.internals[com] = float(status.internals.get(com, 0.) - weight - status.loops.get(node, 0.))#社区内部的度改变
#     print('remove status.internals[%d] :%d'%(com,status.internals[com]))
    status.node2com[node] = -1
#     print('remove status.node2com[%d] :%d'%(node,status.node2com[node]))


def __insert(node, com, weight, status):
    """ Insert node into community and modify status"""
    status.node2com[node] = com
#     print('insert后该node%d所属的社区 :%d'%(node,com))
    status.degrees[com] = (status.degrees.get(com, 0.) + status.gdegrees.get(node, 0.))#和这个社区有关的所有的度和/权值和
#     print('insert status.degrees[%d] :%d'%(com,status.degrees[com]))
    status.internals[com] = float(status.internals.get(com, 0.) + weight + status.loops.get(node, 0.))#和这个社区的内部度和/权值和
#     print('insert status.internals[%d] :%d'%(com,status.internals[com]))

def __modularity(S, resolution):
    """
    Fast compute the modularity of the partition of the graph using
    status precomputed
    """
    cnt=0
    result_sum=0
    for status in S:
        #print('cnt:',cnt)
        
        links = float(status.total_weight)#图中的总边数or总权值
        result = 0.
        for community in set(status.node2com.values()):
            in_degree = status.internals.get(community, 0.)
            degree = status.degrees.get(community, 0.)
            if links > 0:
                 #result += possibility[cnt]*(in_degree * resolution / links -  ((degree / (2. * links)) ** 2))#尝试调整resolution乘的位置
                result += (in_degree * resolution / links -  ((degree / (2. * links)) ** 2))#尝试调整resolution乘的位置
                #result += in_degree / links -  (resolution*(degree / (2. * links)) ** 2)#resolution可能需要改进的地方******
        result_sum +=result
      
    #print('_modularity :',result_sum)
    return result_sum


def __randomize(items, random_state):
    """Returns a List containing a random permutation of items"""
    randomized_items = list(items)#图中点的list
    random_state.shuffle(randomized_items)#打乱list中点的顺序
    return randomized_items

