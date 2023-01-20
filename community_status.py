# coding=utf-8


class Status(object):
    """
    To handle several data in one struct.

    Could be replaced by named tuple, but don't want to depend on python 2.6
    """
    node2com = {}
    total_weight = 0
    internals = {}
    degrees = {}
    gdegrees = {}
    def __init__(self):
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])

    def __str__(self):
        return ("node2com : " + str(self.node2com) + " degrees : "
                + str(self.degrees) + " internals : " + str(self.internals)
                + " total_weight : " + str(self.total_weight))

    def copy(self):
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = self.node2com.copy()
        new_status.internals = self.internals.copy()
        new_status.degrees = self.degrees.copy()
        new_status.gdegrees = self.gdegrees.copy()
        new_status.total_weight = self.total_weight

    def init(self, graph, weight, part=None):
        """Initialize the status of a graph with every node in one community"""
        #print('init')
        count = 0
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.total_weight = graph.size(weight=weight)
        #print('self.total_weight :',self.total_weight)
        if part is None:
            #print('Part is None')
            for node in graph.nodes():
#                 print('node :',node)
                self.node2com[node] = count#节点的顺序Id count,一开始每个节点都在单独的社区中
                #print('self.node2com[%d] :%d'%(node, self.node2com[node]))
#                 print('self.node2com[node=%d] : %d'%(node,self.node2com[node]))
                deg = float(graph.degree(node, weight=weight))
#                 print('点%d的degree: %d'%(node,deg))
                if deg < 0:
                    error = "Bad node degree ({})".format(deg)
                    raise ValueError(error)
                self.degrees[count] = deg#根据社区编号存度数/权值 
#                 print('社区%d的degree: %d'%(count,self.degrees[count]))
                self.gdegrees[node] = deg#根据节点存度数/权值
#                 print('点%d的gdegree: %d'%(node,self.gdegrees[node]))
                edge_data = graph.get_edge_data(node, node, default={weight: 0})#字典
                self.loops[node] = float(edge_data.get(weight, 1))#顶点与顶点本身的度数为0，点与点本身无连接
                self.internals[count] = self.loops[node]
                #print('社区%d内部的degree: %d'%(count,self.internals[count]))
#                 print('self.internals[count=%d] : %d'%(count,self.internals[count]))
                count += 1
        else:
            #print('Part is not None')
            for node in graph.nodes():
                #print('node :',node)
                com = part[node]#节点所在的社区编号
#                 print('节点所在的社区编号part[%d] :%d'%(node,com))
                self.node2com[node] = com#存储现在节点所在的社区的编号
#                 print('存储节点所在的社区编号self.node2com[%d] :%d'%(node, self.node2com[node]))
                deg = float(graph.degree(node, weight=weight))
                self.degrees[com] = self.degrees.get(com, 0) + deg #社区的度，如果该社区已经有度数，那就直接加上现在节点的度数，若是该社区原本没有度数，那就赋值为0 外加内
#                 print('社区%d的degree: %d'%(count,self.degrees[count]))
                self.gdegrees[node] = deg#节点的度 外加内
#                 print('点%d的degree: %d'%(node,self.gdegrees[node]))
                inc = 0.
                for neighbor, datas in graph[node].items():#可得到当前节点的neighbors还有当前节点与此邻居连边的权值
                    edge_weight = datas.get(weight, 1)#看该节点和该邻居节点的边权，这是根据是否为有权图决定的，要是有权图就使用有权图里面赋的值，要是无权图，就使用默认值1
                    if edge_weight <= 0:
                        error = "Bad graph type ({})".format(type(graph))
                        raise ValueError(error)
                    if part[neighbor] == com:#若该节点和邻居节点同在一个社区，则更新度/权值
                        if neighbor == node:
                            inc += float(edge_weight)
                        else:
                            inc += float(edge_weight) / 2.#两个点共享一个社区内的边要除以2
                self.internals[com] = self.internals.get(com, 0) + inc#社区内节点的权重和发生改变
#                 print('社区%d内部的degree: %d'%(com,self.internals[com]))
