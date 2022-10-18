def AAVA(G, l): 
    nodes = list(nx.nodes(G)) 
    degree_li = nx.degree(G) 
    d_max = max([i[1] for i in degree_li]) 
    degree_dic = {}
    rank = [] 
    for i in degree_li: 
         degree_dic[i[0]] = i[1] 
    node_ability = {} 
    
    for item in degree_li: 
         degree_total = 0
         neighbors1 = list(G.neighbors(item[0]))
         for nbrr in neighbors1: 
             degree_total += degree_dic[nbrr]
         degree = item[1] 
         node_ability[item[0]] =degree_total/d_max
              
    node_score = get_node_score2(G, nodes, node_ability, degree_dic, d_max,rank) 
 
    for i in range(l): 
        max_score_node, score = max(node_score.items(), key=lambda x: x[1]) 
        rank.append((max_score_node, score)) 
        node_ability[max_score_node] = 0
       
        node_score.pop(max_score_node) 
        cur_nbrs = list(nx.neighbors(G, rank[-1][0])) 
        next_cur_neigh = [] 
        jac_total1 = 0
        jac_total2 = 0
        for nbr_jac in cur_nbrs:
            jac_total1 += simjkd(nbr_jac,max_score_node) 
        for nbr in cur_nbrs: 
             jac1 = simjkd(nbr,max_score_node)
             nnbr = nx.neighbors(G, nbr) 
             next_cur_neigh.extend(nnbr) 
             node_ability[nbr] *= 1-jac1           
 
        next_cur_neighs = list(set(next_cur_neigh))  
        for ih in rank: 
             if ih[0] in next_cur_neighs: 
                 next_cur_neighs.remove(ih[0]) 
        for i in cur_nbrs: 
             if i in next_cur_neighs: 
                 next_cur_neighs.remove(i) 
 
 
        for nnbr_jac in next_cur_neighs:
             jac_total2 += simjkd(nnbr_jac,max_score_node) 
        for nnbr in next_cur_neighs: 
             jac2 = simjkd(nnbr,max_score_node) 
             node_ability[nnbr] *= 1-jac2
             
        X = [] 
        X.extend(cur_nbrs) 
        X.extend(next_cur_neighs) 
        for nbr in next_cur_neighs:
            nbrs = nx.neighbors(G, nbr) 
            X.extend(nbrs) 

 
        X = list(set(X)) 
        for ih in rank: 
             if ih[0] in X: 
                 X.remove(ih[0])
        new_nodeScore = get_node_score2(G, X, node_ability, degree_dic, d_max, rank) 
        print(new_nodeScore)
        node_score.update(new_nodeScore) 
        
    return rank 
  
  def get_node_score2(G, nodesNeedcalcu, node_ability, degree_dic, d_max,rank): 
     weight = get_weight(G, degree_dic, d_max, rank) 
     node_score = {} 
     for node in nodesNeedcalcu:  
         sum2 = 0
         degree_total=0
         neighbors = list(nx.neighbors(G, node))
         
         for nbrr in neighbors: 
             degree_total += degree_dic[nbrr]
         for nbr in neighbors: 
             sum2 += node_ability[nbr] *weight[(nbr, node)]
         node_score[node] = math.sqrt(sum2)+len(neighbors) * degree_total/d_max
     return node_score 
  
  def get_weight(G, degree_dic,d_max, rank):  
    weight = {}
    nodes = nx.nodes(G)
    rank_list = [i[0] for i in rank]
    for node in nodes:
        sum1 = 0 
        sum2 = 0
        neighbors = list(nx.neighbors(G, node)) 
        neighbors_common_rank = list(set(neighbors) & set(rank_list)) 
        if len(neighbors_common_rank) != 0:  
             for nc in neighbors_common_rank: 
                 weight[(node, nc)] = degree_dic[node]/d_max
        neighbours_without_rank = list(set(neighbors) - set(rank_list))   
        if len(neighbours_without_rank) != 0:  
             for nbr in neighbours_without_rank: 
                 sum2 += simjkd(nbr,node)
             for neigh in neighbours_without_rank: 
                 jac = simjkd(neigh,node)
                 weight[(node, neigh)] = jac/sum2
        else:  
             for neigh in neighbors: 
                 weight[(node, neigh)] = 0
    return weight 

def simjkd(u, v):        
        set_v = set( G.neighbors(v))
        set_v.add(v)
        set_u = set( G.neighbors(u))
        set_u.add(u)
        jac = len(set_v & set_u) * 1.0 / len(set_v | set_u)    
        return jac
