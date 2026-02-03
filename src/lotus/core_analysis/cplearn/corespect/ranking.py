import numpy as np
import networkx as nx
from numba import njit
from numba.typed import List

from ..utils.gen_utils import get_kNN


from ..utils.gen_utils import ensure_list_of_lists,truncate_ng_list

import igraph as ig

import numpy as np

from numba import types
from numba.typed import List





@njit
def init_random_walk(adj_list, init_walk_len=10):
    n = len(adj_list)
    pi = np.ones(n) / n   # uniform start

    for _ in range(init_walk_len):
        new_pi = np.zeros(n)
        for u in range(n):
            neighbors = adj_list[u]
            deg = len(neighbors)
            if deg > 0:
                prob = pi[u] / deg
                for v in neighbors:
                    new_pi[v] += prob
        pi = new_pi
    return pi

@njit
def ascending_walk(ng_list,init_score,times):

    n0=len(init_score)

    walker=[]
    walker_ng_list=[]

    for ell in range(n0):
        x=[]
        for v in ng_list[ell]:
            if init_score[v]>init_score[ell]:
                x.append(v)


        walker.append(x)
        walker_ng_list.append(len(x))


    v_cover_n=np.zeros((n0))

    for rounds in range(times):
        v_cover=np.zeros((n0))

        for u in range(n0):
            v=u
            while walker_ng_list[v]>0:
                v= walker[v][np.random.randint(0, walker_ng_list[v])]

            v_cover[u]=init_score[v]

        v_cover_n=v_cover_n+v_cover


    v_cover_n=v_cover_n/times

    return v_cover_n



def FLOW_rank_optimized(ng_list,init_score,r):



    times=200
    v_cover_n=ascending_walk(ng_list,init_score,times)

    n=len(ng_list)
    final_score={}
    for u in range(n):
        final_score[u]=init_score[u]/(v_cover_n[u]+0.00001)




    return final_score


def FlowRank(knn_list,r,**kwargs):

    print("Initiating FlowRank with r=",r)

    r=int(r)


    ng_list=ensure_list_of_lists(knn_list)
    init_score=init_random_walk(ng_list)

    ng_list_for_ascend=truncate_ng_list(ng_list,r)
    final_score=FLOW_rank_optimized(ng_list_for_ascend,init_score,r)

    return final_score

#FlowRank based on number of hops
@njit
def ascending_walk_count_hops(ng_list,init_score,times):



    n0=len(init_score)
    v_cover=np.zeros((n0))

    walker=[]
    walker_ng_list=[]

    for ell in range(n0):
        x=[]
        for v in ng_list[ell]:
            if init_score[v]>init_score[ell]:
                x.append(v)


        walker.append(x)
        walker_ng_list.append(len(x))



    v_cover_n=np.zeros((n0))

    for rounds in range(times):

        #print("rounds=",rounds)

        v_cover=np.zeros((n0))

        for u in range(n0):
            c=0
            v=u
            while walker_ng_list[v]>0:
                v= walker[v][np.random.randint(0, walker_ng_list[v])]
                c+=1

            v_cover[u]=c

        v_cover_n=v_cover_n+v_cover


    v_cover_n=v_cover_n/times

    return v_cover_n





def FLOW_rank_optimized_hops(X,init_score,r):

    times=200

    n=X.shape[0]
    knn_list, _=get_kNN(X, r)



    # ng_list=[]
    #
    # #Preparing for numba. ng_list is a n*r matrix.
    # for ell in range(n):
    #   x=[]
    #   for v in G.neighbors(ell):
    #     x.append(int(v))
    #
    #   ng_list.append(List(x))

    init_score_numba=np.zeros((n))
    for u in init_score:
        init_score_numba[u]=init_score[u]

    v_cover_n=ascending_walk_count_hops(knn_list,init_score_numba,times)

    final_score={}
    for u in range(n):
        final_score[u]=-v_cover_n[u]




    return final_score



def FlowRank_count_hops(X,ranking_algo_params):
    allowed_params = ['q', 'r']
    for key in ranking_algo_params.keys():
        if key not in allowed_params:
            raise ValueError(f"Unwanted parameter found: {key}")

    q = ranking_algo_params.get('q', 40)
    r = ranking_algo_params.get('r', 20)

    #Get initial density estimation.

    knn_list, _=get_kNN(X, q)
    G=nx.DiGraph()
    for i in range(len(knn_list)):
        for j in knn_list[i]:
            G.add_edge(i,j,weight=1)


    init_score=init_random_walk(G)

    final_score=FLOW_rank_optimized_hops(X,init_score,r)


    return final_score

import numpy as np
from collections import defaultdict

#Standard rankings: Degree and PageRank
def PageRank(X,ranking_algo_params):
    allowed_params = ['q', 'r']
    for key in ranking_algo_params.keys():
        if key not in allowed_params:
            raise ValueError(f"Unwanted parameter found: {key}")

    q = ranking_algo_params.get('q', 40)

    alpha=0.85

    knn_list, _=get_kNN(X, q)
    G=nx.DiGraph()
    for i in range(len(knn_list)):
        for j in knn_list[i]:
            G.add_edge(i,j,weight=1)

    pagerank_scores = nx.pagerank(G,alpha=alpha)

    return pagerank_scores

