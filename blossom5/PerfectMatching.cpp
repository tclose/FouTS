#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "blossom5/PMimplementation.h"
#include "blossom5/MinCost/MinCost.h"

//---------//
//  Duals  //
//---------//

void PerfectMatching::ComputeEpsGlobal() {
    Node* r;
    PriorityQueue<REAL>::Item* q;
    Tree* t;
    Tree* t2;
    TreeEdge* e;
    int i, j, k, N = 0, E = 0;
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        t->id = N;
        N += 2;
        for (k = 0; k < 2; k++)
            for (e = t->first[k]; e; e = e->next[k])
                E += 6;
    }
    DualMinCost<REAL>* m = new DualMinCost<REAL>(N, E);
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        i = t->id;
        m->AddUnaryTerm(i, -1);
        m->SetLowerBound(i, 0);
        m->AddUnaryTerm(i + 1, 1);
        m->SetUpperBound(i + 1, 0);
        
        if (t->eps_delta < PM_INFTY) {
            m->SetUpperBound(i, t->eps_delta);
            m->SetLowerBound(i + 1, -t->eps_delta);
        }
        for (e = t->first[0]; e; e = e->next[0]) {
            t2 = e->head[0];
            if (t2 == NULL)
                continue;
            j = e->head[0]->id;
            if ((q = e->pq01[0].GetMin())) {
                m->AddConstraint(j, i, q->slack - t->eps + t2->eps);
                m->AddConstraint(i + 1, j + 1, q->slack - t->eps + t2->eps);
            }
            if ((q = e->pq01[1].GetMin())) {
                m->AddConstraint(i, j, q->slack - t2->eps + t->eps);
                m->AddConstraint(j + 1, i + 1, q->slack - t2->eps + t->eps);
            }
            if ((q = e->pq00.GetMin())) {
                m->AddConstraint(i + 1, j, q->slack - t->eps - t2->eps);
                m->AddConstraint(j + 1, i, q->slack - t->eps - t2->eps);
            }
        }
    }
    m->Solve();
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        i = t->id;
        t->eps_delta = (m->GetSolution(i) - m->GetSolution(i + 1)) / 2;
    }
    delete m;
}

void PerfectMatching::ComputeEpsSingle() {
    Node* r;
    PriorityQueue<REAL>::Item* q;
    Tree* t;
    Tree* t2;
    TreeEdge* e;
    REAL eps = PM_INFTY;
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        if (eps > t->eps_delta)
            eps = t->eps_delta;
        for (e = t->first[0]; e; e = e->next[0]) {
            t2 = e->head[0];
            if ((q = e->pq00.GetMin()) && 2 * eps > q->slack - t->eps - t2->eps) {
                eps = (q->slack - t->eps - t2->eps) / 2;
            }
        }
    }
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        r->tree->eps_delta = eps;
    }
}

void PerfectMatching::ComputeEpsCC() {
    Node* r;
    PriorityQueue<REAL>::Item* q;
    Tree* t0;
    Tree* t;
    Tree* t2;
    Tree* t_next;
    TreeEdge* e;
    REAL eps, eps2;
    Tree* queue_last;
    int dir;
    Tree* FIXED_TREE = trees - 1;
    int component_num = 0;
    TreeEdge** e_ptr;
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        t0->next = NULL;
    }
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        if (t0->next)
            continue;
        eps = t0->eps_delta;
        
        t0->next = queue_last = t = t0;
        while (1) {
            for (dir = 0; dir < 2; dir++)
                for (e_ptr = &t->first[dir], e = *e_ptr; e; e = *e_ptr) {
                    t2 = e->head[dir];
                    if (t2 == NULL) {
                        *e_ptr = e->next[dir];
                        tree_edges->Delete(e);
                        continue;
                    }
                    e_ptr = &e->next[dir];
                    
                    REAL eps00 =
                            ((q = e->pq00.GetMin())) ? (q->slack - t->eps - t2->eps) : PM_INFTY;
                    if (t2->next && t2->next != FIXED_TREE) {
                        if (2 * eps > eps00)
                            eps = eps00 / 2;
                        continue;
                    }
                    
                    REAL eps01[2];
                    eps01[dir] =
                            ((q = e->pq01[dir].GetMin())) ? (q->slack - t->eps + t2->eps) :
                                                            PM_INFTY;
                    eps01[1 - dir] =
                            ((q = e->pq01[1 - dir].GetMin())) ? (q->slack - t2->eps + t->eps) :
                                                                PM_INFTY;
                    
                    if (t2->next == FIXED_TREE)
                        eps2 = t2->eps_delta;
                    else if (eps01[0] > 0 && eps01[1] > 0)
                        eps2 = 0;
                    else {
                        queue_last->next = t2;
                        queue_last = t2;
                        t2->next = t2;
                        if (eps > eps00)
                            eps = eps00;
                        if (eps > t2->eps_delta)
                            eps = t2->eps_delta;
                        continue;
                    }
                    if (eps > eps00 - eps2)
                        eps = eps00 - eps2;
                    if (eps > eps2 + eps01[dir])
                        eps = eps2 + eps01[dir];
                }
            
            if (t->next == t)
                break;
            t = t->next;
        }
        for (t = t0;; t = t_next) {
            t->eps_delta = eps;
            t_next = t->next;
            t->next = FIXED_TREE;
            if (t_next == t)
                break;
        }
        component_num++;
    }
    //printf("%d CCs ", component_num);
}

void PerfectMatching::ComputeEpsSCC() {
    PriorityQueue<REAL>::Item* q;
    Node* r;
    Tree* t0;
    Tree* t;
    Tree* t2;
    TreeEdge* e;
    TreeEdge** e_ptr;
    REAL eps;
    int c, dir;
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        t0->dfs_parent = NULL;
        
        for (dir = 0; dir < 2; dir++)
            for (e_ptr = &t0->first[dir], e = *e_ptr; e; e = *e_ptr) {
                t2 = e->head[dir];
                if (t2 == NULL) {
                    *e_ptr = e->next[dir];
                    tree_edges->Delete(e);
                    continue;
                }
                e_ptr = &e->next[dir];
            }
    }
    Tree* stack = NULL;
    
    // first DFS
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t0 = r->tree;
        if (t0->dfs_parent)
            continue;
        t = t0;
        e = (t->first[0]) ? t->first[0] : t->first[1];
        t->dfs_parent = (TreeEdge*) trees;
        while (1) {
            if (e == NULL) {
                t->next = stack;
                stack = t;
                
                if (t == t0)
                    break;
                
                e = t->dfs_parent;
                if (t == e->head[0]) {
                    t = e->head[1];
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                } else {
                    t = e->head[0];
                    e = e->next[1];
                }
                continue;
            }
            
            if (e->head[1] == t) {
                if (e->head[0]->dfs_parent || !(q = e->pq01[0].GetMin())
                    || q->slack - t->eps + e->head[0]->eps > 0) {
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                    continue;
                }
                t = e->head[0];
            } else {
                if (e->head[1]->dfs_parent || !(q = e->pq01[1].GetMin())
                    || q->slack - t->eps + e->head[1]->eps > 0) {
                    e = e->next[1];
                    continue;
                }
                t = e->head[1];
            }
            t->dfs_parent = e;
            e = (t->first[0]) ? t->first[0] : t->first[1];
        }
    }
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next)
        r->tree->dfs_parent = NULL;
    
    int component_num = 0;
    while (stack) {
        t0 = stack;
        stack = t0->next;
        if (t0->dfs_parent)
            continue;
        t = t0;
        e = (t->first[0]) ? t->first[0] : t->first[1];
        t->dfs_parent = (TreeEdge*) trees;
        while (1) {
            if (e == NULL) {
                e = t->dfs_parent;
                t->dfs_parent = (TreeEdge*) ((char*) trees + component_num);
                if (t == t0)
                    break;
                if (t == e->head[0]) {
                    t = e->head[1];
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                } else {
                    t = e->head[0];
                    e = e->next[1];
                }
                continue;
            }
            
            if (e->head[1] == t) {
                if (e->head[0]->dfs_parent || !(q = e->pq01[1].GetMin())
                    || q->slack - e->head[0]->eps + t->eps > 0) {
                    e = (e->next[0]) ? e->next[0] : t->first[1];
                    continue;
                }
                t = e->head[0];
            } else {
                if (e->head[1]->dfs_parent || !(q = e->pq01[0].GetMin())
                    || q->slack - e->head[1]->eps + t->eps > 0) {
                    e = e->next[1];
                    continue;
                }
                t = e->head[1];
            }
            t->dfs_parent = e;
            e = (t->first[0]) ? t->first[0] : t->first[1];
        }
        component_num++;
    }
    
    Tree** array = new Tree*[component_num];
    memset(array, 0, component_num * sizeof(Tree*));
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        t->id = (int) ((char*) t->dfs_parent - (char*) trees);
        t->next = array[t->id];
        array[t->id] = t;
    }
    
    for (c = component_num - 1; c >= 0; c--) {
        eps = PM_INFTY;
        for (t = array[c]; t; t = t->next) {
            if (eps > t->eps_delta)
                eps = t->eps_delta;
            FOR_ALL_TREE_EDGES(t, e, dir)
            {
                t2 = e->head[dir];
                REAL eps00 = (q = e->pq00.GetMin()) ? (q->slack - t->eps - t2->eps) : PM_INFTY;
                REAL eps01[2];
                eps01[dir] =
                        ((q = e->pq01[dir].GetMin())) ? (q->slack - t->eps + t2->eps) : PM_INFTY;
                eps01[1 - dir] =
                        ((q = e->pq01[1 - dir].GetMin())) ? (q->slack - t2->eps + t->eps) :
                                                            PM_INFTY;
                if (t2->id < c) {
                    if (eps > eps01[dir])
                        eps = eps01[dir];
                    if (eps > eps00)
                        eps = eps00;
                } else if (t2->id == c) {
                    if (2 * eps > eps00)
                        eps = eps00 / 2;
                } else {
                    if (eps > eps01[dir] + t2->eps_delta)
                        eps = eps01[dir] + t2->eps_delta;
                    if (eps > eps00 - t2->eps_delta)
                        eps = eps00 - t2->eps_delta;
                }
            }
        }
        for (t = array[c]; t; t = t->next)
            t->eps_delta = eps;
    }
    
    delete[] array;
    //printf("%d SCCs ", component_num);
}

void PerfectMatching::CommitEps() {
    printf("CommitEps()\n");
    Node* i;
    Node* j;
    Node* r;
    int dir;
    Edge* a;
    EdgeIterator I;
    Tree* t;
    TreeEdge* e;
    TreeEdge** e_ptr;
    REAL eps, eps2;
    PriorityQueue<REAL>::Item* q;
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        t = r->tree;
        eps = t->eps;
        
        i = r;
        while (1) {
            i->y += eps;
            if (!i->is_tree_root) {
                Node* i0 = i;
                i = ARC_HEAD(i0->match) ;
                if (i->is_blossom)
                    ARC_TO_EDGE_PTR(i0->match) ->slack -= eps;
                else
                    i->y -= eps;
                FOR_ALL_EDGES(i, a, dir, I)
                {
                    GET_OUTER_HEAD(a, dir, j);
                    
                    a->slack += eps;
                    if (j->flag == 0)
                        a->slack -= j->tree->eps;
                }
                i = i0;
            }
            
            MOVE_NODE_IN_TREE(i);
        }
        
        t->pq0.Update(-eps);
        
        PriorityQueue<REAL> pq00 = t->pq00;
        t->pq00.Reset();
        for (q = pq00.GetAndResetFirst(); q; q = pq00.GetAndResetNext()) {
            a = (Edge*) q;
            if (ProcessEdge00(a))
                t->pq00.Add(a);
        }
        
        for (e_ptr = &t->first[0], e = *e_ptr; e; e = *e_ptr) {
            if (e->head[0] == NULL) {
                *e_ptr = e->next[0];
                tree_edges->Delete(e);
                continue;
            }
            e_ptr = &e->next[0];
            
            eps2 = e->head[0]->eps;
            e->pq00.Update(-eps - eps2);
        }
    }
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next)
        r->tree->eps = 0;
}

bool PerfectMatching::UpdateDuals() {
    Node* r;
    
    double start_time = get_time();
    
    ////////////////////////////////////////////////////////////////////////////////////
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        Tree* t = r->tree;
        PriorityQueue<REAL>::Item* q;
        REAL eps = PM_INFTY;
        if ((q = t->pq0.GetMin()))
            eps = q->slack;
        if ((q = t->pq_blossoms.GetMin()) && eps > q->slack)
            eps = q->slack;
        while ((q = t->pq00.GetMin())) {
            if (ProcessEdge00((Edge*) q, false))
                break;
            t->pq00.Remove(q, pq_buf);
        }
        if (q && 2 * eps > q->slack)
            eps = q->slack / 2;
        t->eps_delta = eps - t->eps;
    }
    
    if (tree_num >= options.dual_LP_threshold * node_num) {
        if (options.dual_greedy_update_option == 0)
            ComputeEpsCC();
        else if (options.dual_greedy_update_option == 1)
            ComputeEpsSCC();
        else
            ComputeEpsSingle();
    } else
        ComputeEpsGlobal();
    
    REAL delta = 0;
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        if (r->tree->eps_delta > 0) {
            delta += r->tree->eps_delta;
            r->tree->eps += r->tree->eps_delta;
        }
    }
    
    stat.dual_time += get_time() - start_time;
    
    return (delta > PM_THRESHOLD );
}

//----------//
//  Expand  //
//----------//

inline void PerfectMatching::ProcessSelfloop(Node* b, Edge* a) {
    int dir;
    Node* j;
    Node* prev[2];
    for (dir = 0; dir < 2; dir++) {
        j = a->head[dir];
        GET_PENULTIMATE_BLOSSOM(j);
        prev[dir] = j;
    }
    if (prev[0] != prev[1]) {
        ADD_EDGE(prev[0], a, 1);
        ADD_EDGE(prev[1], a, 0);
        a->slack -= 2 * prev[0]->blossom_eps;
    } else {
        a->next[0] = prev[0]->blossom_selfloops;
        prev[0]->blossom_selfloops = a;
    }
    
}

void PerfectMatching::Expand(Node* b) {
    assert(b->is_blossom);
    assert(b->is_outer);
    assert(b->flag == 1);
    
    double start_time = get_time();
    
    Node* i;
    Node* j;
    Node* k;
    Edge* a;
    EdgeIterator I;
    int dir;
    ExpandTmpItem* tmp_item;
    Tree* t = b->tree;
    REAL eps = t->eps;
    Edge* a_augment = NULL;
    
    GET_TREE_PARENT(b, i);
    a = ARC_TO_EDGE_PTR(b->tree_parent);
    dir = ARC_TO_EDGE_DIR(b->tree_parent);
    
    j = a->head0[1 - dir];
    GET_PENULTIMATE_BLOSSOM(j);
    MOVE_EDGE(b, j, a, dir);
    
    a = ARC_TO_EDGE_PTR(b->match);
    dir = ARC_TO_EDGE_DIR(b->match);
    k = a->head0[1 - dir];
    GET_PENULTIMATE_BLOSSOM(k);
    MOVE_EDGE(b, k, a, dir);
    
    i = ARC_HEAD(k->blossom_sibling) ;
    while (1) {
        tmp_item = expand_tmp_list->New();
        tmp_item->i = i;
        tmp_item->blossom_parent = i->blossom_parent;
        tmp_item->blossom_grandparent = i->blossom_grandparent;
        i->flag = 2;
        
        // blossom_selfloops
        i->is_outer = 1;
        while ((a = i->blossom_selfloops)) {
            i->blossom_selfloops = a->next[0];
            ProcessSelfloop(i, a);
        }
        i->is_outer = 0;
        
        if (i == k)
            break;
        i->match = i->blossom_sibling;
        j = ARC_HEAD(i->match) ;
        tmp_item = expand_tmp_list->New();
        tmp_item->i = j;
        tmp_item->blossom_parent = j->blossom_parent;
        tmp_item->blossom_grandparent = j->blossom_grandparent;
        j->flag = 2;
        
        // blossom_selfloops
        j->is_outer = 1;
        while ((a = j->blossom_selfloops)) {
            j->blossom_selfloops = a->next[0];
            ProcessSelfloop(j, a);
        }
        j->is_outer = 0;
        
        j->match = ARC_REV(i->match);
        i = ARC_HEAD(j->blossom_sibling) ;
    }   
    k->match = b->match;
    i = ARC_TAIL(b->tree_parent) ;
    Arc* aa = i->blossom_sibling;
    i->flag = 1;
    i->tree = b->tree;
    i->y += b->tree->eps;
    i->tree_parent = b->tree_parent;
    if (i != k) {
        Node** i_ptr;
        if (i->match == aa) {
            i = ARC_HEAD(i->match) ;
            i_ptr = &j;
            while ( 1 )
            {   
                aa = i->blossom_sibling;
                i->flag = 0; i->tree = b->tree; i->y -= t->eps;
                *i_ptr = i;
                i_ptr = &i->first_tree_child;
                i->tree_sibling_prev = i;
                i->tree_sibling_next = NULL;
                i = ARC_HEAD(aa);
                i->flag = 1; i->tree = b->tree; i->y += t->eps;
                i->tree_parent = ARC_REV(aa);
                if (i == k) break;
                i = ARC_HEAD(i->match);
            }
            *i_ptr = ARC_HEAD(k->match);
        }
        else
        {   
            i = k;
            j = ARC_HEAD(k->match);
            do
            {   
                i->tree_parent = i->blossom_sibling;
                i->flag = 1; i->tree = b->tree; i->y += b->tree->eps;
                i = ARC_HEAD(i->tree_parent);
                i->flag = 0; i->tree = b->tree; i->y -= b->tree->eps;
                i->first_tree_child = j;
                j = i;
                i->tree_sibling_prev = i;
                i->tree_sibling_next = NULL;
                i = ARC_HEAD(i->match);
            }while ( i->flag != 1 );
        }   
        i = ARC_HEAD(k->match) ;

        j->tree_sibling_prev = i->tree_sibling_prev;
        j->tree_sibling_next = i->tree_sibling_next;
        if (i->tree_sibling_prev->tree_sibling_next)
            i->tree_sibling_prev->tree_sibling_next = j;
        else
            ARC_HEAD(b->tree_parent) ->first_tree_child = j;
        if (i->tree_sibling_next)
            i->tree_sibling_next->tree_sibling_prev = j;
        else
            ARC_HEAD(b->tree_parent) ->first_tree_child->tree_sibling_prev = j;

        i->tree_sibling_prev = i;
        i->tree_sibling_next = NULL;
    }
    
    // go through inner arcs
    i = k;
    while (1) {
        // "-" node
        if (i->is_blossom) {
            a = ARC_TO_EDGE_PTR(i->match);
            REAL tmp = a->slack;
            a->slack = i->y;
            i->y = tmp;
            t->pq_blossoms.Add(a);
        }
        FOR_ALL_EDGES(i, a, dir, I)
        {
            j = a->head[dir];
            if (j->flag != 0)
                a->slack -= eps;
        }
        i->is_processed = 1;
        if (i->tree_parent == b->tree_parent)
            break;
        i = ARC_HEAD(i->tree_parent) ;
        // "+" node
        FOR_ALL_EDGES(i, a, dir, I)
        {
            j = a->head[dir];
            if (j->flag == 2) {
                a->slack += eps;
                t->pq0.Add(a);
            } else if (j->flag == 0 && i < j) {
                a->slack += 2 * eps;
                t->pq00.Add(a);
            }
        }
        i->is_processed = 1;
        i = ARC_HEAD(i->match) ;
    }   

        // go through boundary arcs
    for (tmp_item = expand_tmp_list->ScanFirst(); tmp_item; tmp_item =
            expand_tmp_list->ScanNext()) {
        i = tmp_item->i;
        j = tmp_item->blossom_parent;
        tmp_item->blossom_parent = i->blossom_parent;
        i->blossom_parent = j;
        j = tmp_item->blossom_grandparent;
        tmp_item->blossom_grandparent = i->blossom_grandparent;
        i->blossom_grandparent = j;
    }
    for (dir = 0; dir < 2; dir++) {
        if (!b->first[dir])
            continue;
        b->first[dir]->prev[dir]->next[dir] = NULL;
        
        Edge* a_next;
        for (a = b->first[dir]; a; a = a_next) {
            a_next = a->next[dir];
            i = a->head0[1 - dir];
            GET_PENULTIMATE_BLOSSOM2(i);
            ADD_EDGE(i, a, dir);
            GET_OUTER_HEAD(a, dir, j);
            
            if (i->flag == 1)
                continue;
            
            if (j->flag == 0 && j->tree != t)
                j->tree->pq_current->pq01[1 - j->tree->dir_current].Remove(a, pq_buf);
            
            if (i->flag == 2) {
                a->slack += eps;
                if (j->flag == 0)
                    j->tree->pq0.Add(a);
            } else {
                a->slack += 2 * eps;
                if (j->flag == 2)
                    t->pq0.Add(a);
                else if (j->flag == 0) {
                    if (j->tree != t) {
                        if (!j->tree->pq_current)
                            AddTreeEdge(t, j->tree);
                        if (a->slack <= j->tree->eps + eps)
                            a_augment = a;
                    }
                    j->tree->pq_current->pq00.Add(a);
                } else if (j->tree != t) {
                    if (!j->tree->pq_current)
                        AddTreeEdge(t, j->tree);
                    j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
                }
                
            }
        }
    }
    for (tmp_item = expand_tmp_list->ScanFirst(); tmp_item; tmp_item =
            expand_tmp_list->ScanNext()) {
        i = tmp_item->i;
        i->blossom_parent = tmp_item->blossom_parent;
        i->blossom_grandparent = tmp_item->blossom_grandparent;
        i->is_outer = 1;
    }
    expand_tmp_list->Reset();
    
    b->tree_sibling_next = removed_first;
    removed_first = b;
    removed_num++;
    if (4 * removed_num > node_num)
        FreeRemoved();
    
    blossom_num--;
    stat.expand_count++;
    
    stat.expand_time += get_time() - start_time;
    
    if (a_augment)
        Augment(a_augment);
}

void PerfectMatching::FreeRemoved() {
    Node* i0;
    Node* i;
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        for (i = i0; !i->is_outer && !i->is_marked; i = i->blossom_parent) {
            i->is_marked = 1;
            if (i->blossom_grandparent->is_removed)
                i->blossom_grandparent = i->blossom_parent;
        }
    }
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        for (i = i0; !i->is_outer && i->is_marked; i = i->blossom_parent) {
            i->is_marked = 0;
        }
    }
    
    while ((i = removed_first)) {
        removed_first = i->tree_sibling_next;
        blossoms->Delete(i);
        removed_num--;
    }

    assert(removed_num == 0);
}

//----------//
//  Init  //
//----------//

void PerfectMatching::InitGreedy(bool allocate_trees) {
    Node* i;
    int dir;
    Edge* a;
    EdgeIterator I;
    Tree* t = NULL;
    Node* last_root = &nodes[node_num];
    REAL slack_min;
    
    for (i = nodes; i < nodes + node_num; i++)
        i->y = PM_INFTY;
    for (a = edges; a < edges + edge_num; a++) {
        if (a->head[0]->y > a->slack)
            a->head[0]->y = a->slack;
        if (a->head[1]->y > a->slack)
            a->head[1]->y = a->slack;
    }
    for (a = edges; a < edges + edge_num; a++) {
        i = a->head[0];
        if (!i->is_outer) {
            i->is_outer = 1;
            i->y /= 2;
        }
        a->slack -= i->y;
        i = a->head[1];
        if (!i->is_outer) {
            i->is_outer = 1;
            i->y /= 2;
        }
        a->slack -= i->y;
    }
    
    tree_num = node_num;
    for (i = nodes; i < nodes + node_num; i++) {
        if (i->flag == 2)
            continue;
        slack_min = PM_INFTY;
        FOR_ALL_EDGES(i, a, dir, I)
            if (slack_min > a->slack)
                slack_min = a->slack;
        i->y += slack_min;
        FOR_ALL_EDGES(i, a, dir, I)
        {
            if (a->slack <= slack_min && i->flag == 0 && a->head[dir]->flag == 0) {
                i->flag = 2;
                a->head[dir]->flag = 2;
                i->match = EDGE_DIR_TO_ARC(a, dir);
                a->head[dir]->match = EDGE_DIR_TO_ARC(a, 1-dir);
                tree_num -= 2;
            }
            a->slack -= slack_min;
        }
    }
    if (allocate_trees) {
        if (tree_num > tree_num_max) {
            if (trees)
                free(trees);
            tree_num_max = tree_num;
            trees = (Tree*) malloc(tree_num_max * sizeof(Tree));
        }
        t = trees;
    }
    for (i = nodes; i < nodes + node_num; i++) {
        if (i->flag != 0)
            continue;
        i->is_tree_root = 1;
        i->first_tree_child = NULL;
        i->tree_sibling_prev = last_root;
        last_root->tree_sibling_next = i;
        last_root = i;
        if (allocate_trees) {
            i->tree = t;
            t->root = i;
            t->eps = 0;
            t->first[0] = t->first[1] = NULL;
            t->pq_current = NULL;
            t->pq00.Reset();
            t->pq0.Reset();
            t->pq_blossoms.Reset();
            t++;
        }
    }
    last_root->tree_sibling_next = NULL;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

PerfectMatching::Node* PerfectMatching::FindBlossomRootInit(Edge* a0) {
    Node* i;
    Node* j;
    Node* _i[2];
    Node* r;
    int branch;
    
    _i[0] = ARC_HEAD(a0) ;
    _i[1] = ARC_TAIL(a0) ;
    branch = 0;
    while (1) {
        if (!_i[branch]->is_outer) {
            r = _i[branch];
            j = _i[1 - branch];
            break;
        }
        _i[branch]->is_outer = 0;
        if (_i[branch]->is_tree_root) {
            j = _i[branch];
            i = _i[1 - branch];
            while (i->is_outer) {
                i->is_outer = 0;
                i = ARC_HEAD(i->match) ;
                i->is_outer = 0;
                i = ARC_HEAD(i->tree_parent) ;
            }   
            r = i;
            break;
        }
        i = ARC_HEAD(_i[branch]->match) ;
        i->is_outer = 0;
        _i[branch] = ARC_HEAD(i->tree_parent) ;
        branch = 1 - branch;
    }
    i = r;
    while (i != j) {
        i = ARC_HEAD(i->match) ;
        i->is_outer = 1;
        i = ARC_HEAD(i->tree_parent);
        i->is_outer = 1;
    }   
    return r;
}

void PerfectMatching::ShrinkInit(Edge* a0, Node* tree_root) {
    int branch, flag;
    Node* i;
    Node* j;
    Node* r;
    Arc* a_prev;
    Arc* aa;
    
    tree_root->flag = 2;
    i = tree_root->first_tree_child;
    if (i)
        while (1) {
            ARC_HEAD(i->match) ->flag = 2;
            i->flag = 2;

            MOVE_NODE_IN_TREE(i);
        }

    r = FindBlossomRootInit(a0);
    
    if (!r->is_tree_root) {
        j = ARC_HEAD(r->match) ;
        j->match = aa = j->tree_parent;
        i = ARC_HEAD(aa);
        while ( !i->is_tree_root )
        {   
            j = ARC_HEAD(i->match);
            i->match = ARC_REV(aa);
            j->match = aa = j->tree_parent;
            i = ARC_HEAD(aa);
        }
        i->match = ARC_REV(aa);
    }   

    tree_root->is_tree_root = 0;
    
    branch = 0;
    flag = 0;
    a_prev = EDGE_DIR_TO_ARC(a0, 0);
    i = ARC_HEAD(a_prev) ;
    while (1) {
        Arc* a_next = (flag == 0) ? i->match : i->tree_parent;
        flag = 1 - flag;
        i->flag = 0;
        i->match = NULL;
        if (branch == 0) {
            i->blossom_sibling = a_next;
            if (i == r) {
                branch = 1;
                flag = 0;
                a_prev = ARC_REV(a0);
                i = ARC_HEAD(a_prev) ;
                if (i == r)
                    break;
            } else {
                a_prev = i->blossom_sibling;
                i = ARC_HEAD(a_prev) ;
            }
        }
        else
        {   
            i->blossom_sibling = ARC_REV(a_prev);
            a_prev = a_next;
            i = ARC_HEAD(a_prev);
            if (i == r) break;
        }
    }           
    i->blossom_sibling = ARC_REV(a_prev);
}

void PerfectMatching::ExpandInit(Node* k) {
    Node* i = ARC_HEAD(k->blossom_sibling) ;
    Node* j;

    while ( 1 )
    {   
        i->flag = 2; i->is_outer = 1;
        if (i == k) break;
        i->match = i->blossom_sibling;
        j = ARC_HEAD(i->match);
        j->flag = 2; j->is_outer = 1;
        j->match = ARC_REV(i->match);
        i = ARC_HEAD(j->blossom_sibling);
    }
}   

void PerfectMatching::AugmentBranchInit(Node* i0, Node* r) {
    Node* tree_root_prev = r->tree_sibling_prev;
    Node* i;
    Node* j;
    Arc* aa;
    
    r->flag = 2;
    i = r->first_tree_child;
    if (i)
        while (1) {
            ARC_HEAD(i->match) ->flag = 2;
            i->flag = 2;

            MOVE_NODE_IN_TREE(i);
        }
    i = i0;
    if (!i0->is_tree_root) {
        j = ARC_HEAD(i0->match) ;
        j->match = aa = j->tree_parent;
        i = ARC_HEAD(aa);
        while ( !i->is_tree_root )
        {   
            j = ARC_HEAD(i->match);
            i->match = ARC_REV(aa);
            j->match = aa = j->tree_parent;
            i = ARC_HEAD(aa);
        }
        i->match = ARC_REV(aa);
    }   
    r->is_tree_root = 0;
    tree_root_prev->tree_sibling_next = r->tree_sibling_next;
    if (r->tree_sibling_next)
        r->tree_sibling_next->tree_sibling_prev = tree_root_prev;
    tree_num--;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

// true_slack(a) = slack(a) + ...

// i->flag=0, i->is_processed=1:              true_slack -= eps
// i->flag=1, i->match->head->is_processed=1: true_slack += eps - slack(i->match)

void PerfectMatching::InitGlobal() {
    Node* i;
    Node* j;
    Node* r;
    Node* r2;
    Node* r3 = NULL;    // initialize to prevent compiler warning
    Edge* a;
    EdgeIterator I;
    int dir;
    Tree TREE;
    enum {
        NONE, AUGMENT, SHRINK
    } flag;
    
    InitGreedy();
    
    for (i = nodes; i < nodes + node_num; i++)
        i->best_edge = NULL;
    
    PriorityQueue<REAL> pq;
    
    for (r = nodes[node_num].tree_sibling_next; r;) {
        r2 = r->tree_sibling_next;
        if (r2)
            r3 = r2->tree_sibling_next;
        i = r;
        
        pq.Reset();
        
        r->tree = &TREE;
        
        REAL eps = 0;
        Arc* critical_arc = NULL;
        REAL critical_eps = PM_INFTY;
        flag = NONE;
        Node* branch_root = i;
        
        while (1) {
            i->is_processed = 1;
            i->y -= eps;
            if (!i->is_tree_root)
                ARC_HEAD(i->match) ->y += eps;

            FOR_ALL_EDGES(i, a, dir, I)
            {
                a->slack += eps;
                j = a->head[dir];
                
                if (j->tree == &TREE) {
                    // same tree
                    if (j->flag == 0) {
                        REAL slack = a->slack;
                        if (!j->is_processed)
                            slack += eps;
                        if (2 * critical_eps > slack || critical_arc == NULL) {
                            flag = SHRINK;
                            critical_eps = slack / 2;
                            critical_arc = EDGE_DIR_TO_ARC(a, dir);
                            if (critical_eps <= eps)
                                break;
                            //pq.DecreaseUpperBound(critical_eps);
                        }
                    }
                } else if (j->flag == 0) {
                    // different tree
                    if (critical_eps >= a->slack || critical_arc == NULL) {
                        flag = AUGMENT;
                        critical_eps = a->slack;
                        critical_arc = EDGE_DIR_TO_ARC(a, dir);
                        if (critical_eps <= eps)
                            break;
                        //pq.DecreaseUpperBound(critical_eps);
                    }
                } else {
                    // free node
                    if (a->slack > eps) {
                        if (a->slack < critical_eps) {
                            if (j->best_edge == NULL) {
                                j->best_edge = a;
                                pq.Add(a);
                            } else {
                                if (a->slack < j->best_edge->slack) {
                                    pq.Decrease(j->best_edge, a, pq_buf);
                                    j->best_edge = a;
                                }
                            }
                        }
                    } else {
                        assert(j->flag == 2 && !j->is_blossom && !ARC_HEAD(j->match)->is_blossom);
                        if (j->best_edge)
                            pq.Remove(j->best_edge, pq_buf);
                        j->flag = 1;
                        j->tree = i->tree;
                        j->tree_parent = EDGE_DIR_TO_ARC(a, 1-dir);
                        j = ARC_HEAD(j->match) ;
                        if (j->best_edge)
                            pq.Remove(j->best_edge, pq_buf);
                        ADD_TREE_CHILD(i, j);
                    }
                }
            }
            
            if (dir < 2 && a) {
                Edge* atmp = a;
                int dirtmp = dir;
                CONTINUE_FOR_ALL_EDGES(i, atmp, dirtmp, I)
                    atmp->slack += eps;
                break;
            }
            
            // move i
            if (i->first_tree_child)
                i = i->first_tree_child;
            else {
                while (i != branch_root && !i->tree_sibling_next) {
                    i = ARC_HEAD(i->match) ; i = ARC_HEAD(i->tree_parent);}
                if (i == branch_root)
                {   
                    PriorityQueue<REAL>::Item* q = pq.GetMin();
                    if (q == NULL || q->slack >= critical_eps)
                    {   
                        eps = critical_eps;
                        break;
                    }
                    pq.Remove(q, pq_buf);
                    a = (Edge*)q;
                    dir = (a->head[0]->flag == 2) ? 0 : 1;
                    j = a->head[0];
                    Arc* aa = EDGE_DIR_TO_ARC(a, dir);
                    eps = a->slack;
                    assert(eps < critical_eps);

                    // continue growth
                    i = ARC_TAIL(aa);
                    j = ARC_HEAD(aa);

                    assert(j->flag == 2 && !j->is_blossom && !ARC_HEAD(j->match)->is_blossom);
                    j->flag = 1;
                    j->tree = i->tree;
                    j->tree_parent = ARC_REV(aa);
                    j = ARC_HEAD(j->match);
                    if (j->best_edge) pq.Remove(j->best_edge, pq_buf);
                    ADD_TREE_CHILD(i, j);
                    i = branch_root = j;
                    continue;
                }
                i = i->tree_sibling_next;
            }
        }           

                    // update slacks
        i = r;
        while (1) {
            if (i->is_processed) {
                i->y += eps;
                if (!i->is_tree_root) {
                    j = ARC_HEAD(i->match) ;
                    j->y -= eps;
                    REAL delta = eps - ARC_TO_EDGE_PTR(i->match)->slack;
                    FOR_ALL_EDGES(j, a, dir, I) a->slack += delta;
                    j->best_edge = NULL;
                }   
                FOR_ALL_EDGES(i, a, dir, I)
                {
                    if (!PriorityQueue<REAL>::isReset(a)) {
                        assert(a->head[dir]->flag == 2 && a->head[dir]->best_edge == a);
                        a->head[dir]->best_edge = NULL;
                        PriorityQueue<REAL>::ResetItem(a);
                    }
                    a->slack -= eps;
                }
                
                i->is_processed = 0;
            } else {
                if (!i->is_tree_root)
                    ARC_HEAD(i->match) ->best_edge = NULL;
                }
                i->best_edge = NULL;

                MOVE_NODE_IN_TREE(i);
            }

        i = ARC_TAIL(critical_arc) ;
        j = ARC_HEAD(critical_arc) ;
        if (flag == SHRINK) {
            // shrink
            ShrinkInit(ARC_TO_EDGE_PTR(critical_arc), r);
        } else {
            // augment
            AugmentBranchInit(i, r);
            if (j->is_outer) {
                AugmentBranchInit(j, j);
            } else {
                ExpandInit(j);
                tree_num--;
            }
            i->match = critical_arc;
            j->match = ARC_REV(critical_arc);
        }
        
        r = r2;
        if (r && !r->is_tree_root)
            r = r3;
    }
    
    if (tree_num > tree_num_max) {
        if (trees)
            free(trees);
        tree_num_max = tree_num;
        trees = (Tree*) malloc(tree_num_max * sizeof(Tree));
    }
    Tree* t = trees;
    for (r = nodes; r < nodes + node_num; r++) {
        if (!r->is_outer) {
            ExpandInit(r);
            r->is_tree_root = 1;
            r->flag = 0;
            r->first_tree_child = NULL;
            if (t == trees) {
                nodes[node_num].tree_sibling_next = r;
                r->tree_sibling_prev = &nodes[node_num];
            } else {
                (t - 1)->root->tree_sibling_next = r;
                r->tree_sibling_prev = (t - 1)->root;
            }
            r->tree = t;
            t->root = r;
            t->eps = 0;
            t->first[0] = t->first[1] = NULL;
            t->pq_current = NULL;
            t->pq00.Reset();
            t->pq0.Reset();
            t->pq_blossoms.Reset();
            t++;
        }
    }
    assert(t == trees+tree_num);
    if (t == trees)
        nodes[node_num].tree_sibling_next = NULL;
    else
        (t - 1)->root->tree_sibling_next = NULL;
}

//----------//
//  Interface  //
//----------//

PerfectMatching::PerfectMatching(int nodeNum, int edgeNumMax)
        : trees(NULL), node_num(nodeNum), edge_num(0), edge_num_max(edgeNumMax), tree_num_max(0), removed_first(
                  NULL), blossom_num(0), removed_num(0), first_solve(true) {
    if (node_num & 1) {
        printf("# of nodes is odd: perfect matching cannot exist\n");
        exit(1);
    }
    nodes = (Node*) malloc((node_num + 1) * sizeof(Node));
    edges_orig = (char*) malloc(edge_num_max * sizeof(Edge) + 1);
    edges = (Edge*) ((((POINTER_TYPE) edges_orig) & 1) ? (edges_orig + 1) : edges_orig);
    memset(nodes, 0, (node_num + 1) * sizeof(Node));
    
    blossoms = new DBlock<Node>(256);
    tree_edges = new DBlock<TreeEdge>(256);
    expand_tmp_list = new Block<ExpandTmpItem>(256);
    pq_buf = PriorityQueue<REAL>::AllocateBuf();
}

void PerfectMatching::Save(char* filename, int format) {
    if (!first_solve) {
        printf("Save() cannot be called after Solve()!\n");
        exit(1);
    }
    int e;
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        printf("Can't open %s\n", filename);
        exit(1);
    }
    if (format == 0) {
        fprintf(fp, "p edge %d %d\n", node_num, edge_num);
        for (e = 0; e < edge_num; e++) {
            fprintf(fp, "e %d %d %d\n", 1 + (int) (edges[e].head0[1] - nodes),
                    1 + (int) (edges[e].head0[0] - nodes), (int) edges[e].slack / COST_FACTOR);
        }
    } else {
        fprintf(fp, "%d %d\n", node_num, edge_num);
        for (e = 0; e < edge_num; e++) {
            fprintf(fp, "%d %d %d\n", (int) (edges[e].head0[1] - nodes),
                    (int) (edges[e].head0[0] - nodes), (int) edges[e].slack / COST_FACTOR);
        }
    }
    fclose(fp);
}

PerfectMatching::~PerfectMatching() {
    free(nodes);
    free(edges_orig);
    delete blossoms;
    delete tree_edges;
    delete expand_tmp_list;
    if (trees)
        free(trees);
    PriorityQueue<REAL>::DeallocateBuf(pq_buf);
}

PerfectMatching::EdgeId PerfectMatching::AddEdge(NodeId _i, NodeId _j, REAL cost) {
    if (_i < 0 || _i >= node_num || _j < 0 || _j > node_num || _i == _j) {
        printf("wrong node id's! (%d,%d)\n", _i, _j);
        exit(1);
    }
    if (edge_num >= edge_num_max)
        ReallocateEdges();
    Node* i = nodes + _i;
    Node* j = nodes + _j;
    Edge* a = edges + edge_num;
    
    ADD_EDGE(i, a, 0);
    ADD_EDGE(j, a, 1);
    a->head0[0] = j;
    a->head0[1] = i;
    
    a->slack = cost * COST_FACTOR;
    PriorityQueue<REAL>::ResetItem(a);
    
    return edge_num++;
}

int PerfectMatching::GetSolution(EdgeId e) {
    assert(e>=0 && e<edge_num);
    Edge* a = edges + e;
    return (a->head0[1]->match == EDGE_DIR_TO_ARC(a, 0) ) ? 1 : 0;
}

PerfectMatching::NodeId PerfectMatching::GetMatch(NodeId i) {
    assert(i>=0 && i<node_num);
    return (int) (ARC_HEAD0(nodes[i].match) -nodes);
}   

void PerfectMatching::GetRealEndpoints(Edge* a, Node*& tail, Node*& head) {
    Node* i;
    Node* j;
    int delta = 0;
    
    for (i = a->head0[1]; !i->is_outer; i = i->blossom_parent, delta--) {
    }
    for (j = a->head0[0]; !j->is_outer; j = j->blossom_parent, delta++) {
    }
    if (i == j) {
        i = a->head0[1];
        j = a->head0[0];
        while (delta < 0) {
            i = i->blossom_parent;
            delta++;
        }
        while (delta > 0) {
            j = j->blossom_parent;
            delta--;
        }
        while (i->blossom_parent != j->blossom_parent) {
            i = i->blossom_parent;
            j = j->blossom_parent;
        }
    }
    tail = i;
    head = j;
    assert(
            (i->is_outer && j->is_outer) || (i->blossom_parent==j->blossom_parent && !i->is_outer && !j->is_outer));
}

void PerfectMatching::ReallocateEdges() {
    printf(
            "Warning: reallocating edges. Increasing edge_num_max in the constructor may improve memory efficiency!\n");
    edge_num_max = edge_num_max * 3 / 2 + 16;
    char* edges_orig_old = edges_orig;
    Edge* edges_old = edges;
    edges_orig = (char*) realloc(edges_orig_old, edge_num_max * sizeof(Edge) + 1);
    edges = (Edge*) ((((POINTER_TYPE) edges_orig_old) & 1) ? (edges_orig + 1) : edges_orig);
    if (((POINTER_TYPE) edges) & 1) {
        char* edges_orig_old2 = edges_orig;
        Edge* edges_old2 = edges;
        
        edges_orig = (char*) malloc(edge_num_max * sizeof(Edge) + 1);
        edges = (Edge*) ((((POINTER_TYPE) edges_orig_old) & 1) ? (edges_orig + 1) : edges_orig);
        memcpy(edges, edges_old2, edge_num * sizeof(Edge));
        free(edges_orig_old2);
    }
    
#define UPDATE_EDGE_PTR(ptr) ptr = (Edge*)((char*)(ptr) + ((char*)edges - (char*)edges_old))
#define UPDATE_ARC_PTR(ptr) ptr = (Arc*)((char*)(ptr) + ((char*)edges - (char*)edges_old))
    
    Node* i;
    Edge* a;
    for (a = edges; a < edges + edge_num; a++) {
        if (a->next[0])
            UPDATE_EDGE_PTR(a->next[0]);
        if (a->next[1])
            UPDATE_EDGE_PTR(a->next[1]);
        if (a->prev[0])
            UPDATE_EDGE_PTR(a->prev[0]);
        if (a->prev[1])
            UPDATE_EDGE_PTR(a->prev[1]);
    }
    if (first_solve) {
        for (i = nodes; i < nodes + node_num; i++) {
            if (i->first[0])
                UPDATE_EDGE_PTR(i->first[0]);
            if (i->first[1])
                UPDATE_EDGE_PTR(i->first[1]);
        }
    } else {
        Node* i0;
        for (i0 = nodes; i0 < nodes + node_num; i0++) {
            i = i0;
            while (1) {
                if (i->is_outer) {
                    UPDATE_ARC_PTR(i->match);
                    if (i->first[0])
                        UPDATE_EDGE_PTR(i->first[0]);
                    if (i->first[1])
                        UPDATE_EDGE_PTR(i->first[1]);
                    break;
                }
                UPDATE_ARC_PTR(i->blossom_sibling);
                if (i->first[0])
                    UPDATE_EDGE_PTR(i->first[0]);
                if (i->first[1])
                    UPDATE_EDGE_PTR(i->first[1]);
                
                i = i->blossom_parent;
                if (i->is_outer) {
                    if (i->is_marked)
                        break;
                    i->is_marked = 1;
                } else {
                    if (!i->is_marked)
                        break;
                    i->is_marked = 0;
                }
            }
        }
        for (i0 = nodes; i0 < nodes + node_num; i0++) {
            i = i0;
            while (1) {
                if (i->is_outer)
                    break;
                
                i = i->blossom_parent;
                if (i->is_outer) {
                    if (!i->is_marked)
                        break;
                    i->is_marked = 0;
                } else {
                    if (i->is_marked)
                        break;
                    i->is_marked = 1;
                }
            }
        }
    }
}

int PerfectMatching::GetBlossomNum() {
    return blossom_num;
}

void PerfectMatching::GetDualSolution(int* blossom_parents, REAL* twice_y) {
    int _i0, id = node_num;
    int* child_ptr;
    Node* i0;
    Node* i;
    int* tmp_array = new int[blossom_num];
    
    int* tmp_array_ptr = tmp_array;
    for (_i0 = 0, i0 = nodes; _i0 < node_num; _i0++, i0++) {
        twice_y[_i0] = i0->y;
        if (i0->is_outer) {
            blossom_parents[_i0] = -1;
            continue;
        }
        child_ptr = &blossom_parents[_i0];
        i = i0->blossom_parent;
        while (1) {
            if (i->is_marked) {
                *child_ptr = i->lca_preorder;
                break;
            }
            i->is_marked = 1;
            *tmp_array_ptr++ = i->lca_preorder;
            *child_ptr = i->lca_preorder = id++;
            child_ptr = &blossom_parents[i->lca_preorder];
            twice_y[i->lca_preorder] = i->y;
            if (i->is_outer) {
                *child_ptr = -1;
                break;
            }
            i = i->blossom_parent;
        }
    }
    
    assert(id == node_num+blossom_num && tmp_array_ptr == tmp_array + blossom_num);
    
    tmp_array_ptr = tmp_array;
    for (_i0 = 0, i0 = nodes; _i0 < node_num; _i0++, i0++) {
        if (i0->is_outer)
            continue;
        i = i0->blossom_parent;
        while (1) {
            if (!i->is_marked)
                break;
            i->is_marked = 0;
            i->lca_preorder = *tmp_array_ptr++;
            if (i->is_outer)
                break;
            i = i->blossom_parent;
        }
    }
    
    delete[] tmp_array;
}

//----------//
//  Main  //
//----------//

void PerfectMatching::Finish() {
    
#define IS_VALID_MATCH(i) ((Edge*)(i->match) >= edges && (Edge*)(i->match) < edges + edge_num)
    
    Node* i0;
    Node* i;
    Node* j;
    Node* k;
    Node* b;
    Node* b_prev;
    Node* b_prev_prev;
    
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        if (IS_VALID_MATCH(i0))
            continue;
        b_prev = NULL;
        b = i0;
        do {
            b->blossom_grandparent = b_prev;
            b_prev = b;
            b = b->blossom_parent;
        } while (!IS_VALID_MATCH(b));

        b_prev_prev = b_prev->blossom_grandparent;
        while (1) {
            for (k = ARC_TAIL0(b->match) ; k->blossom_parent!=b; k=k->blossom_parent) {}
            k->match = b->match;
            i = ARC_HEAD(k->blossom_sibling);
            while ( i != k )
            {   
                i->match = i->blossom_sibling;
                j = ARC_HEAD(i->match);
                j->match = ARC_REV(i->match);
                i = ARC_HEAD(j->blossom_sibling);
            }

            b = b_prev;
            if (!b->is_blossom) break;
            b_prev = b_prev_prev;
            b_prev_prev = b_prev->blossom_grandparent;
        }
    }
}           

void PerfectMatching::AddTreeEdge(Tree* t0, Tree* t1) {
    TreeEdge* e = tree_edges->New();
    e->head[0] = t1;
    e->head[1] = t0;
    e->next[0] = t0->first[0];
    t0->first[0] = e;
    e->next[1] = t1->first[1];
    t1->first[1] = e;
    
    e->pq00.Reset();
    e->pq01[0].Reset();
    e->pq01[1].Reset();
    
    t1->pq_current = e;
    t1->dir_current = 0;
}

bool PerfectMatching::ProcessEdge00(Edge* a, bool update_boundary_edge) {
    int dir;
    Node* j;
    Node* prev[2];
    Node* last[2];
    for (dir = 0; dir < 2; dir++) {
        if (a->head[dir]->is_outer) {
            prev[dir] = NULL;
            last[dir] = a->head[dir];
        } else {
            j = a->head[dir];
            GET_PENULTIMATE_BLOSSOM(j);
            prev[dir] = j;
            last[dir] = prev[dir]->blossom_parent;
            //assert(last[dir]->is_outer);
        }
    }
    
    if (last[0] != last[1]) {
        for (dir = 0; dir < 2; dir++) {
            j = a->head[dir];
            if (j != last[dir]) {
                int dir_rev = 1 - dir;
                MOVE_EDGE(j, last[dir], a, dir_rev);
            }
        }
        if (update_boundary_edge)
            a->slack -= 2 * a->head[0]->tree->eps;
        return true;
    }
    
    if (prev[0] != prev[1]) {
        for (dir = 0; dir < 2; dir++) {
            j = a->head[dir];
            if (j != prev[dir]) {
                int dir_rev = 1 - dir;
                MOVE_EDGE(j, prev[dir], a, dir_rev);
            }
        }
        a->slack -= 2 * prev[0]->blossom_eps;
        return false;
    }
    
    for (dir = 0; dir < 2; dir++) {
        j = a->head[1 - dir];
        REMOVE_EDGE(j, a, dir);
    }
    a->next[0] = prev[0]->blossom_selfloops;
    prev[0]->blossom_selfloops = a;
    return false;
}

inline void PerfectMatching::AugmentBranch(Node* i0) {
    int dir;
    Tree* t = i0->tree;
    Node* r = t->root;
    Node* tree_root_prev = r->tree_sibling_prev;
    Node* i;
    Node* j;
    Edge* a;
    EdgeIterator I;
    Arc* aa;
    REAL eps = t->eps;
    PriorityQueue<REAL>::Item* q;
    TreeEdge* e;
    TreeEdgeIterator T;
    Tree* t2;
    
    t = r->tree;
    t->pq_current = t;
    
    FOR_ALL_TREE_EDGES_X(t, e, dir, T)
    {
        t2 = e->head[dir];
        e->head[1 - dir] = NULL;    // mark it for deletion
                
        t2->pq_current = e;
        t2->dir_current = dir;
    }
    
    i = r->first_tree_child;
    if (i)
        while (1) {
            Node* i0 = i;
            i = ARC_HEAD(i->match) ;
            if (i->is_processed) {
                if (i->is_blossom) {
                    a = ARC_TO_EDGE_PTR(i->match);
                    REAL tmp = a->slack;
                    a->slack = i->y;
                    i->y = tmp;
                    PriorityQueue<REAL>::ResetItem(a);
                }
                FOR_ALL_EDGES(i, a, dir, I)
                {
                    GET_OUTER_HEAD(a, dir, j);
                    
                    if (j->flag == 0 && j->is_processed) {
                        if (j->tree != t) {
                            a->slack += eps;
                            if (PriorityQueue<REAL>::isReset(a))
                                j->tree->pq0.Add(a);
                        }
                    } else
                        a->slack += eps;
                }
            }
            
            i = i0;
            MOVE_NODE_IN_TREE(i);
        }
    
    ///////////////////////////////////////////////////////////////////
    
    FOR_ALL_TREE_EDGES(t, e, dir)
    {
        t2 = e->head[dir];
        t2->pq_current = NULL;
        
        e->pq01[1 - dir].Merge(t2->pq0);
        for (q = e->pq00.GetFirst(); q; q = e->pq00.GetNext(q)) {
            q->slack -= eps;
            int dir2;
            for (dir2 = 0; dir2 < 2; dir2++)
                GET_OUTER_HEAD((Edge*)q, dir2, j);
        }
        e->pq00.Merge(t2->pq0);
        for (q = e->pq01[dir].GetAndResetFirst(); q; q = e->pq01[dir].GetAndResetNext()) {
            q->slack -= eps;
            int dir2;
            for (dir2 = 0; dir2 < 2; dir2++)
                GET_OUTER_HEAD((Edge*)q, dir2, j);
        }
    }
    for (q = t->pq0.GetAndResetFirst(); q; q = t->pq0.GetAndResetNext()) {
        q->slack -= eps;
        int dir2;
        for (dir2 = 0; dir2 < 2; dir2++)
            GET_OUTER_HEAD((Edge*)q, dir2, j);
    }
    for (q = t->pq00.GetAndResetFirst(); q; q = t->pq00.GetAndResetNext()) {
        ProcessEdge00((Edge*) q);
    }
    
    ///////////////////////////////////////////////////////////////////
    
    r->flag = 2;
    r->is_processed = 0;
    i = r->first_tree_child;
    r->y += eps;
    if (i)
        while (1) {
            j = ARC_HEAD(i->match) ;
            j->flag = 2;
            i->flag = 2;
            j->is_processed = 0;
            i->is_processed = 0;
            j->y -= eps;
            i->y += eps;

            MOVE_NODE_IN_TREE(i);
        }   

            ///////////////////////////////////////////////////////////////////
            
    i = i0;
    if (!i0->is_tree_root) {
        j = ARC_HEAD(i0->match) ;
        GET_TREE_PARENT(j, i);
        j->match = aa = j->tree_parent;
        while ( !i->is_tree_root )
        {   
            j = ARC_HEAD(i->match);
            i->match = ARC_REV(aa);
            GET_TREE_PARENT(j, i);
            j->match = aa = j->tree_parent;
        }
        i->match = ARC_REV(aa);
    }   
    r->is_tree_root = 0;
    tree_root_prev->tree_sibling_next = r->tree_sibling_next;
    if (r->tree_sibling_next)
        r->tree_sibling_next->tree_sibling_prev = tree_root_prev;
    tree_num--;
}

void PerfectMatching::Augment(Edge* a) {
    Node* j;
    int dir;
    
    for (dir = 0; dir < 2; dir++) {
        GET_OUTER_HEAD(a, dir, j);
        AugmentBranch(j);
        j->match = EDGE_DIR_TO_ARC(a, 1-dir);
    }
    if (options.verbose) {
        int k = 1;
        while (k < tree_num)
            k *= 2;
        if (k == tree_num || tree_num <= 8 || (tree_num <= 64 && (tree_num % 8) == 0)) {
            printf("%d.", tree_num);
            fflush(stdout);
        }
    }
}

inline void PerfectMatching::GrowNode(Node* i) {
    //assert(i->is_outer);
    //assert(i->flag == 0);
    
    Edge* a;
    EdgeIterator I;
    int dir;
    Node* j;
    Tree* t = i->tree;
    REAL eps = t->eps;
    Edge* a_augment = NULL;
    
    FOR_ALL_EDGES(i, a, dir, I)
    {
        GET_OUTER_HEAD(a, dir, j);
        
        if (j->flag == 2) {
            a->slack += eps;
            if (a->slack > 0) {
                t->pq0.Add(a);
            } else {
                j->flag = 1;
                j->tree = i->tree;
                j->tree_parent = EDGE_DIR_TO_ARC(a, 1-dir);
                j->y += eps;
                j = ARC_HEAD(j->match) ;
                j->y -= eps;
                ADD_TREE_CHILD(i, j);
            }
        } else {
            if (j->flag == 0 && j->is_processed) {
                if (!PriorityQueue<REAL>::isReset(a))
                    j->tree->pq0.Remove(a, pq_buf);
                if (a->slack <= j->tree->eps && j->tree != t)
                    a_augment = a;
                a->slack += eps;
                if (!j->tree->pq_current)
                    AddTreeEdge(t, j->tree);
                j->tree->pq_current->pq00.Add(a);
            } else {
                a->slack += eps;
                if (j->flag == 1 && j->tree != t) {
                    if (!j->tree->pq_current)
                        AddTreeEdge(t, j->tree);
                    j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
                }
            }
        }
    }
    
    //assert(!i->is_processed);
    i->is_processed = 1;
    
    if (!i->is_tree_root) {
        j = ARC_HEAD(i->match) ;
        //assert(!j->is_processed);
        j->is_processed = 1;
        if (j->is_blossom)
        {   
            a = ARC_TO_EDGE_PTR(i->match);
            REAL tmp = a->slack; a->slack = j->y; j->y = tmp;
            t->pq_blossoms.Add(a);
        }
    }   

    if (a_augment)
        Augment(a_augment);
    
    stat.grow_count++;
}

void PerfectMatching::GrowTree(Node* r, bool new_subtree) {
    //assert(r->flag == 0);
    
    Node* i = r;
    Node* j;
    Node* stop = r->tree_sibling_next;
    if (new_subtree && r->first_tree_child)
        stop = r->first_tree_child;
    Edge* a;
    EdgeIterator I;
    int dir;
    Tree* t = r->tree;
    REAL eps = t->eps;
    int tree_num0 = tree_num;
    
    while (1) {
        if (!i->is_tree_root) {
            // process "-" node
            i = ARC_HEAD(i->match) ;
            FOR_ALL_EDGES(i, a, dir, I)
            {   
                GET_OUTER_HEAD(a, dir, j);

                if (j->flag == 2) a->slack -= eps;
                else
                {   
                    if (j->flag == 0 && j->is_processed)
                    {   
                        if (!PriorityQueue<REAL>::isReset(a)) j->tree->pq0.Remove(a, pq_buf);
                        a->slack -= eps;
                        if (j->tree != t)
                        {   
                            if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
                            j->tree->pq_current->pq01[1-j->tree->dir_current].Add(a);
                        }
                    }
                    else a->slack -= eps;
                }
            }
            i = ARC_HEAD(i->match);
        }
        // process "+" node
        GrowNode(i);
        if (tree_num != tree_num0) break;

        if (i->first_tree_child) i = i->first_tree_child;
        else
        {   
            while (i != r && !i->tree_sibling_next) {i = ARC_HEAD(i->match); GET_TREE_PARENT(i, i);}
            i = i->tree_sibling_next;
        }
        if (i == stop) break;
    }
}           

void PerfectMatching::Solve(bool finish) {
    Node* i;
    Node* j;
    Node* r;
    Node* r2;
    Node* r3 = NULL;    // initialize to prevent compiler warning
    PriorityQueue<REAL>::Item* q;
    Edge* a;
    Tree* t;
    Tree* t2;
    TreeEdge* e;
    TreeEdgeIterator T;
    int dir;
    REAL eps;
    
    double start_time = get_time();
    
    if (IS_INT) {
        if (options.dual_greedy_update_option == 2) {
            printf("Fixed eps approach can only be used with floating point REAL!\n");
            printf("Change REAL to double in PerfectMatching.h and recompile\n");
            exit(1);
        }
        if (options.dual_LP_threshold > 0) {
            printf("LP approach can only be used with floating point REAL!\n");
            printf("Change REAL to double in PerfectMatching.h and recompile\n");
            exit(1);
        }
    }
    if (options.verbose) {
        printf("perfect matching with %d nodes and %d edges\n", node_num, edge_num);
        fflush(stdout);
    }
    
    if (first_solve) {
        if (options.verbose) {
            printf("    starting init...");
            fflush(stdout);
        }
        if (options.fractional_jumpstart)
            InitGlobal();
        else
            InitGreedy();
        if (options.verbose)
            printf("done [%.3f secs]. ", get_time() - start_time);
        first_solve = false;
    } else if (options.verbose)
        printf("    solving updated problem. ");
    
    if (options.verbose) {
        printf("%d trees\n    .", tree_num);
        fflush(stdout);
    }
    
    memset(&stat, 0, sizeof(Stat));
    
    ///////////////////////////////////////////////////////
    //       first pass - initialize auxiliary graph     //
    ///////////////////////////////////////////////////////
    
    for (r = nodes[node_num].tree_sibling_next; r; r = r->tree_sibling_next) {
        //assert(!r->is_processed);
        t = r->tree;
        //assert(!t->first[0] && !t->first[1]);
        
        EdgeIterator I;
        FOR_ALL_EDGES(r, a, dir, I)
        {
            j = a->head[dir];
            if (j->flag == 2)
                t->pq0.Add(a);
            else if (j->is_processed) {
                //assert(j->flag == 0);
                if (!j->tree->pq_current)
                    AddTreeEdge(t, j->tree);
                j->tree->pq_current->pq00.Add(a);
            }
        }
        r->is_processed = 1;
        FOR_ALL_TREE_EDGES(t, e, dir)
            e->head[dir]->pq_current = NULL;
    }
    
    ///////////////////////////////////////////////////////
    //                  main loop                        //
    ///////////////////////////////////////////////////////
    
    while (1) {
        int tree_num0 = tree_num;
        Stat stat0 = stat;
        REAL delta = 0;
        
        for (r = nodes[node_num].tree_sibling_next; r;) {
            r2 = r->tree_sibling_next;
            if (r2)
                r3 = r2->tree_sibling_next;
            t = r->tree;
            
            int tree_num1 = tree_num;
            
            //////////////////////////////////////////////////////////////////////
            // step 1 - traversing auxiliary graph, setting pq_current pointers //
            //////////////////////////////////////////////////////////////////////
            t->pq_current = t;
            if (options.update_duals_before) {
                eps = PM_INFTY;
                Edge* a_augment = NULL;
                REAL eps_augment = PM_INFTY;
                if ((q = t->pq0.GetMin()))
                    eps = q->slack;
                if ((q = t->pq_blossoms.GetMin()) && eps > q->slack)
                    eps = q->slack;
                while ((q = t->pq00.GetMin())) {
                    if (ProcessEdge00((Edge*) q, false))
                        break;
                    t->pq00.Remove(q, pq_buf);
                }
                if (q && 2 * eps > q->slack)
                    eps = q->slack / 2;
                FOR_ALL_TREE_EDGES_X(t, e, dir, T)
                {
                    t2 = e->head[dir];
                    t2->pq_current = e;
                    t2->dir_current = dir;
                    if ((q = e->pq00.GetMin()) && (!a_augment || eps_augment > q->slack - t2->eps)) {
                        a_augment = (Edge*) q;
                        eps_augment = q->slack - t2->eps;
                    }
                    if ((q = e->pq01[dir].GetMin()) && eps > q->slack + t2->eps)
                        eps = q->slack + t2->eps;
                }
                if (eps > eps_augment)
                    eps = eps_augment;
                if (eps > t->eps) {
                    delta += eps - t->eps;
                    t->eps = eps;
                }
                if (a_augment && eps_augment <= t->eps)
                    Augment(a_augment);
            } else {
                FOR_ALL_TREE_EDGES_X(t, e, dir, T)
                {
                    t2 = e->head[dir];
                    t2->pq_current = e;
                    t2->dir_current = dir;
                    
                    if ((q = e->pq00.GetMin()) && (q->slack - t->eps <= t2->eps)) {
                        Augment((Edge*) q);
                        break;
                    }
                }
            }
            
            /////////////////////////////////
            //   step 2 - growing tree     //
            /////////////////////////////////
            eps = t->eps;
            REAL twice_eps = 2 * eps;
            
            while (tree_num1 == tree_num) {
                if ((q = t->pq0.GetMin()) && q->slack <= t->eps) {
                    a = (Edge*) q;
                    dir = (a->head[1]->flag == 2 && a->head[1]->is_outer) ? 1 : 0;
                    GET_OUTER_HEAD(a, 1-dir, i);
                    j = a->head[dir];
                    //assert(i->flag==0 && j->flag==2 && i->is_outer && j->is_outer && i->tree==t);
                    
                    j->flag = 1;
                    j->tree = i->tree;
                    j->tree_parent = EDGE_DIR_TO_ARC(a, 1-dir);
                    j->y += eps;
                    j = ARC_HEAD(j->match) ;
                    j->y -= eps;
                    ADD_TREE_CHILD(i, j);
                    
                    GrowTree(j, true);
                } else if ((q = t->pq00.GetMin()) && q->slack <= twice_eps) {
                    t->pq00.Remove(q, pq_buf);
                    a = (Edge*) q;
                    if (ProcessEdge00(a))
                        Shrink(a);
                } else if ((q = t->pq_blossoms.GetMin()) && q->slack <= eps) {
                    t->pq_blossoms.Remove(q, pq_buf);
                    a = (Edge*) q;
                    j = (a->head[0]->flag == 1) ? a->head[0] : a->head[1];
                    REAL tmp = a->slack;
                    a->slack = j->y;
                    j->y = tmp;
                    Expand(j);
                } else
                    break;
            }
            
            ///////////////////////////////////////////////////////////////////////
            // step 3 - traversing auxiliary graph, clearing pq_current pointers //
            ///////////////////////////////////////////////////////////////////////
            if (tree_num1 == tree_num) {
                t->pq_current = NULL;
                if (options.update_duals_after) {
                    eps = PM_INFTY;
                    Edge* a_augment = NULL;
                    REAL eps_augment = PM_INFTY;
                    if ((q = t->pq0.GetMin()))
                        eps = q->slack;
                    if ((q = t->pq_blossoms.GetMin()) && eps > q->slack)
                        eps = q->slack;
                    while ((q = t->pq00.GetMin())) {
                        if (ProcessEdge00((Edge*) q, false))
                            break;
                        t->pq00.Remove(q, pq_buf);
                    }
                    if (q && 2 * eps > q->slack)
                        eps = q->slack / 2;
                    FOR_ALL_TREE_EDGES(t, e, dir)
                    {
                        t2 = e->head[dir];
                        e->head[dir]->pq_current = NULL;
                        if ((q = e->pq00.GetMin()) && (!a_augment
                                || eps_augment > q->slack - t2->eps)) {
                            a_augment = (Edge*) q;
                            eps_augment = q->slack - t2->eps;
                        }
                        if ((q = e->pq01[dir].GetMin()) && eps > q->slack + t2->eps)
                            eps = q->slack + t2->eps;
                    }
                    if (eps > eps_augment)
                        eps = eps_augment;
                    bool progress = false;
                    if (eps > t->eps) {
                        delta += eps - t->eps;
                        t->eps = eps;
                        progress = true;
                    }
                    if (a_augment && eps_augment <= t->eps)
                        Augment(a_augment);
                    else if (progress && tree_num >= options.single_tree_threshold * node_num) {
                        // continue with the same tree
                        r = t->root;
                        continue;
                    }
                } else {
                    FOR_ALL_TREE_EDGES(t, e, dir)
                        e->head[dir]->pq_current = NULL;
                }
            }
            
            ///////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////
            
            r = r2;
            if (r && !r->is_tree_root)
                r = r3;
        }
        
        if (tree_num == 0)
            break;
        
        if (tree_num == tree_num0)
        //&& stat.grow_count == stat0.grow_count
        //&& stat.shrink_count == stat0.shrink_count
        //&& stat.expand_count == stat0.expand_count )
        {
            if (!UpdateDuals()) {
                if (!IS_INT && delta <= PM_THRESHOLD )    // for numerical stability
                {
                    //CommitEps();
                    int dual_greedy_update_option = options.dual_greedy_update_option;
                    options.dual_greedy_update_option = 2;
                    UpdateDuals();
                    options.dual_greedy_update_option = dual_greedy_update_option;
                }
            }
        }
    }
    
    if (finish)
        Finish();
    
    if (options.verbose) {
        printf("\ndone [%.3f secs]. %d grows, %d expands, %d shrinks\n", get_time() - start_time,
                stat.grow_count, stat.expand_count, stat.shrink_count);
        printf("    expands: [%.3f secs], shrinks: [%.3f secs], dual updates: [%.3f secs]\n",
                stat.expand_time, stat.shrink_time, stat.dual_time);
        fflush(stdout);
    }
}

//----------//
//  Repair  //
//----------//

struct PerfectMatching::LCATreeX: LCATree {
        LCATreeX(int size)
                : LCATree(size) {
            rev_mapping = new Node*[size];
        }
        ~LCATreeX() {
            delete[] rev_mapping;
        }
        Node** rev_mapping;
};

void PerfectMatching::StartUpdate() {
    Node* i0;
    Node* i;
    Node* j;
    Node* b;
    
    while ((i = removed_first)) {
        removed_first = i->tree_sibling_next;
        blossoms->Delete(i);
        removed_num--;
    }
    
    Edge* a;
    Edge* selfloop_first = NULL;
    Edge* selfloop_last = NULL;
    
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        i0->is_processed = 0;
        if (i0->is_outer)
            continue;
        
        i0->is_tree_root = 0;
        i0->blossom_ptr = NULL;
        i = i0;
        while (1) {
            j = i->blossom_parent;
            j->is_processed = 0;
            if (j->is_outer) {
                j->first_tree_child = i;
                break;
            }
            if (j->is_marked)
                break;
            if ((a = j->blossom_selfloops)) {
                if (selfloop_last)
                    selfloop_last->next[1] = a;
                else
                    selfloop_first = a;
                selfloop_last = a;
                a->next[1] = NULL;
            }
            j->blossom_ptr = i;
            i = j;
        }
        b = (i->blossom_parent->is_outer) ? i->blossom_parent :
                                            i->blossom_parent->blossom_grandparent;
#ifdef LCA_REPAIRS
        if (!b->is_marked) {
            b->lca_size = 1;
            b->is_marked = 1;
        }
#endif
        while (1) {
#ifdef LCA_REPAIRS
            b->lca_size++;
#endif
            ARC_TO_EDGE_PTR(i->blossom_sibling) ->y_saved = i->y;
            i->y += i->blossom_parent->y;
            if (!i->is_blossom)
                break;
            i->is_marked = 1;
            j = i;
            i = i->blossom_ptr;
            j->blossom_grandparent = b;
        }
        i->blossom_grandparent = b;
    }
    
#ifdef LCA_REPAIRS
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        if (i0->is_outer)
            continue;
        b = i0->blossom_grandparent;
        if (!b->is_marked)
            continue;
        b->is_marked = 0;
        LCATreeX* lca = new LCATreeX(b->lca_size);
        b->blossom_ptr = b->first_tree_child;
        i = b;
        while (1) {
            if (i->blossom_ptr)
                i = i->blossom_ptr;
            else {
                while (1) {
                    if (i->is_outer)
                        break;
                    i->lca_preorder = lca->Add(i, i->blossom_parent);
                    lca->rev_mapping[i->lca_preorder] = i;
                    i = ARC_HEAD(i->blossom_sibling) ;
                    if (i != i->blossom_parent->blossom_ptr)
                        break;
                    i = i->blossom_parent;
                }
                if (i->is_outer) {
                    lca->AddRoot(i);
                    break;
                }
            }
        }
        b->lca = lca;
    }
#endif
    
    while ((a = selfloop_first)) {
        selfloop_first = a->next[1];
        do {
            Edge* a_next = a->next[0];
            
#ifdef LCA_REPAIRS
            int _i = a->head0[1]->lca_preorder;
            int _j = a->head0[0]->lca_preorder;
            Node* b = a->head0[1]->blossom_grandparent;
            b->lca->GetPenultimateNodes(_i, _j);
            i = b->lca->rev_mapping[_i];
            j = b->lca->rev_mapping[_j];
#else
            GetRealEndpoints(a, i, j);
#endif
            ADD_EDGE(i, a, 0);
            ADD_EDGE(j, a, 1);
            a->slack -= 2 * i->blossom_eps;
            a = a_next;
        } while (a);
    }
    
    /*
     for (i0=nodes; i0<nodes+node_num; i0++)
     {
     if (i0->is_outer) continue;
     b = i0->blossom_grandparent;
     if (b->lca)
     {
     delete b->lca;
     b->lca = NULL;
     }
     }
     */

    nodes[node_num].first_tree_child = NULL;
}

void PerfectMatching::FinishUpdate() {
    Node* i0;
    Node* i;
    Node* j;
    Edge* a;
    EdgeIterator I;
    int dir;
    Tree* t;
    
    for (i0 = nodes; i0 < nodes + node_num; i0++) {
        if (i0->is_outer)
            continue;
        
#ifdef LCA_REPAIRS
        if (i0->blossom_grandparent->lca) {
            delete i0->blossom_grandparent->lca;
            i0->blossom_grandparent->lca = NULL;
        }
#endif
        
        //////////////////////////////////////////////////////////////
        if (!i0->blossom_grandparent->is_removed) {
            i = i0;
            do {
                i->y = ARC_TO_EDGE_PTR(i->blossom_sibling) ->y_saved;
                i->is_marked = 0;
                i->blossom_selfloops = NULL;
                i = i->blossom_parent;
            } while (i->is_marked);
            continue;
        }
        //////////////////////////////////////////////////////////////
        
        i = i0->blossom_parent;
        while (1) {
            if (i->is_removed && !i->is_outer)
                break;
            REAL y_parent = (i->is_outer) ? 0 : i->blossom_parent->y;
            for (dir = 0; dir < 2; dir++) {
                if (!i->first[dir])
                    continue;
                i->first[dir]->prev[dir]->next[dir] = NULL;
                Edge* a_next;
                for (a = i->first[dir]; a; a = a_next) {
                    a_next = a->next[dir];
                    j = a->head0[1 - dir];
                    ADD_EDGE(j, a, dir);
                    a->slack += j->blossom_parent->y - y_parent;
                }
                i->first[dir] = NULL;
            }
            if (i->is_removed)
                break;
            
            j = i->blossom_parent;
            i->is_removed = 1;
            i->tree_sibling_next = removed_first;
            removed_first = i;
            i = j;
        }
        i0->y = ARC_TO_EDGE_PTR(i0->blossom_sibling) ->y_saved;
        i0->is_outer = 1;
        i0->flag = 2;
        i0->is_tree_root = 1;
    }
    
    Node* blossom_list = nodes[node_num].first_tree_child;
    
    for (i = nodes; i < nodes + node_num; i++) {
        if (!i->is_tree_root)
            continue;
        i->first_tree_child = nodes[node_num].first_tree_child;
        nodes[node_num].first_tree_child = i;
        REAL slack_min = PM_INFTY;
        FOR_ALL_EDGES(i, a, dir, I)
        {
            if (slack_min > a->slack)
                slack_min = a->slack;
        }
        i->y += slack_min;
        FOR_ALL_EDGES(i, a, dir, I)
            a->slack -= slack_min;
    }
    
    tree_num = 0;
    for (i = nodes[node_num].first_tree_child; i != blossom_list; i = i->first_tree_child) {
        tree_num++;
        if (!i->is_tree_root)
            continue;
        FOR_ALL_EDGES(i, a, dir, I)
        {
            j = a->head[dir];
            if (a->slack <= 0 && j->is_tree_root) {
                i->is_tree_root = j->is_tree_root = 0;
                i->match = EDGE_DIR_TO_ARC(a, dir);
                j->match = EDGE_DIR_TO_ARC(a, 1-dir);
                tree_num -= 2;
                break;
            }
        }
    }
    for (; i; i = i->first_tree_child) {
        if (i->is_removed) {
            i->is_tree_root = 0;
            continue;
        }
        tree_num++;
    }
    
    if (tree_num > tree_num_max) {
        if (trees)
            free(trees);
        tree_num_max = tree_num;
        trees = (Tree*) malloc(tree_num_max * sizeof(Tree));
    }
    t = trees;
    
    Node* last_root = &nodes[node_num];
    Node* i_next;
    for (i = nodes; i; i = i_next) {
        if (!i->is_blossom)
            i_next = (i < nodes + node_num) ? (i + 1) : blossom_list;
        else
            i_next = i->first_tree_child;
        if (!i->is_tree_root)
            continue;
        
        i->flag = 0;
        i->first_tree_child = NULL;
        i->tree_sibling_prev = last_root;
        last_root->tree_sibling_next = i;
        last_root = i;
        i->tree = t;
        t->root = i;
        t->eps = 0;
        t->first[0] = t->first[1] = NULL;
        t->pq_current = NULL;
        t->pq00.Reset();
        t->pq0.Reset();
        t->pq_blossoms.Reset();
        t++;
    }
    
    assert(t == trees + tree_num);
    last_root->tree_sibling_next = NULL;
    
    while ((i = removed_first)) {
        removed_first = i->tree_sibling_next;
        blossoms->Delete(i);
        blossom_num--;
    }
}

PerfectMatching::REAL PerfectMatching::GetTwiceSum(NodeId i) {
    assert(i>=0 && i<node_num);
    return nodes[i].y;
}

inline void PerfectMatching::ProcessNegativeEdge(Edge* a) {
    int dir;
    Node* i;
    for (dir = 0; dir < 2; dir++) {
        i = a->head0[dir];
        if (i->is_outer) {
            if (!i->is_tree_root) {
                i->is_tree_root = 1;
                i = ARC_HEAD(i->match) ;
                assert(!i->is_tree_root && i->is_outer);
                i->is_tree_root = 1;
                if (i->is_blossom) {
                    i->first_tree_child = nodes[node_num].first_tree_child;
                    nodes[node_num].first_tree_child = i;
                }
            }
            return;
        }
        if (i->blossom_grandparent->is_removed)
            return;
    }
    
    Node* b = i->blossom_grandparent;
    assert(b->is_outer);
    
    if (!b->is_tree_root) {
        b->is_tree_root = 1;
        i = ARC_HEAD(b->match) ;
        assert(!i->is_tree_root && i->is_outer);
        i->is_tree_root = 1;
        if (i->is_blossom) {
            i->first_tree_child = nodes[node_num].first_tree_child;
            nodes[node_num].first_tree_child = i;
        }
    }
    
    b->is_removed = 1;
    b->tree_sibling_next = removed_first;
    removed_first = b;
}

PerfectMatching::EdgeId PerfectMatching::AddNewEdge(NodeId _i, NodeId _j, REAL cost,
                                                    bool do_not_add_if_positive_slack) {
    assert(_i>=0 && _i<node_num && _j>=0 && _j<node_num && _i!=_j);
    if (edge_num >= edge_num_max)
        ReallocateEdges();
    Node* i = nodes + _i;
    Node* j = nodes + _j;
    Edge* a = edges + edge_num;
    
    a->slack = cost * COST_FACTOR;
    a->head0[0] = j;
    a->head0[1] = i;
    Node* bi = (i->is_outer) ? i : i->blossom_grandparent;
    Node* bj = (j->is_outer) ? j : j->blossom_grandparent;
    if (bi == bj) {
#ifdef LCA_REPAIRS
        int _i = i->lca_preorder;
        int _j = j->lca_preorder;
        bi->lca->GetPenultimateNodes(_i, _j);
        i = bi->lca->rev_mapping[_i];
        j = bi->lca->rev_mapping[_j];
#else
        GetRealEndpoints(a, i, j);
#endif
        a->slack += i->blossom_parent->y + j->blossom_parent->y;
    } else {
        i = bi;
        j = bj;
    }
    a->slack -= a->head0[0]->y + a->head0[1]->y;
    
    if (a->slack >= 0 && do_not_add_if_positive_slack)
        return -1;
    
    ADD_EDGE(i, a, 0);
    ADD_EDGE(j, a, 1);
    PriorityQueue<REAL>::ResetItem(a);
    
    if (a->slack < 0) {
        ProcessNegativeEdge(a);
    }
    
    return edge_num++;
}

void PerfectMatching::UpdateCost(EdgeId e, REAL delta_cost) {
    assert(e>=0 && e<edge_num);
    Edge* a = edges + e;
    a->slack += delta_cost * COST_FACTOR;
    if (a->slack == 0)
        return;
    if (a->slack > 0) {
        Node* i = a->head[1];
        Node* j = a->head[0];
        if (i->is_outer) {
            if (ARC_TO_EDGE_PTR(i->match) != a && ARC_TO_EDGE_PTR(j->match) != a)
                return;
        } else {
            if (ARC_TO_EDGE_PTR(i->blossom_sibling) != a && ARC_TO_EDGE_PTR(j->blossom_sibling)
                    != a)
                return;
        }
    }
    ProcessNegativeEdge(a);
}

//----------//
//  Shrink  //
//----------//

PerfectMatching::Node* PerfectMatching::FindBlossomRoot(Edge* a0) {
    Node* i;
    Node* j;
    Node* _i[2];
    Node* r;
    int branch;
    
    _i[0] = ARC_HEAD(a0) ;
    _i[1] = ARC_TAIL(a0) ;
    branch = 0;
    while (1) {
        if (_i[branch]->is_marked) {
            r = _i[branch];
            j = _i[1 - branch];
            break;
        }
        _i[branch]->is_marked = 1;
        if (_i[branch]->is_tree_root) {
            j = _i[branch];
            i = _i[1 - branch];
            while (!i->is_marked) {
                i->is_marked = 1;
                i = ARC_HEAD(i->match) ;
                GET_TREE_PARENT(i, i);
            }
            r = i;
            break;
        }
        i = ARC_HEAD(_i[branch]->match) ;
        GET_TREE_PARENT(i, _i[branch]);
        branch = 1 - branch;
    }
    i = r;
    while (i != j) {
        i = ARC_HEAD(i->match) ;
        i = ARC_HEAD(i->tree_parent);
        i->is_marked = 0;
    }   
        // clear is_marked and is_outer
    i = ARC_HEAD(a0) ;
    while (i != r) {
        i->is_marked = 0;
        i->is_outer = 0;
        i = ARC_HEAD(i->match) ;
        i->is_outer = 0;
        i = ARC_HEAD(i->tree_parent) ;
    }   
    i = ARC_TAIL(a0) ;
    while (i != r) {
        i->is_marked = 0;
        i->is_outer = 0;
        i = ARC_HEAD(i->match) ;
        i->is_outer = 0;
        i = ARC_HEAD(i->tree_parent) ;
    }   
    r->is_marked = 0;
    r->is_outer = 0;
    
    return r;
}

void PerfectMatching::Shrink(Edge* a0) {
    //assert(a0->head[0]->is_outer && a0->head[1]->is_outer);
    //assert(a0->head[0]->flag == 0 && a0->head[1]->flag == 0);
    
    double start_time = get_time();
    
    int branch, dir;
    Node* r;
    Node* i;
    Node* j;
    Edge* a;
    Edge** a_inner_ptr;
    Arc* a_prev;
    Node* b = blossoms->New();
    Edge* a_augment = NULL;
    Edge* b_match;
    
    b->first[0] = b->first[1] = NULL;
    
    // set is_outer=0 for all nodes in the blossom
    r = FindBlossomRoot(a0);
    Tree* t = r->tree;
    REAL eps = t->eps;
    
    b->first_tree_child = NULL;
    i = ARC_HEAD(a0) ;
    branch = 0;
    while (1) {
        if (i == r && branch)
            break;
        i->is_marked = 1;
        if (i == r) {
            branch = 1;
            i = ARC_TAIL(a0) ;
            continue;
        }
        
        // remove i from the list of children
        REMOVE_FROM_TREE(i);
        
        // move children of i to the list of children of b
        if (i->first_tree_child) {
            j = i->first_tree_child;
            if (!b->first_tree_child)
                b->first_tree_child = j;
            else {
                Node* j_last = j->tree_sibling_prev;
                j->tree_sibling_prev = b->first_tree_child->tree_sibling_prev;
                b->first_tree_child->tree_sibling_prev->tree_sibling_next = j;
                b->first_tree_child->tree_sibling_prev = j_last;
            }
        }
        
        // go to parent
        i = ARC_HEAD(i->match) ;
        i->is_marked = 1;
        if (i->is_blossom) {
            a = ARC_TO_EDGE_PTR(i->match);
            t->pq_blossoms.Remove(a, pq_buf);
            REAL tmp = a->slack;
            a->slack = i->y;
            i->y = tmp;
        }
        i = ARC_HEAD(i->tree_parent) ;
    }   

        // move children of r to the list of children of b
    if (i->first_tree_child) {
        j = i->first_tree_child;
        if (!b->first_tree_child)
            b->first_tree_child = j;
        else {
            Node* j_last = j->tree_sibling_prev;
            j->tree_sibling_prev = b->first_tree_child->tree_sibling_prev;
            b->first_tree_child->tree_sibling_prev->tree_sibling_next = j;
            b->first_tree_child->tree_sibling_prev = j_last;
        }
    }
    
    // init b
    b->is_removed = 0;
    b->is_outer = 1;
    b->flag = 0;
    b->is_blossom = 1;
    b->is_tree_root = r->is_tree_root;
    b->is_processed = 1;
    b->tree = t;
    b->y = -eps;
    b->is_marked = 0;
    
    // replace r with b in the tree
    b->tree_sibling_prev = r->tree_sibling_prev;
    b->tree_sibling_next = r->tree_sibling_next;
    Node* b_parent = NULL;
    if (!b->is_tree_root) {
        b_parent = ARC_HEAD(r->match) ; GET_TREE_PARENT(b_parent, b_parent);
    }   
    if (b->tree_sibling_prev->tree_sibling_next)
        b->tree_sibling_prev->tree_sibling_next = b;
    else
        b_parent->first_tree_child = b;
    if (b->tree_sibling_next)
        b->tree_sibling_next->tree_sibling_prev = b;
    else if (b_parent)
        b_parent->first_tree_child->tree_sibling_prev = b;
    
    if (b->is_tree_root) {
        b->tree->root = b;
        b_match = NULL;
    } else {
        b->match = r->match;
        b_match = ARC_TO_EDGE_PTR(b->match);
    }
    REAL b_match_slack = 0;    // initialize to prevent compiler warning
    if (b_match && ARC_HEAD(b->match) ->is_blossom)
    {   
        b_match_slack = b_match->slack;
        b_match->slack = ARC_HEAD(b->match)->y;
    }

    // second pass over nodes in the blossom
    branch = 0;
    a_prev = EDGE_DIR_TO_ARC(a0, 0);
    i = ARC_HEAD(a_prev) ;
    while (1) {
        // update Arc::next and Arc::head pointers
        if (i->flag == 0)
            i->y += eps;
        else
            i->y -= eps;
        i->is_processed = 0;
        
        if (i->flag == 1) {
            Edge* a_prev;
            for (dir = 0; dir < 2; dir++)
                if (i->first[dir]) {
                    for (a_inner_ptr = &i->first[dir], a = *a_inner_ptr, a_prev = a->prev[dir], a_prev->next[dir] =
                            NULL; a; a = *a_inner_ptr) {
                        Node* j0 = a->head[dir];
                        for (j = j0; !j->is_outer && !j->is_marked; j = j->blossom_parent) {
                        }
                        if (j != j0) { /*assert(j->flag == 0);*/
                            int dir_rev = 1 - dir;
                            MOVE_EDGE(j0, j, a, dir_rev);
                        }
                        if (j->is_marked)    // "inner" arc
                        {
                            a_inner_ptr = &a->next[dir];
                            a->prev[dir] = a_prev;
                            a_prev = a;
                            
                            if (j->flag == 1)
                                a->slack += eps;
                        } else    // "boundary" arc
                        {
                            *a_inner_ptr = a->next[dir];
                            ADD_EDGE(b, a, dir);
                            
                            if (j->flag == 0 && j->tree != t) {
                                j->tree->pq_current->pq01[1 - j->tree->dir_current].Remove(a, pq_buf);
                                if (a->slack + eps <= j->tree->eps)
                                    a_augment = a;
                            }
                            a->slack += 2 * eps;
                            if (j->flag == 2)
                                t->pq0.Add(a);
                            else if (j->flag == 0) {
                                if (!j->tree->pq_current)
                                    AddTreeEdge(t, j->tree);
                                j->tree->pq_current->pq00.Add(a);
                            } else if (j->tree != t) {
                                if (!j->tree->pq_current)
                                    AddTreeEdge(t, j->tree);
                                j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
                            }
                        }
                    }
                    if (i->first[dir]) {
                        a_prev->next[dir] = i->first[dir];
                        i->first[dir]->prev[dir] = a_prev;
                    }
                }
        }
        
        Arc* a_next = (i->flag == 0) ? i->match : i->tree_parent;
        i->blossom_parent = b;
        i->match = NULL;
        i->blossom_grandparent = b;
        i->blossom_selfloops = NULL;
        if (branch == 0) {
            i->blossom_sibling = a_next;
            if (i == r) {
                branch = 1;
                a_prev = ARC_REV(a0);
                i = ARC_HEAD(a_prev) ;
                if (i == r)
                    break;
            } else {
                a_prev = i->blossom_sibling;
                i = ARC_HEAD(a_prev) ;
            }
        }
        else
        {   
            i->blossom_sibling = ARC_REV(a_prev);
            a_prev = a_next;
            i = ARC_HEAD(a_prev);
            if (i == r) break;
        }
    }           
    i->blossom_sibling = ARC_REV(a_prev);
    r->is_tree_root = 0;
    
    for (i = ARC_HEAD(r->blossom_sibling) ;; i = ARC_HEAD(i->blossom_sibling))
    {   
        i->is_marked = 0;
        i->blossom_eps = eps;
        if (i == r) break;
    }

    if (b_match) {
        if (ARC_HEAD(b->match) ->is_blossom)
        {   
            b_match->slack = b_match_slack;
        }
        dir = ARC_TO_EDGE_DIR(b->match);
        //assert(b_match->head[1-dir] == r);
        MOVE_EDGE(r, b, b_match, dir);
    }

    stat.shrink_count++;
    blossom_num++;
    stat.shrink_time += get_time() - start_time;
    
    if (a_augment)
        Augment(a_augment);
}

