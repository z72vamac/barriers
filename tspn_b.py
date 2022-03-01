# TSPN-B

import time
import itertools

import gurobipy as gp
from gurobipy import GRB
from matplotlib.patches import Circle

import auxiliar_functions as af
from data import *
import neighborhood as neigh


def tspn_b(barriers, neighborhoods, prepro=True, log=False, dominant = False, picture=False, time_limit=7200, init=False):

    if not(dominant):
        first_time = time.time()

        vertices_neighborhood = list(itertools.product(range(-len(neighborhoods), 0), range(1)))
        vertices_neighborhood = vertices_neighborhood[::-1]
        # list(itertools.product(vertices_neighborhood, range(len(barriers)), range(2))) + list(itertools.product(range(len(barriers)),
        # range(2), vertices_neighborhood))

        edges_neighborhood = []

        for (a, b) in vertices_neighborhood:
            for c in range(len(barriers)):
                for d in range(2):
                    if prepro:
                        # Point of the barrier to check if is visible by the neighborhood
                        point = barriers[c][d]

                        # Neighborhood to check if it is visible by the point
                        neighborhood = neighborhoods[abs(a) - 1]

                        if af.cansee(point, neighborhood, barriers):
                            # Appending the feasible edges to edges_neighborhood
                            edges_neighborhood.append((a, b, c, d))
                            edges_neighborhood.append((c, d, a, b))
                    else:
                        edges_neighborhood.append((a, b, c, d))
                        edges_neighborhood.append((c, d, a, b))

        vertices_barrier = list(itertools.product(range(len(barriers)), range(2)))

        edges_barrier = []
        for v, i in vertices_barrier:
            for w, j in vertices_barrier:
                if v != w:
                    if prepro:
                        barrier = [barriers[v][i], barriers[w][j]]

                        if any([not (af.intersect(barrieri, barrier)) for barrieri in barriers]):
                            edges_barrier.append((v, i, w, j))
                    else:
                        edges_barrier.append((v, i, w, j))

        vertices_total = vertices_neighborhood + vertices_barrier
        edges_total = edges_neighborhood + edges_barrier

        if log:
            print("vertices_neighborhood = " + str(vertices_neighborhood))
            print("vertices_barrier = " + str(vertices_barrier))

            print("edges_neighborhood = " + str(edges_neighborhood))
            print("edges_barrier = " + str(edges_barrier))

        epsilon_index = []  # epsilon(S / T, B, i) = 1 si (P_{S/T}, P_B^i)\in E_{S/T}

        for a, b, c, d in edges_neighborhood:
            if a < 0:
                epsilon_index.append((a, b, c, d))

        # print(epsilon_index)

        point_index = []
        for a, b in vertices_neighborhood:
            for dim in range(2):
                point_index.append((a, b, dim))

        if log:
            print("epsilon = " + str(epsilon_index))

        y_index = edges_total

        if log:
            print("y = " + str(y_index))

        dist_index = edges_total

        if log:
            print("dist = " + str(dist_index))

        dif_index = []

        for a, b, c, d in edges_total:
            for dim in range(2):
                dif_index.append((a, b, c, d, dim))

        if log:
            print("dif = " + str(dif_index))

        delta_index = []  # delta(S / T, B, i, B') = 1 si (P_{S/T}, P_B^i) y B' do not intersect.

        for a, b, c, d in epsilon_index:
            for e in range(len(barriers)):
                delta_index.append((a, b, c, d, e))

        if log:
            print("delta = " + str(delta_index))

        gamma_index = []  # gamma(S / T, B, i, B') = alpha(S / T, B, i, B', 0, order2)*alpha(S / T, B, i, B', 1, order2)

        for a, b, c, d, e in delta_index:
            for f in range(2):
                gamma_index.append((a, b, c, d, e, f))

        if log:
            print("gamma = " + str(gamma_index))

        # beta(S/T, B, i, B', order) = 1 if sign(P_{S/T} P_B^i | B') is the same (order = 0); sign(B' | P_{S/T} P_B^i) is
        # the same (order = 1);
        beta_index = gamma_index

        beta0_index = []
        beta1_index = []

        for a, b, c, d, e, h in beta_index:
            if h == 0:
                beta0_index.append((a, c, d, e))
            else:
                beta1_index.append((a, c, d, e))

        beta0_index = list(set(beta0_index))
        beta1_index = list(set(beta1_index))

        if log:
            print("beta = " + str(beta_index))

        # alpha(S/T, B, i, B', order1, order2) = 1 si sign(P_{S/T} | B') > 0 (order1 = 0); sign(P_B^i | B') > 0 (order1)
        # = 1;
        alpha_index = []  # list(itertools.product(vertices_neighborhood, range(len(barriers)), range(2), range(len(barriers)),
        # range(2), range(2)))

        for a, b, c, d, e, f in gamma_index:
            for g in range(2):
                alpha_index.append((a, b, c, d, e, f, g))

        if log:
            print("alpha = " + str(alpha_index))

        alpha00_index = []
        alpha10_index = []
        alpha_1_index = []

        for a, b, c, d, e, f, h in alpha_index:
            if f == 0 and h == 0:
                alpha00_index.append((a, e))

            if f == 1 and h == 0:
                alpha10_index.append((c, d, e))

            if h == 1:
                alpha_1_index.append((a, c, d, e, f))

        alpha00_index = list(set(alpha00_index))
        alpha10_index = list(set(alpha10_index))
        alpha_1_index = list(set(alpha_1_index))

        # P_S and P_T: indices of the points in the neighborhoods
        p_index = edges_neighborhood
        # for index in vertices_neighborhood:
        #     for dim in range(2):
        #         p_index.append((index[0], index[1], dim))

        if log:
            print("p_index = " + str(p_index))

        # socp variables:
        d_inside_index = vertices_neighborhood

        dif_inside_index = []

        for (a, b) in vertices_neighborhood:
            for dim in range(2):
                dif_inside_index.append((a, b, dim))

        if log:
            print("d_inside_index = " + str(d_inside_index))
            print("dif_inside_index = " + str(dif_inside_index))

        # z variables:
        # z_index = list(itertools.product(vertices_neighborhood, vertices_neighborhood))
        # print("z_index = " + str(z_index))

        # f variables:
        g_index = y_index

        if log:
            print("g_index = " + str(g_index))

        model = gp.Model('HTSPS_Model')

        epsilon = model.addVars(epsilon_index, vtype=GRB.BINARY, name='epsilon')
        p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
        y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
        dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')
        dif = model.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
        delta = model.addVars(delta_index, vtype=GRB.BINARY, name='delta')
        gamma = model.addVars(gamma_index, vtype=GRB.BINARY, name='gamma')
        # beta = model.addVars(beta_index, vtype = GRB.BINARY, name = 'beta')
        # alpha = model.addVars(alpha_index, vtype = GRB.BINARY, name = 'alpha')

        alpha00 = model.addVars(alpha00_index, vtype=GRB.BINARY, name='alpha00')
        alpha10 = model.addVars(alpha10_index, vtype=GRB.BINARY, name='alpha10')
        alpha_1 = model.addVars(alpha_1_index, vtype=GRB.BINARY, name='alpha_1')

        beta0 = model.addVars(beta0_index, vtype=GRB.BINARY, name='beta0')
        beta1 = model.addVars(beta1_index, vtype=GRB.BINARY, name='beta1')

        # point = model.addVars(p_index, vtype = GRB.CONTINUOUS, lb = 0.1, ub = 99.9, name = 'point')
        point = model.addVars(point_index, vtype=GRB.CONTINUOUS, name='point')

        d_inside = model.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='d_inside')
        dif_inside = model.addVars(dif_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif_inside')
        landa = model.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa')
        # z = model.addVars(z_index, vtype = GRB.BINARY, name = 'z')
        g = model.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name='g')

        model.update()

        if init:
            time_h, objval_h = heuristic(barriers, neighborhoods)
        # Alpha constraint

        # alpha-C for a, b, c, d, e, f, h in alpha.keys(): # print((a, b, c, d, e, f, h)) if f == 0 and h == 0:
        # model.addConstr((1 - alpha[a, b, c, d, e, f, h])*small_m <= af.determinant([point[a, 0, 0], point[a, 0, 1]],
        # barriers[e][0], barriers[e][1])) model.addConstr(af.determinant([point[a, 0, 0], point[a, 0, 1]], barriers[e][
        # 0], barriers[e][1]) <= big_m*alpha[a, b, c, d, e, f, h]) # elif e == 1 and f == 0: #     model.addConstr((1 -
        # alpha[a, b, c, d, e, f])*small_m <= af.determinant(barriers[b][c], barriers[d][0], barriers[d][1])) #
        # model.addConstr(af.determinant(barriers[b][c], barriers[d][0], barriers[d][1]) <= big_m*alpha[a, b, c, d, e, f])
        #
        # elif f == 1 and h == 0: model.addConstr((1 - alpha[a, b, c, d, e, f, h])*abs(af.determinant(barriers[c][d],
        # barriers[e][0], barriers[e][1]))*(-1) <= af.determinant(barriers[c][d], barriers[e][0], barriers[e][1])) if
        # af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0: model.addConstr(alpha[a, b, c, d, e, f,
        # h] >= af.determinant(barriers[c][d], barriers[e][0], barriers[e][1])/abs(af.determinant(barriers[c][d],
        # barriers[e][0], barriers[e][1])))
        #
        # else: model.addConstr((1 - alpha[a, b, c, d, e, f, h])*small_m <= af.determinant(barriers[e][f], [point[a, 0,
        # 0], point[a, 0, 1]], barriers[c][d])) model.addConstr(af.determinant(barriers[e][f], [point[a, 0, 0], point[a,
        # 0, 1]], barriers[c][d]) <= big_m*alpha[a, b, c, d, e, f, h] small_m, big_m = estima_det(neighborhoods[0],
        # barriers[1])

        # model.addConstrs((1 - alpha00[a, e])*(estima_det(neighborhoods[abs(a)-1], barriers[e])[0]) <= af.determinant([
        # point[a, 0, 0], point[a, 0, 1]], barriers[e][0], barriers[e][1]) for a, e in alpha00.keys()) model.addConstrs(
        # af.determinant([point[a, 0, 0], point[a, 0, 1]], barriers[e][0], barriers[e][1]) <= estima_det(neighborhoods[
        # abs(a)-1], barriers[e])[1]*alpha00[a, e] for a, e in alpha00.keys())

        for a, e in alpha00.keys():
            small_m, big_m = af.estima_det(neighborhoods[abs(a) - 1], barriers[e])

            if big_m <= 0:
                alpha00[a, e] = 0

            elif small_m >= 0:
                alpha00[a, e] = 1

            else:

                model.addConstr(
                    (1 - alpha00[a, e]) * small_m <= af.determinant([point[a, 0, 0], point[a, 0, 1]], barriers[e][0],
                                                                    barriers[e][1]))
                model.addConstr(
                    af.determinant([point[a, 0, 0], point[a, 0, 1]], barriers[e][0], barriers[e][1]) <= big_m * alpha00[
                        a, e])

        for c, d, e in alpha10.keys():
            if af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0:
                alpha10[c, d, e] = (af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]) / abs(
                    af.determinant(barriers[c][d], barriers[e][0], barriers[e][1])) >= 0) * 1

            # model.addConstrs((1 - alpha10[c, d, e])*abs(af.determinant(barriers[c][d], barriers[e][0], barriers[e][
            # 1]))*(-1) <= af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys()
            # if af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0) model.addConstrs(abs(
            # af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]))*alpha10[c, d, e] >= af.determinant(
            # barriers[c][d], barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys() if af.determinant(barriers[
            # c][d], barriers[e][0], barriers[e][1]) != 0)

        # model.addConstrs((1 - alpha10[c, d, e])*abs(af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]))*(
        # -1) <= af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys() if
        # af.determinant(barriers[c][d], barriers[e][0], barriers[e][1]) != 0) model.addConstrs(abs(af.determinant(
        # barriers[c][d], barriers[e][0], barriers[e][1]))*alpha10[c, d, e] >= af.determinant(barriers[c][d],
        # barriers[e][0], barriers[e][1]) for c, d, e in alpha10.keys() if af.determinant(barriers[c][d], barriers[e][0],
        # barriers[e][1]) != 0)

        # model.addConstrs((1 - alpha_1[a, c, d, e, f])*estima_det(neighborhoods[abs(a)-1], [barriers[e][f], barriers[c][
        # d]])[0] <= af.determinant(barriers[e][f], [point[a, 0, 0], point[a, 0, 1]], barriers[c][d]) for a, c, d, e,
        # f in alpha_1.keys()) model.addConstrs(af.determinant(barriers[e][f], [point[a, 0, 0], point[a, 0, 1]],
        # barriers[c][d]) <= estima_det(neighborhoods[abs(a)-1], [barriers[e][f], barriers[c][d]])[1]*alpha_1[a, c, d, e,
        # f] for a, c, d, e, f in alpha_1.keys())

        small_m = -20000
        big_m = 20000

        model.addConstrs(
            (1 - alpha_1[a, c, d, e, f]) * small_m <= af.determinant(barriers[e][f], [point[a, 0, 0], point[a, 0, 1]],
                                                                     barriers[c][d]) for
            a, c, d, e, f in alpha_1.keys())
        model.addConstrs(
            af.determinant(barriers[e][f], [point[a, 0, 0], point[a, 0, 1]], barriers[c][d]) <= big_m * alpha_1[
                a, c, d, e, f] for
            a, c, d, e, f in alpha_1.keys())

        # model.addConstrs((1 - alpha_1[a, c, d, e, f])*estima_det(neighborhoods[abs(a)-1], [barriers[e][f], barriers[c][
        # d]])[0] <= af.determinant(barriers[e][f], [point[a, 0, 0], point[a, 0, 1]], barriers[c][d]) for a, c, d, e,
        # f in alpha_1.keys()) model.addConstrs(af.determinant(barriers[e][f], [point[a, 0, 0], point[a, 0, 1]],
        # barriers[c][d]) <= estima_det(neighborhoods[abs(a)-1], [barriers[e][f], barriers[c][d]])[1]*alpha_1[a, c, d, e,
        # f] for a, c, d, e, f in alpha_1.keys())

        # model.addConstr(alpha00[-2,0] == 0)
        # model.addConstr(alpha00[-2,1] == 1)

        # beta-C for a, b, c, d, e, h in beta.keys(): model.addConstr(beta[a, b, c, d, e, h] == 2*gamma[a, b, c, d, e,
        # h] - alpha[a, b, c, d, e, 0, h]-alpha[a, b, c, d, e, 1, h]+1) for a, b, c, d, e, h in beta.keys(): if h == 0:
        # model.addConstr(beta[a, b, c, d, e, h] == 2*gamma[a, b, c, d, e, h] - alpha00[a, e] - alpha10[c, d,
        # e] + 1) else: model.addConstr(beta[a, b, c, d, e, h] == 2*gamma[a, b, c, d, e, h] - alpha_1[a, c, d, e,
        # 0] - alpha_1[a, c, d, e, 1] + 1)

        # model.addConstrs(beta0[a, c, d, e] == 2*alpha00[a, e]*alpha10[c, d, e] - alpha00[a, e] - alpha10[c, d,
        # e] + 1 for a, c, d, e in beta0.keys()) model.addConstrs(beta1[a, c, d, e] == 2*alpha_1[a, c, d, e, 0]*alpha_1[
        # a, c, d, e, 1] - alpha_1[a, c, d, e, 0] - alpha_1[a, c, d, e, 1] + 1 for a, c, d, e in beta1.keys())

        model.addConstrs(
            beta0[a, c, d, e] == 2 * gamma[a, 0, c, d, e, 0] - alpha00[a, e] - alpha10[c, d, e] + 1 for a, c, d, e in
            beta0.keys())
        model.addConstrs(
            beta1[a, c, d, e] == 2 * gamma[a, 0, c, d, e, 1] - alpha_1[a, c, d, e, 0] - alpha_1[a, c, d, e, 1] + 1 for
            a, c, d, e in beta1.keys())

        # gamma-C
        for a, b, c, d, e, h in gamma.keys():
            if h == 0:
                model.addConstr(gamma[a, b, c, d, e, h] <= alpha00[a, e])
                model.addConstr(gamma[a, b, c, d, e, h] <= alpha10[c, d, e])
                model.addConstr(gamma[a, b, c, d, e, h] >= alpha00[a, e] + alpha10[c, d, e] - 1)

            if h == 1:
                model.addConstr(gamma[a, b, c, d, e, h] <= alpha_1[a, c, d, e, 0])
                model.addConstr(gamma[a, b, c, d, e, h] <= alpha_1[a, c, d, e, 1])
                model.addConstr(gamma[a, b, c, d, e, h] >= alpha_1[a, c, d, e, 0] + alpha_1[a, c, d, e, 1] - 1)

                # delta-C
        for a, c, d, e in beta0.keys():
            model.addConstr((beta0[a, c, d, e] + beta1[a, c, d, e] <= 2 * delta[a, 0, c, d, e]))
            model.addConstr(delta[a, 0, c, d, e] <= 2 * (beta0[a, c, d, e] + beta1[a, c, d, e]))

        # epsilon-C and y-C
        for a, b, c, d in epsilon.keys():
            # print((a, b, c, d))
            model.addConstr((delta.sum(a, b, c, d, '*') - len(barriers)) + 1 <= epsilon[a, b, c, d])
            model.addConstr(len(barriers) * epsilon[a, b, c, d] <= delta.sum(a, b, c, d, '*'))

            # if a < 0 and b >= 0 and c >= 0:
            model.addConstr(y[a, b, c, d] + y[c, d, a, b] <= 2 * epsilon[a, b, c, d])

            # if a < 0 and b < 0 and c == 0:
            #     model.addConstr(y[a, b, c] + y[b, a, c] <= 2*epsilon[a, b, c])

            # if no_ve(neighborhoods[abs(a)-1], barriers[b][0]) and no_ve(neighborhoods[abs(a)-1], barriers[b][1]):
            #     model.addConstr(epsilon[a, b, c] <= 0)

        # NS and NT constraints
        for a, b, dim in point.keys():
            neighborhood = neighborhoods[abs(a) - 1]

            if type(neighborhood) is neigh.Circle:
                model.addConstr(dif_inside[a, b, dim] >= point[a, b, dim] - neighborhoods[abs(a) - 1].center[dim])
                model.addConstr(dif_inside[a, b, dim] >= neighborhoods[abs(a) - 1].center[dim] - point[a, b, dim])

                model.addConstr(
                    gp.quicksum(dif_inside[a, b, dim] * dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b] *
                    d_inside[a, b])
                model.addConstr(d_inside[a, b] <= neighborhoods[abs(a) - 1].radii)

            if type(neighborhood) is neigh.Poligonal:
                model.addConstrs(
                    point[a, b, dim] == landa[a, b] * neighborhood.vertices_total[0][dim] + (1 - landa[a, b]) *
                    neighborhood.vertices_total[1][dim] for
                    dim in range(2))

        # model.addConstr(y[-1, 2, 1] >= 0.5)
        # model.addConstr(y[-2, 2, 1] >= 0.5)
        # model.addConstr(point[-1, 0] == 30)
        # model.addConstr(point[-1, 1] == 10)

        # model.addConstr(epsilon[-4, 3, 1] >= 0.5)
        # model.addConstr(epsilon[3, 1, -4] >= 0.5)

        # model.addConstr(y[3, 1, -1] >= 0.5)

        # dist constraints
        for a, b, c, d, dim in dif_index:
            if (a, b, c, d) in edges_barrier:
                model.addConstr(dist[a, b, c, d] == np.linalg.norm(np.array(barriers[a][b]) - np.array(barriers[c][d])))

            if (a, b, c, d) in edges_neighborhood:
                if a < 0:
                    model.addConstr(dif[a, b, c, d, dim] >= point[a, b, dim] - barriers[c][d][dim])
                    model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + barriers[c][d][dim])
                    model.addConstr(
                        gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] *
                        dist[a, b, c, d])
                if c < 0:
                    model.addConstr(dif[a, b, c, d, dim] >= point[c, d, dim] - barriers[a][b][dim])
                    model.addConstr(dif[a, b, c, d, dim] >= - point[c, d, dim] + barriers[a][b][dim])
                    model.addConstr(
                        gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] *
                        dist[a, b, c, d])

            # if (a, b, c, d) in ENN: model.addConstr(dif[a, b, c, d, dim] >=   point[a, b, dim] - point[c, d,
            # dim]) model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + point[c, d, dim]) model.addConstr(
            # gp.quicksum(dif[a, b, c, d, dim]*dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d]*dist[a, b,
            # c, d])

        l_out = 0
        u_out = 10000

        # p constraints
        for a, b, c, d in p.keys():

            if a < 0:
                neighborhood = neighborhoods[abs(a) - 1]
                punto = barriers[c][d]
                l_out = af.estima_L(neighborhood, punto)
                u_out = af.estima_U(neighborhood, punto)

            if c < 0:
                neighborhood = neighborhoods[abs(c) - 1]
                punto = barriers[a][b]
                l_out = af.estima_L(neighborhood, punto)
                u_out = af.estima_U(neighborhood, punto)

            # print((l_out, u_out))
            model.addConstr(p[a, b, c, d] >= l_out * y[a, b, c, d])
            model.addConstr(p[a, b, c, d] >= dist[a, b, c, d] - u_out * (1 - y[a, b, c, d]))

            # model.addConstr(p[a, b, c, d] <= dist[a, b, c, d]* u_out)

        # model.addConstrs(z[v, v] == 0 for v in vertices_total

        # Restriccion 1
        model.addConstrs(gp.quicksum(y[v, i, vertices_neighborhood, j] for v, i in vertices_total if
                                     (v, i, vertices_neighborhood, j) in edges_neighborhood) >= 1 for
                         vertices_neighborhood, j in vertices_neighborhood)

        # Restriccion 2
        for v, i in vertices_total:
            model.addConstr(gp.quicksum(y[v, i, vertices_neighborhood, j] for vertices_neighborhood, j in vertices_total if
                                        (v, i, vertices_neighborhood, j) in edges_total) == gp.quicksum(
                y[vertices_neighborhood, j, v, i] for vertices_neighborhood, j in vertices_total if
                (vertices_neighborhood, j, v, i) in edges_total))

        # Restriccion 3
        for vertices_neighborhood, i in vertices_neighborhood:
            if vertices_neighborhood <= -2:
                model.addConstr(gp.quicksum(g[vertices_neighborhood, i, v, j] for v, j in vertices_total if
                                            (vertices_neighborhood, i, v, j) in edges_neighborhood) - gp.quicksum(
                    g[v, j, vertices_neighborhood, i] for v, j in vertices_total if
                    (v, j, vertices_neighborhood, i) in edges_neighborhood) == 1)

        # Restriccion 4
        for vertices_barrier, i in vertices_barrier:
            model.addConstr(gp.quicksum(g[(w, j, vertices_barrier, i)] for w, j in vertices_total if
                                        (w, j, vertices_barrier, i) in edges_total) - gp.quicksum(
                g[(vertices_barrier, i, w, j)] for w, j in vertices_total if
                (vertices_barrier, i, w, j) in edges_total) == 0)

        # Restriccion 5
        model.addConstrs(g[a, b, c, d] <= (len(neighborhoods) - 1) * y[a, b, c, d] for a, b, c, d in g.keys())
        # model.addConstrs(gp.quicksum(y[v, i, vertices_neighborhood] for v, i in vertices_barrier) == 1 for
        # vertices_neighborhood in vertices_neighborhood) model.addConstrs(gp.quicksum(y[vertices_neighborhood, v,
        # i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in vertices_neighborhood)
        #
        # model.addConstrs(gp.quicksum(g[v, i, vertices_neighborhood] for v, i in vertices_barrier) - gp.quicksum(g[
        # vertices_neighborhood, v, i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in
        # vertices_neighborhood) model.addConstrs(gp.quicksum(g[index] for index in edges_barrier if index[2] == v and
        # index[3] == i) + gp.quicksum(g[index] for index in edges_neighborhood if index[1] == v and index[2] == i) -
        # gp.quicksum(g[index] for index in edges_barrier if index[0] == v and index[1] == i) - gp.quicksum(g[index] for
        # index in edges_neighborhood if index[0] == v and index[1] == i) == 0 for v, i in vertices_barrier)

        # model.addConstrs(gp.quicksum(f[-1, w, k] for w in vertices_neighborhood if w <= -2) == 1 for k in
        # vertices_neighborhood if k <= -2) model.addConstrs(gp.quicksum(f[v, w, w] ))

        # model.addConstrs(gp.quicksum(z[w, v] for v in vertices_neighborhood if w != v) == 1 for w in
        # vertices_neighborhood)

        # flow conservation constraints for index in y_index: if len(index) == 3: model.addConstrs(gp.quicksum(y[tupla]
        # for tupla in edges_neighborhood if tupla[0] == v) == 1 for v in vertices_neighborhood)
        #
        # for v, i in vertices_barrier: tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[0] == v
        # and tupla[1] == i]) + gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v and tupla[1] ==
        # i]) tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[2] == v and tupla[3] == i]) +
        # gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[1] == v and tupla[2] == i])
        #
        #     model.addConstr(tuplas_salen - tuplas_entran == 0)
        #
        # for v in vertices_neighborhood:
        #     tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v])
        #     tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[2] == v])
        #
        #     model.addConstr(tuplas_salen - tuplas_entran == 0)
        #
        # model.addConstrs(gp.quicksum(y[tupla] for tupla in edges_neighborhood if tupla[2] == w) == 1 for w in
        # vertices_neighborhood)

        model.update()

        objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(
            dist[index] * y[index] for index in edges_barrier)
        model.setObjective(objective, GRB.MINIMIZE)

        second_time = time.time()

        time_elapsed = second_time - first_time

        model.update()

        model.Params.Threads = 6
        model.Params.timeLimit = time_limit - time_elapsed
        # model.Params.LazyConstraints = 1
        # model.Params.NumericFocus = 1
        # model.Params.NonConvex = 2

        model.write('prueba.lp')
        model.write('prueba.mps')

        model.optimize()

        results = [len(neighborhoods), len(barriers), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

        if init:
            try:
                results[-2] = time_h
                results[-1] = objval_h
            except:
                print('No solution obtained by the heuristic')

        if model.Status == 3:
            model.computeIIS()
            model.write('infeasible_constraints.ilp')
            return results

        if model.SolCount == 0:
            return results

        model.write('solution.sol')

        results[2] = model.getAttr('MIPGap')
        results[3] = model.Runtime + time_elapsed
        results[4] = model.getAttr('NodeCount')
        results[5] = model.ObjVal

        y_indices = []

        for index in edges_total:
            if y[index].X > 0.5:
                y_indices.append(index)

        if log:
            print(y_indices)

        g_indices = []

        for index in g_index:
            if g[index].X > 0.5:
                g_indices.append(g[index])

        if log:
            print(g_indices)

        if picture:
            fig, ax = plt.subplots()

            for b in barriers:
                ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')

            for n in neighborhoods:
                ax.add_artist(n.artist)

            p_vals = model.getAttr('x', point)
            print(p_vals)

            points = []
            for keys, vals in p_vals.items():
                points.append(vals)

            points = np.array(points).reshape((len(neighborhoods), 2))
            print(points)

            for i in points:
                ax.scatter(i[0], i[1], s=10, c='black')

            # print(points)

            segments = []

            for a, b, c, d in y_indices:
                if (a, b, c, d) in edges_neighborhood:
                    if a < 0:
                        segments.append(
                            [points[abs(a) - 1][0], barriers[c][d][0], points[abs(a) - 1][1], barriers[c][d][1]])
                    if c < 0:
                        segments.append(
                            [barriers[a][b][0], points[abs(c) - 1][0], barriers[a][b][1], points[abs(c) - 1][1]])
                if (a, b, c, d) in edges_barrier:
                    segments.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])

                # if (a, b, c, d) in ENN: segments.append([points[abs(a)-1][0], points[abs(c)-1][0], points[abs(a)-1][
                # 1], points[abs(c)-1][1]])

            # print(segments)
            for segment in segments:
                ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
                         head_width=1, length_includes_head=True, color='black')

            # plt.axis([-5, 105, -5, 105])
            plt.axis([0, 100, 0, 100])

            ax.set_aspect('equal')
            plt.show()
    else:
        first_time = time.time()

        vertices_neighborhood = list(itertools.product(range(-len(neighborhoods), 0), range(1)))
        vertices_neighborhood = vertices_neighborhood[::-1]

        edges_neighborhood = []

        for (a, b) in vertices_neighborhood:
            for c in range(len(barriers)):
                for d in range(2):
                    if prepro:
                        # Point of the barrier to check if is visible by the neighborhood
                        point = barriers[c][d]

                        # Neighborhood to check if it is visible by the point
                        neighborhood = neighborhoods[abs(a) - 1]

                        if af.cansee(point, neighborhood, barriers):
                            # Appending the feasible edges to edges_neighborhood
                            edges_neighborhood.append((a, b, c, d))
                            edges_neighborhood.append((c, d, a, b))
                    else:
                        edges_neighborhood.append((a, b, c, d))
                        edges_neighborhood.append((c, d, a, b))

        vertices_barrier = list(itertools.product(range(len(barriers)), range(2)))

        edges_barrier = []
        for v, i in vertices_barrier:
            for w, j in vertices_barrier:
                if v != w:
                    if prepro:
                        barrier = [barriers[v][i], barriers[w][j]]

                        if any([not (af.intersect(barrieri, barrier)) for barrieri in barriers]):
                            edges_barrier.append((v, i, w, j))
                    else:
                        edges_barrier.append((v, i, w, j))

        vertices_total = vertices_neighborhood + vertices_barrier
        edges_total = edges_neighborhood + edges_barrier

        if log:
            print("vertices_neighborhood = " + str(vertices_neighborhood))
            print("vertices_barrier = " + str(vertices_barrier))

            print("edges_neighborhood = " + str(edges_neighborhood))
            print("edges_barrier = " + str(edges_barrier))

        point_index = []
        for a, b in vertices_neighborhood:
            for dim in range(2):
                point_index.append((a, b, dim))

        y_index = edges_total

        if log:
            print("y = " + str(y_index))

        dist_index = edges_total

        if log:
            print("dist = " + str(dist_index))

        dif_index = []

        for a, b, c, d in edges_total:
            for dim in range(2):
                dif_index.append((a, b, c, d, dim))

        # P_S and P_T: indices of the points in the neighborhoods
        p_index = edges_neighborhood

        if log:
            print("p_index = " + str(p_index))

        dom_set = af.dominant_set(neighborhoods, barriers)

        # z variables:
        z_index = list(dom_set.keys())
        print("z_index = " + str(z_index))

        # f variables:
        g_index = y_index

        if log:
            print("g_index = " + str(g_index))


        model = gp.Model('HTSPN_Model_dominant')

        p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
        y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
        dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')
        dif = model.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')

        # point = model.addVars(p_index, vtype = GRB.CONTINUOUS, lb = 0.1, ub = 99.9, name = 'point')
        point = model.addVars(point_index, vtype=GRB.CONTINUOUS, name='point')

        z = model.addVars(z_index, vtype = GRB.BINARY, name = 'z')
        g = model.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name='g')

        model.update()

        if init:
            time_h, objval_h = heuristic(barriers, neighborhoods)

        # NS and NT constraints
        for a, b, dim in point.keys():

            model.addConstr(point[a, b, dim] == gp.quicksum(dom_set[index][dim]*z[index] for index in dom_set.keys() if index[0] == a))


        for a, c, d, e, f in z.keys():
            model.addConstr(z[a, c, d, e, f] <= y[c, d, a, 0])
            model.addConstr(z[a, c, d, e, f] <= y[a, 0, e, f])
            model.addConstr(z[a, c, d, e, f] >= y[c, d, a, 0] + y[a, 0, e, f] - 1)

        # z constraints
        # for
        # model.addConstr(y[-1, 2, 1] >= 0.5)
        # model.addConstr(y[-2, 2, 1] >= 0.5)
        # model.addConstr(point[-1, 0] == 30)
        # model.addConstr(point[-1, 1] == 10)

        # model.addConstr(epsilon[-4, 3, 1] >= 0.5)
        # model.addConstr(epsilon[3, 1, -4] >= 0.5)

        # model.addConstr(y[3, 1, -1] >= 0.5)

        # dist constraints
        for a, b, c, d, dim in dif_index:
            if (a, b, c, d) in edges_barrier:
                model.addConstr(dist[a, b, c, d] == np.linalg.norm(np.array(barriers[a][b]) - np.array(barriers[c][d])))

            if (a, b, c, d) in edges_neighborhood:
                if a < 0:
                    model.addConstr(dif[a, b, c, d, dim] >= point[a, b, dim] - barriers[c][d][dim])
                    model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + barriers[c][d][dim])
                    model.addConstr(
                        gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] *
                        dist[a, b, c, d])
                if c < 0:
                    model.addConstr(dif[a, b, c, d, dim] >= point[c, d, dim] - barriers[a][b][dim])
                    model.addConstr(dif[a, b, c, d, dim] >= - point[c, d, dim] + barriers[a][b][dim])
                    model.addConstr(
                        gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] *
                        dist[a, b, c, d])

            # if (a, b, c, d) in ENN: model.addConstr(dif[a, b, c, d, dim] >=   point[a, b, dim] - point[c, d,
            # dim]) model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + point[c, d, dim]) model.addConstr(
            # gp.quicksum(dif[a, b, c, d, dim]*dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d]*dist[a, b,
            # c, d])

        l_out = 0
        u_out = 10000

        # p constraints
        for a, b, c, d in p.keys():

            if a < 0:
                neighborhood = neighborhoods[abs(a) - 1]
                punto = barriers[c][d]
                l_out = af.estima_L(neighborhood, punto)
                u_out = af.estima_U(neighborhood, punto)

            if c < 0:
                neighborhood = neighborhoods[abs(c) - 1]
                punto = barriers[a][b]
                l_out = af.estima_L(neighborhood, punto)
                u_out = af.estima_U(neighborhood, punto)

            # print((l_out, u_out))
            model.addConstr(p[a, b, c, d] >= l_out * y[a, b, c, d])
            model.addConstr(p[a, b, c, d] >= dist[a, b, c, d] - u_out * (1 - y[a, b, c, d]))

            # model.addConstr(p[a, b, c, d] <= dist[a, b, c, d]* u_out)

        # model.addConstrs(z[v, v] == 0 for v in vertices_total

        # Restriccion 1
        model.addConstrs(gp.quicksum(y[v, i, vertices_neighborhood, j] for v, i in vertices_total if
                                     (v, i, vertices_neighborhood, j) in edges_neighborhood) >= 1 for
                         vertices_neighborhood, j in vertices_neighborhood)

        # Restriccion 2
        for v, i in vertices_total:
            model.addConstr(gp.quicksum(y[v, i, vertices_neighborhood, j] for vertices_neighborhood, j in vertices_total if
                                        (v, i, vertices_neighborhood, j) in edges_total) == gp.quicksum(
                y[vertices_neighborhood, j, v, i] for vertices_neighborhood, j in vertices_total if
                (vertices_neighborhood, j, v, i) in edges_total))

        # Restriccion 3
        for vertices_neighborhood, i in vertices_neighborhood:
            if vertices_neighborhood <= -2:
                model.addConstr(gp.quicksum(g[vertices_neighborhood, i, v, j] for v, j in vertices_total if
                                            (vertices_neighborhood, i, v, j) in edges_neighborhood) - gp.quicksum(
                    g[v, j, vertices_neighborhood, i] for v, j in vertices_total if
                    (v, j, vertices_neighborhood, i) in edges_neighborhood) == 1)

        # Restriccion 4
        for vertices_barrier, i in vertices_barrier:
            model.addConstr(gp.quicksum(g[(w, j, vertices_barrier, i)] for w, j in vertices_total if
                                        (w, j, vertices_barrier, i) in edges_total) - gp.quicksum(
                g[(vertices_barrier, i, w, j)] for w, j in vertices_total if
                (vertices_barrier, i, w, j) in edges_total) == 0)

        # Restriccion 5
        model.addConstrs(g[a, b, c, d] <= (len(neighborhoods) - 1) * y[a, b, c, d] for a, b, c, d in g.keys())
        # model.addConstrs(gp.quicksum(y[v, i, vertices_neighborhood] for v, i in vertices_barrier) == 1 for
        # vertices_neighborhood in vertices_neighborhood) model.addConstrs(gp.quicksum(y[vertices_neighborhood, v,
        # i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in vertices_neighborhood)
        #
        # model.addConstrs(gp.quicksum(g[v, i, vertices_neighborhood] for v, i in vertices_barrier) - gp.quicksum(g[
        # vertices_neighborhood, v, i] for v, i in vertices_barrier) == 1 for vertices_neighborhood in
        # vertices_neighborhood) model.addConstrs(gp.quicksum(g[index] for index in edges_barrier if index[2] == v and
        # index[3] == i) + gp.quicksum(g[index] for index in edges_neighborhood if index[1] == v and index[2] == i) -
        # gp.quicksum(g[index] for index in edges_barrier if index[0] == v and index[1] == i) - gp.quicksum(g[index] for
        # index in edges_neighborhood if index[0] == v and index[1] == i) == 0 for v, i in vertices_barrier)

        # model.addConstrs(gp.quicksum(f[-1, w, k] for w in vertices_neighborhood if w <= -2) == 1 for k in
        # vertices_neighborhood if k <= -2) model.addConstrs(gp.quicksum(f[v, w, w] ))

        # model.addConstrs(gp.quicksum(z[w, v] for v in vertices_neighborhood if w != v) == 1 for w in
        # vertices_neighborhood)

        # flow conservation constraints for index in y_index: if len(index) == 3: model.addConstrs(gp.quicksum(y[tupla]
        # for tupla in edges_neighborhood if tupla[0] == v) == 1 for v in vertices_neighborhood)
        #
        # for v, i in vertices_barrier: tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[0] == v
        # and tupla[1] == i]) + gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v and tupla[1] ==
        # i]) tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_barrier if tupla[2] == v and tupla[3] == i]) +
        # gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[1] == v and tupla[2] == i])
        #
        #     model.addConstr(tuplas_salen - tuplas_entran == 0)
        #
        # for v in vertices_neighborhood:
        #     tuplas_salen = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[0] == v])
        #     tuplas_entran = gp.quicksum([y[tupla] for tupla in edges_neighborhood if tupla[2] == v])
        #
        #     model.addConstr(tuplas_salen - tuplas_entran == 0)
        #
        # model.addConstrs(gp.quicksum(y[tupla] for tupla in edges_neighborhood if tupla[2] == w) == 1 for w in
        # vertices_neighborhood)

        model.update()

        objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(
            dist[index] * y[index] for index in edges_barrier)
        model.setObjective(objective, GRB.MINIMIZE)

        second_time = time.time()

        time_elapsed = second_time - first_time

        model.update()

        model.Params.Threads = 6
        model.Params.timeLimit = time_limit - time_elapsed
        # model.Params.LazyConstraints = 1
        model.Params.NumericFocus = 1
        # model.Params.NonConvex = 2

        model.write('prueba.lp')
        model.write('prueba.mps')

        model.optimize()

        results = [len(neighborhoods), len(barriers), np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

        if init:
            try:
                results[-2] = time_h
                results[-1] = objval_h
            except:
                print('No solution obtained by the heuristic')



        if model.SolCount == 0:
            return results

        model.write('solution.sol')

        results[2] = model.getAttr('MIPGap')
        results[3] = model.Runtime + time_elapsed
        results[4] = model.getAttr('NodeCount')
        results[5] = model.ObjVal

        y_indices = []

        for index in edges_total:
            if y[index].X > 0.5:
                y_indices.append(index)

        if log:
            print(y_indices)

        g_indices = []

        for index in g_index:
            if g[index].X > 0.5:
                g_indices.append(g[index])

        if log:
            print(g_indices)

        if picture:
            fig, ax = plt.subplots()

            for b in barriers:
                ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')

            for n in neighborhoods:
                ax.add_artist(n.artist)

            p_vals = model.getAttr('x', point)
            print(p_vals)

            points = []
            for keys, vals in p_vals.items():
                points.append(vals)

            points = np.array(points).reshape((len(neighborhoods), 2))
            print(points)

            for i in points:
                ax.scatter(i[0], i[1], s=10, c='black')

            # print(points)

            segments = []

            for a, b, c, d in y_indices:
                if (a, b, c, d) in edges_neighborhood:
                    if a < 0:
                        segments.append(
                            [points[abs(a) - 1][0], barriers[c][d][0], points[abs(a) - 1][1], barriers[c][d][1]])
                    if c < 0:
                        segments.append(
                            [barriers[a][b][0], points[abs(c) - 1][0], barriers[a][b][1], points[abs(c) - 1][1]])
                if (a, b, c, d) in edges_barrier:
                    segments.append([barriers[a][b][0], barriers[c][d][0], barriers[a][b][1], barriers[c][d][1]])

                # if (a, b, c, d) in ENN: segments.append([points[abs(a)-1][0], points[abs(c)-1][0], points[abs(a)-1][
                # 1], points[abs(c)-1][1]])

            # print(segments)
            for segment in segments:
                ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
                         head_width=1, length_includes_head=True, color='black')

            # plt.axis([-5, 105, -5, 105])
            plt.axis([0, 100, 0, 100])

            ax.set_aspect('equal')
            plt.show()

        pass

    return results
