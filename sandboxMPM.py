# -*- coding: utf-8 -*-
# @Author: chavez
# @Date:   2018-08-20 18:36:09
# @Last Modified by:   chavez
# @Last Modified time: 2018-09-13 15:17:55
import math
import numpy as np
from scipy.linalg import polar
n = 64
dt = 1
frame_dt = 1e-3
dx = 1.0 / n
inv_dx = 1.0 / dx
particle_mass = 1.0
vol = 1.0
hardening = 10.0
E = 1e4
nu = 0.2
mu_0 = E / (2 * (1 + nu))
lambda_0 = E * nu / ((1 + nu) * (1 - 2 * nu))

# class Particle(object):
#   def __init__(self, arg):
#       self.arg = arg

# def advance(dt):
#     for p in particles:
#         base_coord = (p.x*inv_dx - )
#         fx = p.x * inv_dx - base_coord

#         w = []
#         e = math.exp(hardening * (1.0 - p.Jp))
#         mu = mu_0 * e
#         lambdaA = lambda_0 * e
#         J = determinant(p.F)
#         r, s = polar_decomp(p.F, r, s)
#         stress = -4*inv_dx*inv_dx*dt*vol*(2*mu*(p.F-r) * transposed(p.F)+lambdaA*(J-1)*J);
#         for i in range(0, 3):
#             for j in range(0, 3): 
#                 dpos = ([i,j] - fx) * dx
#                 mv = (p.v * particle_mass, particle_mass)
#                 grid[base_coord.x + i][base_coord + j] +=
#                 w[i].x * w[j].y * (mv + [stress+particle_mass*p.C] *dpos,0)

#         for i in range(0, n):
#             for j in range(0, n):
#                 g = grid[i][j]
#                 if g[2] > 0{
#                     g /= g[2]
#                     g += dt * [0, -100, 0]
#                     boundary = 0.05
#                     x = i/n
#                     y = j/n
#                     if (x < boundary || x > 1-boundary || y > 1-boundary)
#                         g = [0,0,0]
#                     if (y < boundary)
#                         g[1] = math.max(0.0, g[1])

#         for p in particles:
#             base_coord = (p.x * inv_dx - [0.5, 0.5])
#             fx = p.x * inv_dx - base_coord
#             w = [Vec(0.5) * sqr(Vec(1.5) - fx), Vec(0.75) - sqrt(fx - Vec(1.0)),
#              Vec(0.5) * sqr(fx - Vec(0.5))]

#             p.C = Mat(0)
#             p.v = Vec(0)

#             for i in range(0, 3):
#                 for j in range(0, 3):
#                     dpos = [i,j] - fx
#                     grid_v = grid[base_coord.x + i][base_coord.y + i]
#                     weight = w[i].x * w[j].y
#                     p.v += weight * grid_v
#                     p.C += 4 * inv_dx * Mat.outer_prodocut(weight * grid_v, dpos)

#             p.x += dt * p.v
#             F = (Mat(1) + dt * p.C) * p.F
#             Mat svd_u, sig, svd_v
#             svd(F, svd_u, sig, svd_v)
#             for i in range(0, 2):
#                 sig[i][j] = clamp(sig[i][i], 1.0 - 2.5e-2, 1.0 + 7.5e-3)
#             oldJ = determinant(F)
#             F = svd_u * sig * transposed(svd_v)
#             Jp_new = clamp(p.Jp * oldJ / determinant(F), 0.6, 20.0)
#             p.Jp  = Jp_new
#             p.F = F

#                 }
#                 

grid = np.zeros(shape=(64, 64))    
def advance(particles):
    for p in np.nditer(particles):
        base_coord = np.array(p['pos']) #* inv_dx
        fx = p['pos'] * inv_dx - base_coord
        e = math.exp(hardening * (1.0 - p['Jp']))
        mu=mu_0 * e
        J = np.linalg.det(p['F'])
        r,s = polar(p['F'])
        stress = -4 *  inv_dx * inv_dx * dt*vol*(2*mu*(p['F']-r) * np.transpose(p['F'])+ (E * nu / ((1 + nu) * (1 - 2 * nu)) * (J-1)*J ))
        #w = np.array([0.5] * math.sqrt())
        #w = np.array(np.array([0.5, 0.5]) * np.sqrt(np.array([0.5, 0.5] - fx)))
        #print(p['pos'])
       # print fx - np.array([0.5, 0.5])
        #print(np.sqrt(fx - np.array([0.5, 0.5])))
        print (fx)

        #print (fx - np.array([1.0, 1.0]))
        #exit()
        #print(np.sqrt(np.array(fx - np.array([1.0, 1.0]))))
        # w = np.array([np.array([0.5, 0.5]) * np.sqrt(np.array([1.5, 1.5] - fx)),
        #              np.array([0.75, 0.75]) - np.sqrt((fx - np.array([1.0, 1.0])))#, 
        #             #np.array([0.5, 0.5]) * np.sqrt(np.array(fx - np.array([0.5, 0.5])))
        #     ])
        # for i in range(0, 3):
        #     for j in range(0, 3):
        #         dpos = ([i,j] - fx) * dx
        #         #weak
        #         mv = [p['v'] * particle_mass, particle_mass]
        #         #grid[base_coord[0] + i][base_coord[1] + j] += w[i][0]*w[j][1] * (mv + 

numOfParticles = 5
def add_object(x,y):
    d = np.dtype([('pos', float, 2), ('v', float, 2), ('F', float, (2, 2)),('C', float, (2, 2)), ('Jp', float)])
    #arr = np.array(1000, dtype=d).view(np.recarray)
    arr = np.rec.array(np.zeros(numOfParticles, dtype=d)) 
    arr.pos = [[(np.random.random_sample() * 2.0-1) * 0.08 + x, (np.random.random_sample() * 2.0-1) * 0.08 + y] for i in range(0, numOfParticles)]
    #print(arr.pos)
    return arr

def main():
    particles = add_object(0,0)
    advance(particles)

main()
#struct Particle { Vec x, v; Mat F, C; real Jp;
#  Particle(Vec x, Vec v=Vec(0)) : x(x), v(v), F(1), C(0), Jp(1) {}