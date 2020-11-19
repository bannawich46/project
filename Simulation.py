import anuga
import numpy
import math
#import matplotlib
from anuga.structures import internal_boundary_functions
from anuga.structures.internal_boundary_functions import pumping_station_function
from anuga.geometry.polygon_function import Polygon_function
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.shallow_water.shallow_water_domain import Domain
from anuga import distribute, myid, numprocs, finalize, barrier
#import matplotlib.pyplot as plt

import zmq
import pickle
import time

context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:7256")

scale = 1.0
length = 306.0 * scale
width = 20.0 * scale
dx = dy = 1.0
waterset = -100

pool = scale * 20.0
px1 = scale * 0.0
py1 = scale * 39.0
px2 = scale * 23.0
py2 = scale * 26.0
px3 = scale * 46.0
py3 = scale * 19.0
px4 = scale * 70.0
py4 = scale * 19.0
px5 = scale * 78.4
py5 = scale * 15.0
px6 = scale * 86.8
py6 = scale * 12.7
px7 = scale * 95.2
py7 = scale * 11.9
px8 = scale * 112.0
py8 = scale * 7.0
px9 = scale * 146.0
py9 = scale * 4.0
px10 = scale * 168.0
py10 = scale * 4.0
px11 = scale * 173.5
py11 = scale * 3.0
px12 = scale * 184.5
py12 = scale * 3.0
px13 = scale * 200.0
py13 = scale * 1.0
px14 = scale * 280.8
py14 = scale * 0.0
px15 = scale * 286.0
py15 = scale * -1.0
m1 = (py2 - py1)/(px2 - px1)
m2 = (py3 - py2)/(px3 - px2)
m3 = (py4 - py3)/(px4 - px3)
m4 = (py5 - py4)/(px5 - px4)
m5 = (py6 - py5)/(px6 - px5)
m6 = (py7 - py6)/(px7 - px6)
m7 = (py8 - py7)/(px8 - px7)
m8 = (py9 - py8)/(px9 - px8)
m9 = (py10 - py9)/(px10 - px9)
m10 = (py11 - py10)/(px11 - px10)
m11 = (py12 - py11)/(px12 - px11)
m12 = (py13 - py12)/(px13 - px12)
m13 = (py14 - py13)/(px14 - px13)
m14 = (py15 - py14)/(px15 - px14)

wallposition1 = pool
wallposition2 = px10
wallheight1 = 60.0
wallheight2 = 60.0
gap = 1.0
endpoint1 = [[wallposition1 - (gap * 2),width/2],[wallposition1,width/2]]
endpoint2 = [[wallposition2 - (gap * 2),width/2],[wallposition2,width/2]]
duration2 = 2

points,vertices,boundary = anuga.rectangular_cross(int(length/dx),int(width/dy),len1=length,len2=width)

riverWall = {'centralWall':[[wallposition1 - gap,0.0,wallheight1],[wallposition1 - gap,width,wallheight1]
                            ,[wallposition2 - gap,width,wallheight2],[wallposition2 - gap,0.0,wallheight2]]}

def topography(x,y):
    z = x
    N = len(x)
    for i in range(N):
        if 0 < x[i] <= pool:
            z[i] = 0
        if pool + px1 < x[i] <= pool + px2:
            z[i] = (x[i]*m1) - py1 - (m1*(pool + px2)) + py2
        if pool + px2 < x[i] <= pool + px3:
            z[i] = (x[i]*m2) - py1 - (m2*(pool + px3)) + py3
        if pool + px3 < x[i] <= pool + px4:
            z[i] = (x[i]*m3) - py1 - (m3*(pool + px4)) + py4
        if pool + px4 < x[i] <= pool + px5:
            z[i] = (x[i]*m4) - py1 - (m4*(pool + px5)) + py5
        if pool + px5 < x[i] <= pool + px6:
            z[i] = (x[i]*m5) - py1 - (m5*(pool + px6)) + py6
        if pool + px6 < x[i] <= pool + px7:
            z[i] = (x[i]*m6) - py1 - (m6*(pool + px7)) + py7
        if pool + px7 < x[i] <= pool + px8:
            z[i] = (x[i]*m7) - py1 - (m7*(pool + px8)) + py8
        if pool + px8 < x[i] <= pool + px9:
            z[i] = (x[i]*m8) - py1 - (m8*(pool + px9)) + py9
        if pool + px9 < x[i] <= pool + px10:
            z[i] = (x[i]*m9) - py1 - (m9*(pool + px10)) + py10
        if pool + px10 < x[i] <= pool + px11:
            z[i] = (x[i]*m10) - py1 - (m10*(pool + px11)) + py11
        if pool + px11 < x[i] <= pool + px12:
            z[i] = (x[i]*m11) - py1 - (m11*(pool + px12)) + py12
        if pool + px12 < x[i] <= pool + px13:
            z[i] = (x[i]*m12) - py1 - (m12*(pool + px13)) + py13
        if pool + px13 < x[i] <= pool + px14:
            z[i] = (x[i]*m13) - py1 - (m13*(pool + px14)) + py14
        if pool + px14 < x[i] <= pool + px15:
            z[i] = (x[i]*m14) - py1 - (m14*(pool + px15)) + py15
    return z

def stagefun(x, y):
    stg = x
    M = len(x)
    scale = 0.4
    pool = scale * 20.0
    px1 = scale * 0.0
    py1 = scale * 39.0
    px2 = scale * 23.0
    py2 = scale * 26.0
    px3 = scale * 46.0
    py3 = scale * 19.0
    px4 = scale * 70.0
    py4 = scale * 19.0
    px5 = scale * 78.4
    py5 = scale * 15.0
    px6 = scale * 86.8
    py6 = scale * 12.7
    px7 = scale * 95.2
    py7 = scale * 11.9
    px8 = scale * 112.0
    py8 = scale * 7.0
    px9 = scale * 146.0
    py9 = scale * 4.0
    px10 = scale * 168.0
    py10 = scale * 4.0
    px11 = scale * 173.5
    py11 = scale * 3.0
    px12 = scale * 184.5
    py12 = scale * 3.0
    px13 = scale * 200.0
    py13 = scale * 1.0
    px14 = scale * 280.8
    py14 = scale * 0.0
    px15 = scale * 286.0
    py15 = scale * -1.0
    m1 = (py2 - py1)/(px2 - px1)
    m2 = (py3 - py2)/(px3 - px2)
    m3 = (py4 - py3)/(px4 - px3)
    m4 = (py5 - py4)/(px5 - px4)
    m5 = (py6 - py5)/(px6 - px5)
    m6 = (py7 - py6)/(px7 - px6)
    m7 = (py8 - py7)/(px8 - px7)
    m8 = (py9 - py8)/(px9 - px8)
    m9 = (py10 - py9)/(px10 - px9)
    m10 = (py11 - py10)/(px11 - px10)
    m11 = (py12 - py11)/(px12 - px11)
    m12 = (py13 - py12)/(px13 - px12)
    m13 = (py14 - py13)/(px14 - px13)
    m14 = (py15 - py14)/(px15 - px14)
    for j in range(M):
        if 0 < x[j] <= pool:
            stg[j] = 0 + waterset
        if pool + px1 < x[j] <= pool + px2:
            stg[j] = (x[j]*m1) - py1 - (m1*(pool + px2)) + py2 + waterset
        if pool + px2 < x[j] <= pool + px3:
            stg[j] = (x[j]*m2) - py1 - (m2*(pool + px3)) + py3 + waterset
        if pool + px3 < x[j] <= pool + px4:
            stg[j] = (x[j]*m3) - py1 - (m3*(pool + px4)) + py4 + waterset
        if pool + px4 < x[j] <= pool + px5:
            stg[j] = (x[j]*m4) - py1 - (m4*(pool + px5)) + py5 + waterset
        if pool + px5 < x[j] <= pool + px6:
            stg[j] = (x[j]*m5) - py1 - (m5*(pool + px6)) + py6 + waterset
        if pool + px6 < x[j] <= pool + px7:
            stg[j] = (x[j]*m6) - py1 - (m6*(pool + px7)) + py7 + waterset
        if pool + px7 < x[j] <= pool + px8:
            stg[j] = (x[j]*m7) - py1 - (m7*(pool + px8)) + py8 + waterset
        if pool + px8 < x[j] <= pool + px9:
            stg[j] = (x[j]*m8) - py1 - (m8*(pool + px9)) + py9 + waterset
        if pool + px9 < x[j] <= pool + px10:
            stg[j] = (x[j]*m9) - py1 - (m9*(pool + px10)) + py10 + waterset
        if pool + px10 < x[j] <= pool + px11:
            stg[j] = (x[j]*m10) - py1 - (m10*(pool + px11)) + py11 + waterset
        if pool + px11 < x[j] <= pool + px12:
            stg[j] = (x[j]*m11) - py1 - (m11*(pool + px12)) + py12 + waterset
        if pool + px12 < x[j] <= pool + px13:
            stg[j] = (x[j]*m12) - py1 - (m12*(pool + px13)) + py13 + waterset
        if pool + px13 < x[j] <= pool + px14:
            stg[j] = (x[j]*m13) - py1 - (m13*(pool + px14)) + py14 + waterset
        if pool + px14 < x[j] <= pool + px15:
            stg[j] = (x[j]*m14) - py1 - (m14*(pool + px15)) + py15 + waterset
    return stg

domain = anuga.Domain(points,vertices,boundary)
domain.set_name('Simulation')
domain.set_flow_algorithm('DE0')
domain.set_store_vertices_uniquely(True)

domain.set_quantity('elevation',topography,location = 'centroids')
domain.set_quantity('stage', stagefun,location = 'centroids')

domain.set_quantity('friction',0.001,location = 'centroids')
domain.set_quantity('stage',expression = 'elevation',location = 'centroids')

domain.riverwallData.create_riverwalls(riverWall,verbose = False) 

domain = distribute(domain,verbose = False)
domain.set_store_vertices_uniquely(False)

pump_function1 = anuga.pumping_station_function(domain = domain,
                                                pump_capacity = 1600.0,
                                                hw_to_start_pumping = 0.0,
                                                hw_to_stop_pumping = 0.0,
                                                initial_pump_rate = 1600.0, 
                                                pump_rate_of_increase = 0.0, 
                                                pump_rate_of_decrease = 0.0, 
                                                verbose = False)
pump_function2 = anuga.pumping_station_function(domain = domain,
                                                pump_capacity = 1600.0,
                                                hw_to_start_pumping = 0.0,
                                                hw_to_stop_pumping = 0.0,
                                                initial_pump_rate = 1600.0,
                                                pump_rate_of_increase = 0.0, 
                                                pump_rate_of_decrease = 0.0, 
                                                verbose = False)

pump1 = anuga.Internal_boundary_operator(domain,pump_function1,
                                        width = width,
                                        height = 0.0,
                                        apron = 0.0,
                                        end_points = endpoint1,
                                        verbose = False)
pump2 = anuga.Internal_boundary_operator(domain,pump_function2,
                                        width = width,
                                        height = 0.0,
                                        apron = 0.0,
                                        end_points = endpoint2,
                                        verbose = False)

Bi = anuga.Dirichlet_boundary([0.0,0.0,0.0])
Bo = anuga.Dirichlet_boundary([-100.0,0.0,0.0])
Br = anuga.Reflective_boundary(domain)
domain.set_boundary({'left': Bi,'right': Bo,'top': Br,'bottom': Br})

while True:
    for j in range(60):
        print("Step :", j)
        come1 = socket.recv()
        informationcome1 = pickle.loads(come1)

        if informationcome1 == "close":
            break
        else:
            [open1, open2] = informationcome1

            #print open1, open2
            print("Gate1 :", open1)
            print("Gate2 :", open2)

            pump_function1.pump_rate = open1 * 100
            pump_function2.pump_rate = open2 * 100

        checkstage1 = domain.get_quantity('stage').get_values(interpolation_points = [[wallposition1/2,width/2]])
        checkelevation1 = domain.get_quantity('elevation').get_values(interpolation_points = [[wallposition1/2,width/2]])
        checkstage2 = domain.get_quantity('stage').get_values(interpolation_points = [[(wallposition2 - wallposition1)/2,width/2]])
        checkelevation2 = domain.get_quantity('elevation').get_values(interpolation_points = [[(wallposition2 - wallposition1)/2,width/2]])
        checkstage3 = domain.get_quantity('stage').get_values(interpolation_points = [[(length - wallposition2)/2,width/2]])
        checkelevation3 = domain.get_quantity('elevation').get_values(interpolation_points = [[(length - wallposition2)/2,width/2]])

        check4 = checkstage1 - checkelevation1
        check5 = checkstage2 - checkelevation2
        check6 = checkstage3 - checkelevation3

        check1 = float(numpy.asarray(check4))
        check2 = float(numpy.asarray(check5))
        check3 = float(numpy.asarray(check6))

        mean = (1 + check2 + check3)/3
        SD = math.sqrt((1 - mean)*(1 - mean) + (check2 - mean)*(check2 - mean) + (check3 - mean)*(check3 - mean))/2
        if (SD == 0):
            reward = 0
        else:
            reward = (0.1/SD )
        #if 1.5 < check2 < 2 and 1.5 < check3 < 2:
        #    reward = 1
        #else:  
        #    reward = 0
        
        print("Reward :", reward)
        print check1, check2, check3
        #print reward

        informationgo1 = [reward, check1, check2, check3]

        #print(informationgo1)

        go1 = pickle.dumps(informationgo1)
        socket.send(go1)

        for t in domain.evolve(yieldstep = 0.1,duration = duration2):
            if 0 <= t <= 5 * duration2:
                Bi = anuga.Dirichlet_boundary([11.8988,0.0,0.0])
            if 5 * duration2 < t <= 10 * duration2:
                Bi = anuga.Dirichlet_boundary([0.4124,0.0,0.0])
            if 10 * duration2 < t <= 15 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0636,0.0,0.0])
            if 15 * duration2 < t <= 20 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0000,0.0,0.0])
            if 20 * duration2 < t <= 25 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0003,0.0,0.0])
            if 25 * duration2 < t <= 30 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0000,0.0,0.0])
            if 30 * duration2 < t <= 35 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0000,0.0,0.0])
            if 35 * duration2 < t <= 40 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0000,0.0,0.0])
            if 40 * duration2 < t <= 45 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0000,0.0,0.0])
            if 45 * duration2 < t <= 50 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0052,0.0,0.0])
            if 50 * duration2 < t <= 55 * duration2:
                Bi = anuga.Dirichlet_boundary([0.0116,0.0,0.0])
            if 55 * duration2 < t <= 60 * duration2:
                Bi = anuga.Dirichlet_boundary([0.8765,0.0,0.0])

            domain.set_boundary({'left': Bi,'right': Bo,'top': Br,'bottom': Br})
            #print t