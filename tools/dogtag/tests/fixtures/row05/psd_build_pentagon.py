#!/usr/bin/env python3
# pySecDec build of the 1-loop massless scalar pentagon (build step 1:
# generate + compile the integral library ONCE; the sampler then just calls it at
# many phase-space points).
import os, sys
from pySecDec.loop_integral import LoopIntegralFromGraph, loop_package
from pySecDec.code_writer import sum_package  # noqa

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

li = LoopIntegralFromGraph(
    internal_lines=[['0',[1,2]],['0',[2,3]],['0',[3,4]],['0',[4,5]],['0',[5,1]]],
    external_lines=[['p1',1],['p2',2],['p3',3],['p4',4],['p5',5]],
    replacement_rules=[
        ('p1*p1',0),('p2*p2',0),('p3*p3',0),('p4*p4',0),('p5*p5',0),
        ('p1*p2','v1/2'),('p2*p3','v2/2'),('p3*p4','v3/2'),('p4*p5','v4/2'),('p5*p1','v5/2'),
        ('p1*p3','(v4-v1-v2)/2'),('p2*p4','(v5-v2-v3)/2'),('p3*p5','(v1-v3-v4)/2'),
        ('p1*p4','(v2-v4-v5)/2'),('p2*p5','(v3-v5-v1)/2'),
    ])

loop_package(
    name='pentagon1L',
    loop_integral=li,
    real_parameters=['v1','v2','v3','v4','v5'],
    requested_orders=[2],          # through eps^2 (covers w<=2 finite content)
    form_optimization_level=2,
    decomposition_method='geometric',
)
print("pySecDec pentagon library generated.")
