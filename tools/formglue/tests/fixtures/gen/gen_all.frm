#-
* Generation of router test inputs.
* The Model QCD block below (Particle/Vertex declarations) is the example model of
* the FORM 5 reference manual, chapter "Diagram generation" (the gg->qqbar diagrams_
* demo): J.A.M. Vermaseren et al., form-dev/form, GPL-3.0,
* https://github.com/form-dev/form. It uses only FORM's documented Model/diagrams_
* syntax; no FORM source or manual text is included. The .out fixtures beside this
* file are the output of running FORM 5.0.1 on this script (gen_all.log), and
* tests/fixtures/enum_1loop.out is the same 30-diagram listing that FORM prints for
* the manual's demo; FORM is obtained upstream and run as a separate program.
* (Same QCD model as the shipped harness QCD model.)
Off statistics;
V p1,p2,q1,q2,k1,...,k12;
Set kint: k1,...,k12;

Model QCD;
  Particle qua,QUA,-2;
  Particle gho,GHO,-1;
  Particle glu,+3;
  Vertex QUA,qua,glu:g;
  Vertex GHO,gho,glu:g;
  Vertex glu,glu,glu:g;
  Vertex glu,glu,glu,glu:g^2;
EndModel;

* ---- (1) demo's 30-diagram on-shell 1-loop set, node form (byte-comparable
*      to the reference fixtures/enum_1loop.out) ----
L T2 = diagrams_(QCD,{glu,glu},{qua,QUA},{p1,p2,q1,q2},kint,1,
    `OnShell_'+`NoTadpole_'+`NoSnail_');
.sort
#$n2 = termsin_(T2);
#write "COUNT onshell_nodes = `$n2'"
Format nospaces;
#write <enum_onshell_nodes.out> "%E", T2
#close <enum_onshell_nodes.out>
Drop T2;
.sort

* ---- (2) same set WithEdges_ ----
L T2e = diagrams_(QCD,{glu,glu},{qua,QUA},{p1,p2,q1,q2},kint,1,
    `OnShell_'+`NoTadpole_'+`NoSnail_'+`WithEdges_');
.sort
#$n2e = termsin_(T2e);
#write "COUNT onshell_edges = `$n2e'"
Format nospaces;
#write <enum_onshell_edges.out> "%E", T2e
#close <enum_onshell_edges.out>
Drop T2e;
.sort

* ---- (3) census form: TopologiesOnly_ + WithEdges_ ----
L T2t = diagrams_(QCD,{glu,glu},{qua,QUA},{p1,p2,q1,q2},kint,1,
    `OnShell_'+`NoTadpole_'+`NoSnail_'+`TopologiesOnly_'+`WithEdges_');
.sort
#$n2t = termsin_(T2t);
#write "COUNT census_topo = `$n2t'"
Format nospaces;
#write <census_topo.out> "%E", T2t
#close <census_topo.out>
Drop T2t;
.sort

* ---- (4) QCD 3-loop vacuum (battery C count: 74), node form ----
L V3 = diagrams_(QCD,{},{},{p1,p2,q1,q2},kint,3,0);
.sort
#$nv = termsin_(V3);
#write "COUNT vac3_nodes = `$nv'"
Format nospaces;
#write <vac3_nodes.out> "%E", V3
#close <vac3_nodes.out>
Drop V3;
.sort

* ---- (5) same vacuum set WithEdges_ ----
L V3e = diagrams_(QCD,{},{},{p1,p2,q1,q2},kint,3,`WithEdges_');
.sort
#$nve = termsin_(V3e);
#write "COUNT vac3_edges = `$nve'"
Format nospaces;
#write <vac3_edges.out> "%E", V3e
#close <vac3_edges.out>
Drop V3e;
.sort
.end
