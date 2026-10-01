#-
* 1-loop diagram enumeration for g g -> q qbar with the FORM 5
* integrated generator (grcc), QCD model with ghosts, one quark flavor.
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

* ---- tree (control: must be 3) ----
L T0 = diagrams_(QCD,{glu,glu},{qua,QUA},{p1,p2,q1,q2},kint,0,0);
.sort
#$n0 = termsin_(T0);
#write "COUNT gg_qqbar_tree = `$n0'"
Drop T0;
.sort

* ---- 1 loop, all connected ----
L T1 = diagrams_(QCD,{glu,glu},{qua,QUA},{p1,p2,q1,q2},kint,1,0);
.sort
#$n1 = termsin_(T1);
#write "COUNT gg_qqbar_1loop_all = `$n1'"
Drop T1;
.sort

* ---- 1 loop, amputated on-shell set (no tadpoles, no snails,
*      no external self-energies): the virtual-correction set ----
L T2 = diagrams_(QCD,{glu,glu},{qua,QUA},{p1,p2,q1,q2},kint,1,
    `OnShell_'+`NoTadpole_'+`NoSnail_');
.sort
#$n2 = termsin_(T2);
#write "COUNT gg_qqbar_1loop_onshell = `$n2'"
Format nospaces;
#write <enum_1loop.out> "%E", T2
#close <enum_1loop.out>
Drop T2;
.sort

* ---- 1 loop, 1PI only (for reference) ----
L T3 = diagrams_(QCD,{glu,glu},{qua,QUA},{p1,p2,q1,q2},kint,1,`OnePI_');
.sort
#$n3 = termsin_(T3);
#write "COUNT gg_qqbar_1loop_1PI = `$n3'"
Drop T3;
.sort
.end
