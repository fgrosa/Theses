% selfenergyquark.mp -- do not edit, generated automatically by MasterThesis_Grosa.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 90*1.0pt#, 70*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i);
vright (__o);
vconnect ("fermion", __i, __v1);
vconnect ("gluon,right,tension=0", __v1, __v2);
vconnect ("fermion", __v2, __o);
vconnect ("fermion", __v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
