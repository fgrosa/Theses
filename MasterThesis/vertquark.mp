% vertquark.mp -- do not edit, generated automatically by MasterThesis_Grosa.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 90*1.0pt#, 60*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i);
vright (__o1, __o2);
vconnect ("gluon", __i, __v1);
vconnect ("fermion", __o1, __v2, __v1, __v3, __o2);
vconnect ("gluon,tension=0", __v2, __v3);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
