% selfenergygluon.mp -- do not edit, generated automatically by MasterThesis_Grosa_diff.tex
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
vconnect ("gluon,tension=0.5", __i, __v, __o);
vconnect ("gluon,left=90", __v, __v);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
