% vaccumgluon.mp -- do not edit, generated automatically by MasterThesis_Grosa_diff.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 90*1.0pt#, 60*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i);
vright (__o);
vconnect ("gluon", __i, __v1);
vconnect ("gluon", __v2, __o);
vconnect ("gluon,right,tension=0.6", __v1, __v2, __v1);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
