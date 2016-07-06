% fmfdefault.mp -- do not edit, generated automatically by tesi_Grosa.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 100*1.0pt#, 75*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1);
vright (__o1, __o2);
vconnect ("plain", __i1, __v1);
vconnect ("fermion", __o1, __v1, __o2);
vlabel ("$\pi ^- (k)$", __i1);
vlabel ("$\overline {\nu }_e (p)$", __o1);
vlabel ("$e^- (q)$", __o2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
