(TeX-add-style-hook
 "COVID"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amssymb"
    "amsmath"
    "CJKutf8"
    "color"
    "enumerate"
    "fancyhdr"
    "geometry"
    "graphicx"
    "indentfirst"
    "latexsym"
    "mathrsfs"
    "subfigure"
    "textcomp"
    "amsthm"
    "algorithm"
    "algorithmic"
    "flafter"
    "booktabs"
    "longtable"
    "pxfonts"
    "cite")
   (TeX-add-symbols
    '("ShowRevised" 1)
    '("DRLN" 1)
    '("curve" 1)
    '("avg" 1)
    "cProj"
    "cProjLH"
    "dif"
    "dt"
    "Dim"
    "ebold"
    "ibold"
    "jbold"
    "nPoly"
    "Div"
    "Grad"
    "Proj"
    "Lapl"
    "MARS"
    "Zero"
    "One"
    "Int"
    "mathbb")
   (LaTeX-add-labels
    "fig:numbers"
    "fig:spline"
    "fig:fitting"
    "fig:data")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-amsthm-newtheorems
    "theorem"
    "proposition"
    "lemma"
    "corollary"
    "definition"
    "question"))
 :latex)

