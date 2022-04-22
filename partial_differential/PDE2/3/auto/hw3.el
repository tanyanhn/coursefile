(TeX-add-style-hook
 "hw3"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "10pt" "a4paper")))
   (TeX-run-style-hooks
    "latex2e"
    "num/Problem1"
    "num/Problem2"
    "num/Problem3"
    "num/Problem4"
    "num/Problem5"
    "article"
    "art10"
    "enumerate"
    "geometry"
    "CJKutf8"
    "amsfonts"
    "amsmath"
    "amssymb"
    "amsthm"
    "graphicx"
    "layout"
    "multicol"
    "mathrsfs"
    "fancyhdr"
    "subfigure"
    "tcolorbox"
    "tikz-cd"
    "mathtools"
    "float"
    "bm"
    "booktabs")
   (TeX-add-symbols
    "solname"
    "tbbint"
    "dbbint"
    "bbint"
    "LP"
    "RP"
    "LI"
    "RI"
    "dif")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-amsthm-newtheorems
    "pro"
    "defn"
    "lem"))
 :latex)

