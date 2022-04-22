(TeX-add-style-hook
 "hw1"
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
    "num/Problem6"
    "num/Problem7"
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
    "booktabs"
    "color"
    "listings")
   (TeX-add-symbols
    '("pdfFrac" 2)
    '("difFrac" 2)
    "LP"
    "RP"
    "dif"
    "solname"
    "tbbint"
    "dbbint"
    "bbint")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0))
   (LaTeX-add-amsthm-newtheorems
    "pro"
    "defn")
   (LaTeX-add-xcolor-definecolors
    "mygreen"
    "mygray"
    "mymauve"))
 :latex)

