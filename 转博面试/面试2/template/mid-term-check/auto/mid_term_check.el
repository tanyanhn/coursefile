(TeX-add-style-hook
 "mid_term_check"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "CJKutf8"
    "amsfonts"
    "amsmath"
    "amssymb"
    "amsthm"
    "enumerate"
    "graphicx"
    "layout"
    "mathrsfs"
    "fancyhdr"
    "subfigure"
    "tcolorbox"
    "tikz-cd"
    "color"
    "pifont"
    "verbatim"
    "mathtools"
    "float"
    "bm")
   (TeX-add-symbols
    "dif")
   (LaTeX-add-labels
    "eq:INSE"
    "eq:GePUP")
   (LaTeX-add-amsthm-newtheorems
    "thm"))
 :latex)

