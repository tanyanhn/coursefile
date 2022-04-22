(TeX-add-style-hook
 "Problem1"
 (lambda ()
   (LaTeX-add-labels
    "eq:1"
    "eq:2"
    "eq:3"
    "eq:4"
    "eq:5"
    "eq:6")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

