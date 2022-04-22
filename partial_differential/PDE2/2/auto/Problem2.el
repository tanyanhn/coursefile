(TeX-add-style-hook
 "Problem2"
 (lambda ()
   (LaTeX-add-labels
    "eq:4"
    "eq:5"
    "eq:6")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

