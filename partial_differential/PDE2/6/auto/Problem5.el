(TeX-add-style-hook
 "Problem5"
 (lambda ()
   (LaTeX-add-labels
    "eq:8")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

