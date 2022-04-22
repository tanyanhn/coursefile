(TeX-add-style-hook
 "Problem3"
 (lambda ()
   (LaTeX-add-labels
    "eq:characteristics")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

