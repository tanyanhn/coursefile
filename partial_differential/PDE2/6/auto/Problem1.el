(TeX-add-style-hook
 "Problem1"
 (lambda ()
   (LaTeX-add-labels
    "eq:1")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

