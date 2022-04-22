(TeX-add-style-hook
 "Problem2"
 (lambda ()
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

