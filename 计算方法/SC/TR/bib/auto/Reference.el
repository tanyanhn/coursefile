(TeX-add-style-hook
 "Reference"
 (lambda ()
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :bibtex)

