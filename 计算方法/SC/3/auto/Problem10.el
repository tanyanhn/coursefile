(TeX-add-style-hook
 "Problem10"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (LaTeX-add-labels
    "eq:Bsplines"
    "eq:1"
    "prop:1"
    "prop:4"
    "eq:2"
    "eq:3"
    "eq:4"
    "eq:5")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

