(TeX-add-style-hook
 "Problem2"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (LaTeX-add-labels
    "eq:1"
    "eq:2"
    "eq:3")
   (LaTeX-add-environments
    '("sol" LaTeX-env-args ["argument"] 0)))
 :latex)

