;(compile-file "")
;(require ')
;(load "ccd.cl")

(defmacro pass? (test value)
  `(format t "~a  ~:[FAIL~;pass~]~%" ',test (equal ,test ,value)))

(pass? (amplitude? '(t)) t)
;(pass? (amplitudes? '((t) (t))) t)
;(pass? (amplitudes? '((v) (t))) nil)

(pass? (mapreplace #'evenp '(1 2 3))
       '((1 t 3)))
(pass? (mapreplace #'evenp '(1 1 3))
       '())
(pass? (mapreplace* (lambda (x) (mapreplace #'evenp x)) '((1) (2 3)))
       '((1) (t 3)))
(pass? (mapreplace* (lambda (x) (mapreplace #'evenp x)) '((1) (2) (3)))
       '((1) (t) (3)))

(pass? (sort-pair-ops '((a a 3) (b b 1) (c c 4) (d d 2)))
       '((B B 1) (D D 2) (A A 3) (C C 4)))



; contract with a t3 amplitude
(pass? (contract-hole-amp 9 '(t (h+ p 1) (h+ p 2) (4 p 3)))
       '(t (9 p 1) (h+ p 2) (4 p 3)))
(pass? (contract-particle-amp 9 '(t (h p 1) (h 4 2) (h p 3)))
       '(t (h 9 1) (h 4 2) (h p 3)))
(pass? (contract-particle-amp 9 '(t (h 4 1)))
       nil)
(pass? (contract-hole-amps 9 '((t (h+ p 1)) (t (4 p 2) (h+ p 3))))
       '(((t (9 p 1)) (t (4 p 2) (h+ p 3)))
         ((t (h+ p 1)) (t (4 p 2) (9 p 3)))))
(pass? (contract-particle-amps 9 '((t (h+ p+ 1)) (t (h+ p+ 2))))
       '(((t (h+ 9 1)) (t (h+ p+ 2)))))

(pass? (contract-h2e-amps '(g (h- p- 8) (h- p- 9)) '((t (h+ p+ 1))))
       '())
(pass? (contract-h2e-amps '(g (h- p- 8) (h- p+ 9))
                          '((t (h+ p+ 1)) (t (h+ p+ 2))))
       '(((t (8 8 1)) (t (9 p+ 2))) ((t (8 p+ 1)) (t (9 8 2)))))
(pass? (contract-h2e-amps '(g (h+ p- 8) (h- p+ 9))
                          '((t (h+ p+ 1) (h+ p+ 2))))
       '(((T (9 8 1) (H+ P+ 2))) ((T (H+ 8 1) (9 P+ 2)))))
(pass? (contract-h2e-amps '(g (h- p+ 8) (h+ p+ 9))
                          '((t (h+ p+ 1) (h+ p+ 2) (h+ p+ 3))))
       '(((t (8 p+ 1) (h+ p+ 2) (h+ p+ 3)))))
(pass? (contract-h2e-amps '(g (h- p- 8) (h- p- 9))
                          '((t (h+ p+ 1) (h+ p+ 2) (h+ p+ 3))))
       '(((t (8 8 1) (9 9 2) (h+ p+ 3))) ((t (8 8 1) (9 p+ 2) (h+ 9 3)))
         ((t (8 9 1) (9 8 2) (h+ p+ 3))) ((t (8 p+ 1) (9 8 2) (h+ 9 3)))
         ((t (8 9 1) (h+ 8 2) (9 p+ 3))) ((t (8 p+ 1) (h+ 8 2) (9 9 3)))))
(pass? (contract-h2e-amps '(g (h- p- 8) (h- p- 9))
                          '((t (h+ p+ 1)) (t (h+ p+ 3))))
       '(((t (8 8 1)) (t (9 9 3))) ((t (8 9 1)) (t (9 8 3)))))
(pass? (contract-h2e-amps '(g (h- p- 8) (h- p+ 9))
                          '((t (h+ p+ 1) (h+ p+ 2)) (t (h+ p+ 3) (h+ p+ 4))))
       '(((t (8 8 1) (9 p+ 2)) (t (h+ p+ 3) (h+ p+ 4)))
         ((t (8 8 1) (h+ p+ 2)) (t (9 p+ 3) (h+ p+ 4)))
         ((t (8 p+ 1) (9 8 2)) (t (h+ p+ 3) (h+ p+ 4)))
         ((t (8 p+ 1) (h+ 8 2)) (t (9 p+ 3) (h+ p+ 4)))
         ((t (8 p+ 1) (9 p+ 2)) (t (h+ 8 3) (h+ p+ 4)))
         ((t (8 p+ 1) (h+ p+ 2)) (t (9 8 3) (h+ p+ 4)))
         ((t (8 p+ 1) (h+ p+ 2)) (t (h+ 8 3) (9 p+ 4)))))
(pass? (contract-h2e-amps '(g (h- p- 8) (h- p- 9))
                          '((t (h+ p+ 1) (h+ p+ 2)) (t (h+ p+ 3) (h+ p+ 4))))
       '(((t (8 8 1) (9 9 2)) (t (h+ p+ 3) (h+ p+ 4)))
         ((t (8 8 1) (9 p+ 2)) (t (h+ 9 3) (h+ p+ 4)))
         ((t (8 8 1) (h+ 9 2)) (t (9 p+ 3) (h+ p+ 4)))
         ((t (8 8 1) (h+ p+ 2)) (t (9 9 3) (h+ p+ 4)))
         ((t (8 8 1) (h+ p+ 2)) (t (9 p+ 3) (h+ 9 4)))
         ((t (8 9 1) (9 8 2)) (t (h+ p+ 3) (h+ p+ 4)))
         ((t (8 p+ 1) (9 8 2)) (t (h+ 9 3) (h+ p+ 4)))
         ((t (8 9 1) (h+ 8 2)) (t (9 p+ 3) (h+ p+ 4)))
         ((t (8 p+ 1) (h+ 8 2)) (t (9 9 3) (h+ p+ 4)))
         ((t (8 p+ 1) (h+ 8 2)) (t (9 p+ 3) (h+ 9 4)))
         ((t (8 9 1) (9 p+ 2)) (t (h+ 8 3) (h+ p+ 4)))
         ((t (8 p+ 1) (9 9 2)) (t (h+ 8 3) (h+ p+ 4)))
         ((t (8 p+ 1) (9 p+ 2)) (t (h+ 8 3) (h+ 9 4)))
         ((t (8 9 1) (h+ p+ 2)) (t (9 8 3) (h+ p+ 4)))
         ((t (8 p+ 1) (h+ 9 2)) (t (9 8 3) (h+ p+ 4)))
         ((t (8 p+ 1) (h+ p+ 2)) (t (9 8 3) (h+ 9 4)))
         ((t (8 9 1) (h+ p+ 2)) (t (h+ 8 3) (9 p+ 4)))
         ((t (8 p+ 1) (h+ 9 2)) (t (h+ 8 3) (9 p+ 4)))
         ((t (8 p+ 1) (h+ p+ 2)) (t (h+ 8 3) (9 9 4)))))

;(pass? (commutable? '(h+ 1) '(h- 2)) t)
;(pass? (norm-ordering? '((h+ 1) (h- 2) (p- 1) (p+ 2))) nil)
;
;(pass? (contract-a-pseudop '(p- 4)
;                           '((1 (h+ 3)) (p+ 1) (p+ 2)))
;       '(((1 (kr 4) (h+ 3)) (kr 4) (p+ 2))
;         ((-1 (p- 4) (h+ 3)) (p+ 1) (p+ 2))))
;(pass? (contract-a-pseudop '(h+ 4)
;                           '((1 (h+ 3)) (p+ 1) (p+ 2)))
;       '(((1 (h+ 4) (h+ 3)) (p+ 1) (p+ 2))))
;(pass? (contract-a-pseudop '(h- 4)
;                           '((1 (h+ 3)) (h+ 1) (h+ 2)))
;       '(((1 (h- 4) (h+ 3)) (h+ 1) (h+ 2))
;         ((1 (kr 4) (h+ 3)) (kr 4) (h+ 2))
;         ((-1 (kr 4) (h+ 3)) (h+ 1) (kr 4))))
;(pass? (contract-a-pseudop '(p- 4)
;                           '((1) (kr 1) (p+ 2)))
;       '(((1 (kr 4)) (kr 1) (kr 4))
;         ((-1 (p- 4)) (kr 1) (p+ 2))))
;
;(pass? (contract-as-pseudops '((p- 4) (h- 5))
;                             '(((1 (h+ 3)) (p+ 1) (h+ 2))))
;       '(((1 (kr 4) (h- 5) (h+ 3)) (kr 4) (h+ 2))
;         ((1 (p- 4) (h- 5) (h+ 3)) (p+ 1) (h+ 2))
;         ((-1 (kr 4) (kr 5) (h+ 3)) (kr 4) (kr 5))
;         ((1 (p- 4) (kr 5) (h+ 3)) (p+ 1) (kr 5))))
;(pass? (contract-as-pseudops '((p- 4) (h- 5))
;                             '(((1) (kr 1) (h+ 2))))
;       '(((1 (p- 4) (h- 5)) (kr 1) (h+ 2))
;         ((1 (p- 4) (kr 5)) (kr 1) (kr 5))))
;
;(pass? (norm-order-forward-op-op '(v ((h- 3)) ((p- 4)))
;                                 '(t ((p+ 1)) ((h+ 2))))
;       '((-1 (v ((h- 3)) ((kr 4))) (t ((kr 4)) ((h+ 2))))
;         (1 (v ((h- 3)) ((p- 4))) (t ((p+ 1)) ((h+ 2))))
;         (1 (v ((kr 3)) ((kr 4))) (t ((kr 4)) ((kr 3))))
;         (-1 (v ((kr 3)) ((p- 4))) (t ((p+ 1)) ((kr 3))))))
;
;(pass? (norm-order-forward-op-op '(v ((h- 3)) ((kr 4)))
;                                 '(t ((kr 1)) ((h+ 2))))
;       '((-1 (v ((h- 3)) ((kr 4))) (t ((kr 1)) ((h+ 2))))
;         (1 (v ((kr 3)) ((kr 4))) (t ((kr 1)) ((kr 3))))))
;
;(pass? (norm-order-backward-op-op '(v ((h- 3)) ((p- 4)))
;                                  '(t ((p+ 1)) ((h+ 2))))
;       '((-1 (v ((kr 2)) ((p- 4))) (t ((p+ 1)) ((kr 2))))
;         (1 (v ((kr 2)) ((kr 1))) (t ((kr 1)) ((kr 2))))
;         (1 (v ((h- 3)) ((p- 4))) (t ((p+ 1)) ((h+ 2))))
;         (-1 (v ((h- 3)) ((kr 1))) (t ((kr 1)) ((h+ 2))))))
;
;(pass? (norm-order-backward-op-op '(v ((h- 3)) ((kr 4)))
;                                  '(t ((kr 1)) ((h+ 2))))
;       '((1 (v ((kr 2)) ((kr 4))) (t ((kr 1)) ((kr 2))))
;         (-1 (v ((h- 3)) ((kr 4))) (t ((kr 1)) ((h+ 2))))))
;
;(pass? (remove-uncontract-cell
;         (norm-order-forward-op-op '(v ((h- 1) (h- 2)) ((p- 4) (p- 3)))
;                                   '(t ((p+ 5) (p+ 6)) ((h+ 8) (h+ 7)))))
;       '((1 (v ((kr 1) (kr 2)) ((kr 4) (kr 3)))
;            (t ((kr 3) (kr 4)) ((kr 2) (kr 1))))
;         (-1 (v ((kr 1) (kr 2)) ((kr 4) (kr 3)))
;             (t ((kr 3) (kr 4)) ((kr 1) (kr 2))))))
;
;(pass? (remove-uncontract-cell
;         (norm-order-backward-op-op '(v ((h- 1) (h- 2)) ((p- 4) (p- 3)))
;                                    '(t ((p+ 5) (p+ 6)) ((h+ 8) (h+ 7)))))
;       '((1 (v ((kr 8) (kr 7)) ((kr 5) (kr 6)))
;            (t ((kr 5) (kr 6)) ((kr 8) (kr 7))))
;         (-1 (v ((kr 8) (kr 7)) ((kr 6) (kr 5)))
;             (t ((kr 5) (kr 6)) ((kr 8) (kr 7))))))
;
;(pass? (norm-order-op-cell '(v (p- 3 1))
;                           '((t (p- 1 1)) (t (p+ 2 -1))))
;       '(((v (kr 2 1)) (t (p- 1 -1)) (t (kr 2 -1)))
;         ((v (p- 3 1)) (t (p- 1 -1)) (t (p+ 2 1)))))
;
;(pass? (norm-order-op-cell '(v (p- 4 1) (h- 5 1))
;                           '((t (p- 1 1)) (t (p+ 2 -1) (h+ 3 -1))))
;       '(((v (kr 2 1) (kr 3 1)) (t (p- 1 1)) (t (kr 2 1) (kr 3 -1)))
;         ((v (p- 4 1) (kr 3 1)) (t (p- 1 1)) (t (p+ 2 -1) (kr 3 -1)))
;         ((v (kr 2 1) (h- 5 1)) (t (p- 1 1)) (t (kr 2 1) (h+ 3 1)))
;         ((v (p- 4 1) (h- 5 1)) (t (p- 1 1)) (t (p+ 2 -1) (h+ 3 -1)))))
;
;(pass? (norm-order-op-cells '(t (h- 5 1))
;                            '(((t (p+ 1 1))) ((t (h+ 2 -1)) (v (kr 3 1)))))
;       '(((t (h- 5 1)) (t (p+ 1 -1)))
;         ((t (kr 2 1)) (t (kr 2 -1)) (v (kr 3 1)))
;         ((t (h- 5 1)) (t (h+ 2 1)) (v (kr 3 1)))))
;
;(pass? (norm-order-cell-cell '((v (p- 4 1)) (t (h- 5 1)))
;                             '((t (p+ 1 1)) (t (h+ 3 -1))))
;       '(((v (kr 1 1)) (t (kr 3 1)) (t (kr 1 -1)) (t (kr 3 -1)))
;         ((v (p- 4 1)) (t (kr 3 1)) (t (p+ 1 1)) (t (kr 3 -1)))
;         ((v (kr 1 1)) (t (h- 5 -1)) (t (kr 1 -1)) (t (h+ 3 1)))
;         ((v (p- 4 1)) (t (h- 5 -1)) (t (p+ 1 1)) (t (h+ 3 -1)))))
;
;(pass? (norm-order-cell-cells '((v (h- 4 1)) (t (h- 5 1)))
;                              '(((t (h+ 1 1)) (v (h+ 2 -1)))
;                                ((t (kr 3 1)) (v (h+ 2 1)))))
;       '(((v (kr 2 1)) (t (kr 1 1))  (t (kr 1 1))  (v (kr 2 -1)))
;         ((v (h- 4 1)) (t (kr 1 1))  (t (kr 1 1))  (v (h+ 2 1)))
;         ((v (kr 1 1)) (t (kr 2 1))  (t (kr 1 -1)) (v (kr 2 -1)))
;         ((v (h- 4 1)) (t (kr 2 1))  (t (h+ 1 1))  (v (kr 2 -1)))
;         ((v (kr 1 1)) (t (h- 5 -1)) (t (kr 1 -1)) (v (h+ 2 1)))
;         ((v (kr 2 1)) (t (h- 5 -1)) (t (h+ 1 1))  (v (kr 2 1)))
;         ((v (h- 4 1)) (t (h- 5 -1)) (t (h+ 1 1))  (v (h+ 2 -1)))
;         ((v (h- 4 1)) (t (kr 2 1))  (t (kr 3 1))  (v (kr 2 1)))
;         ((v (kr 2 1)) (t (h- 5 -1)) (t (kr 3 1))  (v (kr 2 -1)))
;         ((v (h- 4 1)) (t (h- 5 -1)) (t (kr 3 1))  (v (h+ 2 1)))))
;
;(pass? (remove-uncontract-cell
;         (norm-order-cell-cell '((v (h- 1 1) (p- 2 1) (h- 3 1) (p- 4 1)))
;                      '((t (h+ 5 1) (p+ 6 1) (h+ 7 1) (p+ 8 1)))))
;       '(((v (kr 5 1) (kr 8 1) (kr 7 1) (kr 6 1))
;          (t (kr 5 -1) (kr 6 1) (kr 7 1) (kr 8 1)))
;         ((v (kr 5 1) (kr 6 1) (kr 7 1) (kr 8 1))
;          (t (kr 5 -1) (kr 6 1) (kr 7 -1) (kr 8 1)))))
;
;(pass? (linked? '((t (kr 4 -1) (kr 3 1)) (v (kr 4 1) (kr 3 1))))
;       t)
;(pass? (linked? '((t (kr 4 -1) (kr 3 1)) (v (p+ 4 1) (p- 3 1))))
;       nil)
;
;(pass? (exist-uncontracted? 
;         '((t (kr 4 -1) (kr 3 1)) (v (kr 4 1) (kr 3 1))))
;       nil)
;(pass? (exist-uncontracted? 
;        '((t (p+ 1 1) (kr 3 1)) (v (p- 4 1) (kr 3 1))))
;       t)