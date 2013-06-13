;;;;
;;;; File: ccd.cl
;;;; Author
;;;; Description:

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; when compile with sbcl.cl -c, uncomment these lines
;(eval-when (:compile-toplevel :load-toplevel :execute)
;  (load "utility.fasl"))
(load "utility.cl")
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun attach-tag (tag x)
  (cons tag x))
(defmacro tag-of (x)
  `(car ,x))
(defmacro content-of (x)
  `(cdr ,x))
;;; commutator : p q^\dagger + q^\dagger p = \delta_{pq}

;;; A node: factor * type_index = (type index factor)
;;; e.g. -1a_i^\dagger = (h- i -1)
;;; type:
;;; h+   ->-|
;;; p+   -<-|
;;; h-   |->-
;;; p-   |-<-
;;; kr   Kronecker delta
;;;
;;; (t_{ijkl}^{abcd} ijkl abcd)
;;; As : list of nodes
;;; '((h+ 1) (p+ 1) (h+ 2) (p+ 2) (h+ 3) (p+ 3))
(defun make-pair-op (hole particle ind)
  (list hole particle ind))
(defun hole-of (a)
  (car a))
(defun particle-of (a)
  (cadr a))
(defun index-of (a)
  (caddr a))
(defun op-w/o-index (a)
  (list (hole-of a) (particle-of a)))

(defun amplitude? (a)
  (eql (tag-of a) 't))
;(defun amplitudes? (a)
;  (every #'amplitude? a))
(defun h2e? (a)
  (eql (tag-of a) 'g))

;;; almost like maplist
;;; if (func lst) return single result
(defun mapreplace (func lst)
  (if (null lst)
      '()
      (let* ((item1 (car lst))
             (item-rest (cdr lst))
             (res1 (funcall func item1))
             (res-rest (mapcar (lambda (a) (cons item1 a))
                               (mapreplace func item-rest))))
        (if (null res1)
            res-rest
            (cons (cons res1 item-rest)
                  res-rest)))))
;;; in case (func lst) return a list of results, flatten the return results
(defun mapreplace* (func lst)
  (if (null lst)
      '()
      (let* ((item1 (car lst))
             (item-rest (cdr lst))
             (res1 (funcall func item1))
             (res-rest (mapcar (lambda (a) (cons item1 a))
                               (mapreplace* func item-rest))))
        (if (null res1)
            res-rest
            (append (mapcar (lambda (a) (cons a item-rest)) res1)
                    res-rest)))))

(defun sort-pair-ops (ops)
  (sort ops (lambda (op1 op2)
              (< (index-of op1) (index-of op2)))))
;;; all possible contraction of one hole and an amplitude
;; return the sum of a list of amp
;; save last status, to avoid double counting of symmetric amplitudes
(defun contract-hole-amp (hole-ind amp)
  (let ((last-pair '()))
    (flet ((commute (pair-op)
             (let ((this-pair (op-w/o-index pair-op)))
               (if (or (numberp (hole-of pair-op))
                       (equal last-pair this-pair))
                   nil
                   (progn (setq last-pair this-pair)
                          (make-pair-op hole-ind
                                        (particle-of pair-op)
                                        (index-of pair-op)))))))
      (mapcar (lambda (ops)
                (attach-tag 't ops))
              (mapreplace #'commute (content-of amp))))))
(defun contract-particle-amp (particle-ind amp)
  (let ((last-pair '()))
    (flet ((commute (pair-op)
             (let ((this-pair (op-w/o-index pair-op)))
               (if (or (numberp (particle-of pair-op))
                       (equal last-pair this-pair))
                   nil
                   (progn (setq last-pair this-pair)
                          (make-pair-op (hole-of pair-op)
                                        particle-ind
                                        (index-of pair-op)))))))
      (mapcar (lambda (ops)
                (attach-tag 't ops))
              (mapreplace #'commute (content-of amp))))))
; return the sum of a list of ampprod
(defun amp-w/o-index (amp)
  (mapcar #'op-w/o-index (cdr amp)))
(defun contract-hole-ampprod (hole-ind ampprod)
  (let ((last-amp '()))
    (mapreplace* (lambda (amp)
                   (let ((this-amp (amp-w/o-index amp)))
                     (if (equal last-amp this-amp)
                         nil
                         (progn (setq last-amp this-amp)
                                (contract-hole-amp hole-ind amp)))))
                 ampprod)))
(defun contract-particle-ampprod (particle-ind ampprod)
  (let ((last-amp '()))
    (mapreplace* (lambda (amp)
                   (let ((this-amp (amp-w/o-index amp)))
                     (if (equal last-amp this-amp)
                         nil
                         (progn (setq last-amp this-amp)
                                (contract-particle-amp particle-ind amp)))))
                 ampprod)))
    
; return a list of amplitudes
(defun contract-h2e-ampprod (h2e ampprod)
  (let* ((pair-e1 (first (content-of h2e)))
         (pair-e2 (second (content-of h2e)))
         (h2e-ops `((,(hole-of pair-e1) ,(index-of pair-e1))
                    (,(particle-of pair-e1) ,(index-of pair-e1))
                    (,(hole-of pair-e2) ,(index-of pair-e2))
                    (,(particle-of pair-e2) ,(index-of pair-e2)))))
    (flet ((contract-iter (h2e-op ampprod)
             (let ((op-symb (first h2e-op))
                   (op-index (second h2e-op)))
               (cond ((eql op-symb 'h-)
                      (contract-hole-ampprod op-index ampprod))
                     ((eql op-symb 'p-)
                      (contract-particle-ampprod op-index ampprod))
                     (t (list ampprod))))))
      (reduce (lambda (ampprod-lst op) ; apply op to every product of amp
                (mapcan (lambda (ampprod) (contract-iter op ampprod))
                        ampprod-lst))
              h2e-ops :initial-value (list ampprod)))))

(defun connected-amp? (amp)
  (notevery (lambda (pair-op)
              (and (eql (hole-of pair-op) 'h+)
                   (eql (particle-of pair-op) 'p+)))
            (cdr amp)))
(defun connected-ampprod? (ampprod)
  (every #'connected-amp? ampprod))
(defun remove-unconnected (ampprod-lst)
  (remove-if-not #'connected-ampprod? ampprod-lst))

;(defun remove-h2e-symmetric-double-counting)

;;;;;;;;;;;;;;;;;;;
;(defun linked? (cell)
;  (labels ((linked-ops? (ops)
;             (cond ((null ops) t)
;                   ((member 'kr (flatten (car ops)))
;                    (linked-ops? (cdr ops)))
;                   (t nil))))
;    (linked-ops? (cdr cell))))
;(defun remove-unlink-cell (cells)
;  (remove-if-not #'linked? cells))
;
;(defun exist-uncontracted? (cell)
;  (cond ((null cell) t)
;        ((intersection (flatten cell) '(h+ h- p+ p-)) t)
;        (t nil)))
;(defun remove-uncontract-cell (cells)
;  (remove-if #'exist-uncontracted? cells))
;
;
;;;; energy expression
;;;; H_0 + [H_0, T_2] + [[H_0, T_2], T_2]
;;;;     + [[[H_0, T_2], T_2], T_2] + [[[[H_0, T_2], T_2], T_2], T_2]
;;;;non-symmetric v
;(defvar vsym '((1 (v (p+ 1) (h+ 2) (p+ 3) (h+ 4)))
;               (1 (v (h- 1) (h+ 2) (h- 3) (h+ 4)))
;               (1 (v (p+ 1) (p- 2) (p+ 3) (p- 4)))
;               (1 (v (h- 1) (p- 2) (h- 3) (p- 4)))))
;(defvar vasym '((1 (v (p+ 1) (h+ 2) (h- 3) (h+ 4)))
;                (1 (v (p+ 1) (h+ 2) (p+ 3) (p- 4)))
;                (1 (v (p+ 1) (h+ 2) (h- 3) (p- 4)))
;                (1 (v (h- 1) (h+ 2) (p+ 3) (p- 4)))
;                (1 (v (h- 1) (h+ 2) (h- 3) (p- 4)))
;                (1 (v (p+ 1) (p- 2) (h- 3) (p- 4)))))
;(defparameter h0 '((1 (v (p+ 1) (h+ 2) (p+ 3) (h+ 4)))
;             (1 (v (p+ 1) (h+ 2) (h- 3) (h+ 4)))
;             (1 (v (p+ 1) (h+ 2) (p+ 3) (p- 4)))
;             (1 (v (p+ 1) (h+ 2) (h- 3) (p- 4)))
;             (1 (v (h- 1) (h+ 2) (p+ 3) (h+ 4)))
;             (1 (v (h- 1) (h+ 2) (h- 3) (h+ 4)))
;             (1 (v (h- 1) (h+ 2) (p+ 3) (p- 4)))
;             (1 (v (h- 1) (h+ 2) (h- 3) (p- 4)))
;             (1 (v (p+ 1) (p- 2) (p+ 3) (h+ 4)))
;             (1 (v (p+ 1) (p- 2) (h- 3) (h+ 4)))
;             (1 (v (p+ 1) (p- 2) (p+ 3) (p- 4)))
;             (1 (v (p+ 1) (p- 2) (h- 3) (p- 4)))
;             (1 (v (h- 1) (p- 2) (p+ 3) (h+ 4)))
;             (1 (v (h- 1) (p- 2) (h- 3) (h+ 4)))
;             (1 (v (h- 1) (p- 2) (p+ 3) (p- 4)))
;             (1 (v (h- 1) (p- 2) (h- 3) (p- 4)))))
;(defparameter t1 '((1 (t (p+ 5) (h+ 6)))))
;(defparameter t2 '((1 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8)))))
;(defparameter t2s '((1 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8)))
;              (1 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8))
;                 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8)))
;              (1 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8))
;                 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8))
;                 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8)))
;              (1 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8))
;                 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8))
;                 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8))
;                 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8)))))
;(defparameter t2s1 '((1 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8)))
;              (1 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8))
;                 (t (p+ 5) (h+ 6) (p+ 7) (h+ 8)))))
;;TODO; make-index
;(defun energy-expression (ts)
;  (remove-uncontract-cell (norm-order-cells-cells h0 ts)))
;(energy-expression t2)
;(energy-expression t2s1)
;(energy-expression t2s)
;
;;;; equations for amplitudes
;(defun amplitude-equation (ts)
;  (remove-uncontract-cell
;    (norm-order-cells-cells '((1 (t (h- 5) (p- 6) (h- 7) (p- 8))))
;                            (remove-unlink-cell
;                              (norm-order-cells-cells h0 ts)))))
;(amplitude-equation t2)
;(amplitude-equation t2s)
;(remove-unlink-cell (norm-order-cells-cells h0 t2s1))
;(amplitude-equation t2s1)
;
;
;;;; vim: ft=lisp
