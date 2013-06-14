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

;;; type:
;;; h+   ->-|   '(+ ? ?)
;;; p+   -<-|   '(? + ?)
;;; h-   |->-   '(- ? ?)
;;; p-   |-<-   '(? - ?)
;;; kr   Kronecker delta
;;;
;;; (t_{ijk}^{abc} (i a 1) (j b 2) (k c 3))

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
               (unless (or (numberp (hole-of pair-op))
                           (equal last-pair this-pair))
                 (setq last-pair this-pair)
                 (make-pair-op hole-ind
                               (particle-of pair-op)
                               (index-of pair-op))))))
      (mapcar (lambda (ops)
                (attach-tag 't ops))
              (mapreplace #'commute (content-of amp))))))
(defun contract-particle-amp (particle-ind amp)
  (let ((last-pair '()))
    (flet ((commute (pair-op)
             (let ((this-pair (op-w/o-index pair-op)))
               (unless (or (numberp (particle-of pair-op))
                           (equal last-pair this-pair))
                 (setq last-pair this-pair)
                 (make-pair-op (hole-of pair-op)
                               particle-ind
                               (index-of pair-op))))))
      (mapcar (lambda (ops)
                (attach-tag 't ops))
              (mapreplace #'commute (content-of amp))))))
; return the sum of a list of ampprod
(defun amp-w/o-index (amp)
  (mapcar #'op-w/o-index (content-of amp)))
(defun contract-hole-ampprod (hole-ind ampprod)
  (let ((last-amp '()))
    (mapreplace* (lambda (amp)
                   (let ((this-amp (amp-w/o-index amp)))
                     (unless (equal last-amp this-amp)
                       (setq last-amp this-amp)
                       (contract-hole-amp hole-ind amp))))
                 ampprod)))
(defun contract-particle-ampprod (particle-ind ampprod)
  (let ((last-amp '()))
    (mapreplace* (lambda (amp)
                   (let ((this-amp (amp-w/o-index amp)))
                     (unless (equal last-amp this-amp)
                       (setq last-amp this-amp)
                       (contract-particle-amp particle-ind amp))))
                 ampprod)))
    
; return a list of amplitudes
(defun contract-ops-ampprods (pair-ops ampprod-lst)
  (flet ((link-h (pair-op ampprod)
           (if (eql (hole-of pair-op) '-)
               (contract-hole-ampprod (index-of pair-op) ampprod)
               (list ampprod)))
         (link-p (pair-op ampprod)
           (if (eql (particle-of pair-op) '-)
               (contract-particle-ampprod (index-of pair-op) ampprod)
               (list ampprod)))
         (contract-iter (func pair-ops ampprod-lst)
           (reduce (lambda (ts-lst pair-op)
                     (mapcan (lambda (ampprod) (funcall func pair-op ampprod))
                             ts-lst))
                   pair-ops :initial-value ampprod-lst)))
    (contract-iter #'link-p pair-ops
                   (contract-iter #'link-h pair-ops ampprod-lst))))
(defun contract-h2e-ampprod (h2e ampprod)
  (contract-ops-ampprods (content-of h2e) (list ampprod)))

; if the 2e-operator is symmetric, remove those duplicated amplitude products
(defun contract-symm-h2e-ampprod (filter h2e ampprod)
  (let ((ampprod-lst (contract-h2e-ampprod h2e ampprod)))
    (remove-if-not filter ampprod-lst)))
(defun ordering-h2e-symm-hole (ampprod)
  (flet ((get-hole-index (amp)
           (mapcan (lambda (pair-op)
                     (let ((hole-idx (hole-of pair-op)))
                       (and (numberp hole-idx) (list hole-idx))))
                   (content-of amp))))
    (format t "~a~%" (mapcan #'get-hole-index ampprod))
    (apply #'< (mapcan #'get-hole-index ampprod))))

(defun ordering-h2e-symm-particle (ampprod)
  (flet ((get-hole-index (amp)
           (mapcan (lambda (pair-op)
                     (let ((particle-idx (particle-of pair-op)))
                       (and (numberp particle-idx) (list particle-idx))))
                   (content-of amp))))
    (apply #'< (mapcan #'get-hole-index ampprod))))

(defun id-pair-op (pair-op &optional (oplist '(1 2 +)))
  (let ((shift (length oplist))
        (h-id (position (hole-of pair-op) oplist))
        (p-id (position (particle-of pair-op) oplist)))
    (+ (* shift h-id) p-id)))

; id of amplitude in terms of id-pair-op
(defun id-amp (amp &optional (oplist '(1 2 +)))
  (sort (mapcar (lambda (x) (id-pair-op x oplist))
                (content-of amp))
        #'<))

(defun list-lt (l1 l2)
  (cond ((null l1) t)
        ((null l2) nil)
        ((eq (car l1) (car l2))
         (list-lt (cdr l1) (cdr l2)))
        (t (< (car l1) (car l2)))))
(defun id-ampprod (ampprod &optional (oplist '(1 2 +)))
  (sort (mapcar (lambda (x) (id-amp x oplist)) ampprod)
        (lambda (x y)
          (let ((len-x (length x))
                (len-y (length y)))
            (if (eq len-x len-y)
                (list-lt x y)
                (< len-x len-y))))))


(defun remove-h2e-symm-dup (ampprod-lst oplist)
  (flet ((ampprod-eql (ts1 ts2)
           (let ((swaplist `(,(second oplist) ,(first oplist) ,@(cddr oplist))))
             (equal (id-ampprod ts1 oplist)
                    (id-ampprod ts2 swaplist)))))
    (remove-duplicates ampprod-lst :test #'ampprod-eql :from-end t)))



(defun connected-amp? (amp)
  (notevery (lambda (pair-op)
              (and (eql (hole-of pair-op) 'h+)
                   (eql (particle-of pair-op) 'p+)))
            (content-of amp)))
(defun connected-ampprod? (ampprod)
  (every #'connected-amp? ampprod))
(defun remove-unconnected (ampprod-lst)
  (remove-if-not #'connected-ampprod? ampprod-lst))


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
