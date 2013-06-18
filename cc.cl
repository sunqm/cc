;;; -*- Mode: Common-Lisp; -*-
;;;; Description: based on Goldstone diagram

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
;;; h+   ->-|
;;; p+   -<-|
;;; h-   |->-
;;; p-   |-<-
;;; kr   Kronecker delta
;;;
;;; (t_{ijk}^{abc} (i j k) (a b c))

(defun op-copy (n op)
  (loop repeat n
       collect op))
(defun make-amp (holes particles)
  (list t holes particles))
(defun holes-of (a)
  (cadr a))
(defun particles-of (a)
  (caddr a))
; n-excitation amplitude without linking
(defun make-new-amp (n)
  (make-amp (op-copy n '+) (op-copy n '+)))
; t1 = (make-amp '(+) '(+))
; t2 = (make-amp (op-copy 2 '+) (op-copy 2 '+))

;;; contraction of one hole and an amplitude
(defun contract-hole-amp (amp)
  (let ((h-amp (replace-once* (lambda (op)
                                (unless (eql op 'k)
                                  'k))
                              (holes-of amp))))
    (if h-amp
        (make-amp h-amp (particles-of amp)))))
(defun contract-particle-amp (amp)
  (let ((p-amp (replace-once* (lambda (op)
                                (unless (eql op 'k)
                                  'k))
                              (particles-of amp))))
    (if p-amp
        (make-amp (holes-of amp) p-amp))))

;(defun amp-w/o-index (amp)
;  (mapcar #'op-w/o-index (content-of amp)))
(defun contract-hole-ampprod (ampprod)
  (let ((last-amp '()))
    (mapreplace (lambda (amp)
                  (unless (equal last-amp amp)
                    (setq last-amp amp)
                    (contract-hole-amp amp)))
                ampprod)))
(defun contract-particle-ampprod (ampprod)
  (let ((last-amp '()))
    (mapreplace (lambda (amp)
                  (unless (equal last-amp amp)
                    (setq last-amp amp)
                    (contract-particle-amp amp)))
                ampprod)))

(defun reduce-h-p (hole-func particle-func h-ops init-value)
  )
(defun transform-ops (ops)
  (let ((hs (mapcar (lambda (symb)
                      (cond ((eql symb '-) 'h-)
                            ((eql symb '+) 'h+)
                            (t symb)))
                    (holes-of ops)))
        (ps (mapcar (lambda (symb)
                      (cond ((eql symb '-) 'p-)
                            ((eql symb '+) 'p+)
                            (t symb)))
                    (particles-of ops))))
    (append hs ps)))
; return a list of amplitudes
(defun contract-ops-ampprods (ops ampprod-lst)
  (flet ((contract-iter (op ampprod)
           (cond ((eql op 'h-)
                  (contract-hole-ampprod ampprod))
                 ((eql op 'p-)
                  (contract-particle-ampprod ampprod))
                 (t (list ampprod)))))
    (reduce (lambda (ts-lst op)
              (mapcan (lambda (ampprod) (contract-iter op ampprod))
                      ts-lst))
            (transform-ops ops) :initial-value ampprod-lst)))

(defun contract-ops-ampprods-uniq (ops ampprod-lst)
  (remove-duplicates (contract-ops-ampprods ops ampprod-lst)
                     :test #'equal))

(defun connected-amp? (amp)
  (flet ((conn? (ops)
           (notevery (lambda (op) (eql op '+)) ops)))
    (or (conn? (holes-of amp))
        (conn? (particles-of amp)))))
(defun connected-ampprod? (ampprod)
  (every #'connected-amp? ampprod))

(defun contract-h2e-ampprod (h2e ampprod)
  (let ((diagrams (contract-ops-ampprods-uniq h2e (list ampprod))))
    (mapcar (lambda (amps) (cons h2e amps))
            (remove-if-not #'connected-ampprod? diagrams))))

(defparameter *h2e-ops*
  '((h (- +) (+ +))
    (h (+ +) (- +))
    (h (- -) (- +))
    (h (- +) (- -))
    (h (- +) (- +))
    (h (- -) (+ +))
    (h (+ +) (- -))
    (h (- -) (- -))))

(defun gen-amps-list (n-lst)
  (if (null n-lst)
      '()
      (let* ((tnew (make-new-amp (car n-lst)))
             (t1 (list tnew))
             (t2 (cons tnew t1))
             (t3 (cons tnew t2))
             (t4 (cons tnew t3))
             (amps-lst (gen-amps-list (cdr n-lst)))
             (new-lst (mapcan (lambda (tx)
                                (mapcar (lambda (amps) (append tx amps))
                                        amps-lst))
                              (list t1 t2 t3)))
             (mod-lst (remove-if (lambda (ts) (> (length ts) 4))
                                 new-lst)))
        (append (list t1 t2 t3 t4) mod-lst))))

(defun count-tot-lines (symbol ampprod)
  (apply #'+ (mapcar (lambda (amp)
                       (count symbol (flatten (content-of amp))))
                     ampprod)))

(defun count-hole-lines (ampprod)
  (apply #'+ (mapcar (lambda (amp)
                       (length (holes-of amp)))
                     ampprod)))

(defun count-loops (ampprod))

; return nil if not replaced
(defun label-amp-hole (symb idx amp)
  (let ((new-ops (replace-once* (lambda (op)
                                  (if (eql op symb) idx))
                                (holes-of amp))))
    (if new-ops
        (make-amp new-ops (particles-of amp)))))
(defun label-amp-particle (symb idx amp)
  (let ((new-ops (replace-once* (lambda (op)
                                  (if (eql op symb) idx))
                                (particles-of amp))))
    (if new-ops
        (make-amp (holes-of amp) new-ops))))
(defun label-contract-line-once (label-func idx ampprod)
  (replace-once (lambda (amp) (funcall label-func 'k idx amp))
                ampprod))

;;; index-the-lines
;;; return a product, sum over the indices which appear twice
(defun label-lines (ampprod)
  (let ((reg 0))
    (flet ((label-iter (op ts label-func)
             (let* ((idx (incf reg))
                    (h-ops (car ts))
                    (new-h (funcall label-func op idx h-ops)))
               (if (eql op '-)
                    (cons h-ops (label-contract-line-once label-func idx (cdr ts)))
                    (cons h-ops (cdr ts))))))
      (reduce (lambda (ts op) (label-iter op ts #'label-contract-hole-once))
              (holes-of (car ampprod))
              :initial-value ampprod))

(defun gen-amp-diagrams (n-lst)
  (flet ((ext-lines (op ampprod)
           (- (count-tot-lines '+ (cons op ampprod)))
              (count-tot-lines '- (list op))))
    (let ((amps (gen-amps-list n-lst))
          (dex-lines (mapcar (lambda (x) (+ x x)) n-lst)))
      (mapcan (lambda (ampprod)
                (mapcan (lambda (h2e-op)
                          (if (member (ext-lines h2e-op ampprod) dex-lines)
                              (contract-h2e-ampprod h2e-op ampprod)))
                        *h2e-ops*))
              amps))))


;;;; vim: ft=lisp
