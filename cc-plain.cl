;;; -*- Mode: Common Lisp; -*-
;;;
;;; Description: based on Goldstone diagram

(load "utility.cl")

(defun attach-tag (tag x)
  (cons tag x))
(defmacro tag-of (x)
  `(car ,x))
(defmacro content-of (x)
  `(cdr ,x))
;;; commutator : p q^\dagger + q^\dagger p = \delta_{pq}

;;; type:
;;; he+   ->-|
;;; pe+   -<-|
;;; he-   |->-
;;; pe-   |-<-
;;; hi+   contracted ->-|
;;; pi+   contracted -<-|
;;; hi-   contracted |->-
;;; pi-   contracted |-<-
;;;
;;; (t_{ijk}^{abc} i a j b k c)

(defun op-copy (n op)
  (loop repeat n
       collect op))

; amp == (t op1 op2 ...)
; n-excitation amplitude without linking
(defun make-new-amp (n)
  (attach-tag 't (flatten (op-copy n '(he+ pe+)))))
; t1 = (t he+ pe+)
; t2 = (t he+ pe+ he+ pe+)

;;; contraction of one hole and an amplitude
(defun contract-hole-amp (amp)
  (let ((ops (replace-once* (lambda (op) (if (eql op 'he+) 'hi+))
                            (content-of amp))))
    (if ops
        (attach-tag 't ops))))
(defun contract-particle-amp (amp)
  (let ((ops (replace-once* (lambda (op) (if (eql op 'pe+) 'pi+))
                            (content-of amp))))
    (if ops
        (attach-tag 't ops))))

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

(defun contract-ops-ampprods (ops ampprod-lst)
  (flet ((hole-func (ampprod-lst op)
           (if (eql op 'he-)
               (let ((res (mapcan #'contract-hole-ampprod ampprod-lst)))
                 (if res res ampprod-lst))
               ampprod-lst))
         (particle-func (ampprod-lst op)
           (if (eql op 'pe-)
               (let ((res (mapcan #'contract-particle-ampprod ampprod-lst)))
                 (if res res ampprod-lst))
               ampprod-lst)))
    (reduce #'particle-func ops ; first reduce all holes to avoid sorting
            :initial-value (reduce #'hole-func ops :initial-value ampprod-lst))))

(defun id-amp (amp &optional (oplist '(hi- pi- hi+ pi+ he- pe- he+ pe+)))
  (let ((shift (length oplist)))
    (reduce (lambda (v x) (+ x (* shift v)))
            (mapcar (lambda (x) (position x oplist)) (content-of amp))
            :initial-value 0)))
(defun sort-ampprod (ampprod)
  (sort ampprod
        (lambda (x y)
          (let ((len-x (length x))
                (len-y (length y)))
            (if (eq len-x len-y)
                (< (id-amp x) (id-amp y))
                (< len-x len-y))))))

;;; result list of contract-ops-ampprods may need to be sorted
(defun contract-ops-ampprods-uniq (ops ampprod-lst)
  (remove-duplicates (contract-ops-ampprods ops ampprod-lst)
                     :test #'equal))

(defun connected-amp? (amp)
  (some (lambda (op) (member op '(hi- pi- hi+ pi+)))
        (content-of amp)))
(defun connected-ampprod? (ampprod)
  (every #'connected-amp? ampprod))

(defun ext-line-to-in-line (op)
  (cond ((eql op 'he-) 'hi-)
        ((eql op 'pe-) 'pi-)
        ((eql op 'he+) 'hi+)
        ((eql op 'pe+) 'pi+)
        (t op)))
(defun contract-h2e-ampprod (h2e ampprod)
  (let ((ctr-h2e (maptree #'ext-line-to-in-line h2e))
        (ctr-amps (contract-ops-ampprods-uniq (content-of h2e)
                                              (list ampprod))))
    (mapcar (lambda (amps) (cons ctr-h2e amps))
            (remove-if-not #'connected-ampprod? ctr-amps))))

(defparameter *h2e-ops*
  '((h he- he+ he+ pe+)
    (h pe- pe+ he+ pe+)
    (h he- he+ he- he+)
    (h pe- pe+ pe- pe+)
    (h he- he+ pe- pe+)
    (h he- pe- he- he+)
    (h he- pe- pe- pe+)
    (h he- pe- he- pe-)))    

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

(defun count-tot-lines (pred ampprod)
  (apply #'+ (mapcar (lambda (amp) (count-if pred (content-of amp)))
                     ampprod)))

(defun count-hole-lines (ampprod)
  (count-tot-lines (lambda (x) (member x '(he- he+ hi- hi+)))
                   ampprod))

; return nil if not replaced
(defun label-amp (test idx amp)
  (let ((ops (replace-once* (lambda (op) (if (funcall test op) idx))
                            (content-of amp))))
    (if ops
        (attach-tag (tag-of amp) ops))))
(defun label-amp-contract-hole (idx amp)
  (label-amp (lambda (op) (or (eql op 'hi-) (eql op 'hi+))) idx amp))
(defun label-amp-contract-particle (idx amp)
  (label-amp (lambda (op) (or (eql op 'pi-) (eql op 'pi+))) idx amp))

;;; if not contracted, return the input ampprod
(defun label-line-once (label-func idx ampprod)
  (replace-once (lambda (amp) (funcall label-func idx amp))
                ampprod))

;;; index-the-lines
;;; return a product, sum over the indices which appear twice
(defun label-lines (ampprod)
  (let ((reg 0))
    (flet ((label-ctr (ts op)
             (cond ((eql op 'hi-)
                    (let ((idx (incf reg)))
                      (cons (label-amp-contract-hole idx (car ts))
                            (label-line-once #'label-amp-contract-hole idx (cdr ts)))))
                   ((eql op 'pi-)
                    (let ((idx (incf reg)))
                      (cons (label-amp-contract-particle idx (car ts))
                            (label-line-once #'label-amp-contract-particle idx (cdr ts)))))
                   (t ts)))
             (label-ext (amp)
               (let ((pairs (mapcar (lambda (op)
                                      (if (or (eql op 'he+) (eql op 'pe+))
                                          (incf reg)
                                          op))
                                    (content-of amp))))
                 (attach-tag (tag-of amp) pairs))))
      (let ((ctr (reduce #'label-ctr (content-of (car ampprod))
                           :initial-value ampprod)))
        (mapcar #'label-ext ctr)))))

; todo
(defun count-loops (ampprod))

;(defun gen-amp-diagrams (n-lst)
;  (flet ((ext-lines (op ampprod)
;           (- (count-tot-lines '+ (cons op ampprod)))
;              (count-tot-lines '- (list op))))
;    (let ((amps (gen-amps-list n-lst))
;          (dex-lines (mapcar (lambda (x) (+ x x)) n-lst)))
;      (mapcan (lambda (ampprod)
;                (mapcan (lambda (h2e-op)
;                          (if (member (ext-lines h2e-op ampprod) dex-lines)
;                              (contract-h2e-ampprod h2e-op ampprod)))
;                        *h2e-ops*))
;              amps))))
